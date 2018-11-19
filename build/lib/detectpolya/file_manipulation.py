#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import collections
import HTSeq
import numpy as np
import progressbar
import sys
import detectpolya
from detectpolya.internals import *

def _retrieveFeaturesFromGTF_(gtf_filename, nlines = None):

	"""
	Parses GTF file and outputs exon and transcript features.
	Exon features are used for gene counts and is a HTSeq.GenomicArrayOfSets
	Transcript features are used for annotations and is a dictionnary with the
	key being the gene id and the value a list of lists [[start, end, length]].
	"""

	# setup progress bar
	widgets =  [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()] 
	if nlines == None:
		nlines = sum(1 for line in open(gtf_filename)) # get number of rows in file
	# bar = progressbar.ProgressBar(widgets = widgets, maxval = nlines) # init progress bar
	# bar.start()

	# init array of features
	exon_features = HTSeq.GenomicArrayOfSets("auto", stranded=True) 
	transcript_features = collections.defaultdict(list)
	
	# iterate through gtf file
	i = 0 # line index
	for feature in HTSeq.GFF_Reader(gtf_filename):

		# update progress bar
		i += 1
		if nlines: 
			if i > nlines:
				break
		# bar.update(i)

		# get information from transcript
		if feature.type == "transcript":
			transcript_features[feature.attr["gene_id"]] += [[feature.iv.start, feature.iv.end, feature.iv.length]]

		# get information from exon
		if feature.type == "exon": 
			exon_features[feature.iv] += feature.attr["gene_id"]

	# bar.finish()

	# format information for transcripts
	transcript_features = {k: np.array(v) for k, v in transcript_features.items()} # convert to array
	transcript_features = {k: [min(v[:,0]), max(v[:,1]), np.mean(v[:,2])] for k, v in transcript_features.items()}

	return transcript_features, exon_features

def analyseFile(filename,
	filetype = "fa",
	paired_ends = False,
	gtf_filename = None,
	polya_min_len = 5,
	polya_max_prop_non_a = 0.2,
	polya_seed_len = 3,
	polya_method = "seed",
	primer_seq = None,
	primer_min_len = 10,
	primer_max_dist = 1,
	nlines = None,
	gtf_nlines = None):

	"""
	Detects poly-adenylation in a file containing reads. Supports FASTA, FASTQ, SAM and BAM

    Args:
    	filename (str): Path of file to be analyzed.
	    qual (str): Type of file. Can be "fa" for fasta, "fq" for fastq, "sam" or "bam".
	    gtf_filename (str): GTF filename.
		polya_method (str): Detection algorithm can be `seed` or `window`.
		polya_min_len (int): Minimum length of a poly-adelynated tail.
		polya_max_prop_non_a (float): Maximum proportion of non-adenosines a poly-adelynated tail may contain.
		polya_seed_len (int): Length of seed for seed algorithm. 
	    primer_seq (str): Primer sequence to be detected.
		primer_min_len (int): Minimum length of primer.
		primer_max_l_dist (float): Maximum Levenshtein distance; mismatched nucleotides between reference and primer.
		nlines (int): Maximum number of lines to be read.
		gtf_nlines (int): Maximum number of lines to be read in GTF file.

    Returns:

	"""

	# import pdb

	if filetype == "fa":
		fasta = True
		paired_ends = False
 		def reader(x): 
 			return SeqIO.parse(x, "fasta")

	elif filetype == "fq":
		fasta = True
		paired_ends = False
 		def reader(x): 
 			return SeqIO.parse(x, "fastq")

	elif filetype == "sam":
		fasta = False
		if not paired_ends:
			reader = HTSeq.SAM_Reader
		if paired_ends:
			def reader(x):
				return HTSeq.pair_SAM_alignments(HTSeq.SAM_Reader(x), bundle=True)

	elif filetype == "bam":
		fasta = False
		if not paired_ends:
			reader = HTSeq.BAM_Reader
		if paired_ends:
			def reader(x):
				return HTSeq.pair_SAM_alignments(HTSeq.BAM_Reader(x), bundle=True)
	else:
		raise NotImplementedError("Choose fasta, fastq, sam or bam")


	# parse GTF file
	if gtf_filename != None:
		transcript_features, exon_features = _retrieveFeaturesFromGTF_(gtf_filename, nlines = gtf_nlines)
	else:
		transcript_features = None
		exon_features = None

	# inialize counters and results dictionnary
	total_counts = collections.Counter()
	results      = collections.defaultdict(list)

	# iterate reads in bam files
	if not nlines:
		nlines = sum(1 for line in reader(filename)) # get number of rows in file

	# progress bar
	widgets =  [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()] 
	bar = progressbar.ProgressBar(widgets = widgets, maxval = nlines) 
	bar.start()

	# if fasta file or no exon feature is provided, then gene id is always unmapped
	if filetype == "fasta" or exon_features == None:
		gene_id = "_unmapped"

	# if not paired end, always ignore second
	first_ignore = False
	if not paired_ends:
		second_ignore = True
		second_read = None
	else:
		second_ignore = False

	if primer_seq == None:
		first_primer = None
		second_primer = None

	# iterate through file
	i = 0
	min_len = min(polya_min_len, primer_min_len)
	for bundles in reader(filename):

		# update progress bar
		i += 1
		if nlines: 
			if i > nlines:
				break
		bar.update(i)

		if paired_ends:

			if len(bundles) == 0:
				continue

			bundle = bundles[0] # only use first mapping

			# extract pair
			first_read, second_read = bundle

			if first_read == None or second_read == None: 
				print "ERROR: missing alignment pair"
				break

		else:
			first_read = bundles

		# get info from reads
		if fasta:
			first_seqinfo  = detectpolya.getSeqInfoSeqIO(first_read)

		else:
			first_seqinfo  = detectpolya.getSeqInfoHTSeq(first_read)

			if paired_ends:
				second_seqinfo = detectpolya.getSeqInfoHTSeq(second_read)

		if exon_features:
			gene_id = getGeneID(first_read, second_read, exon_features)

		# add read to count 
		total_counts[gene_id] += 1


		# check if enough clipped nucleotides for there to be a match
		# if not, ignore
		if filetype == "sam" or filetype == "bam":
			first_ignore = len(first_seqinfo["cigar_operations"]) - first_seqinfo["cigar_operations"].count("M") < min_len
			if paired_ends:
				second_ignore = len(second_seqinfo["cigar_operations"]) - second_seqinfo["cigar_operations"].count("M") < min_len

		# remove nucleotides that correspond to matches
		if fasta:
			seq11 = first_seqinfo["seq"]
		else:
			if not first_ignore:
				seq11 = first_seqinfo["clipped_seq"]
			if not second_ignore:
				seq21 = second_seqinfo["clipped_seq"]

		# detect poly-adenlynation (we only look at the 3')
		if not first_ignore:
			first_polya  = detectpolya.detectPolyA(seq11, 
				qual = first_seqinfo["qual"], \
				min_len = polya_min_len, \
				max_prop_non_a = polya_max_prop_non_a, \
				seed_len = polya_seed_len, \
				method = polya_method)

		if not second_ignore:
			second_polya = detectpolya.detectPolyA(seq21, 
				qual = second_seqinfo["qual"], \
				min_len = polya_min_len, \
				max_prop_non_a = polya_max_prop_non_a, \
				seed_len = polya_seed_len, \
				method = polya_method)

		# detect primer in sequence (primer is already reversed complement in function)
		if primer_seq != None:
			if not first_ignore:
				pdb.set_trace()
				first_primer = detectpolya.detectSubSequence(seq11, primer_seq, min_len = primer_min_len, max_l_dist = primer_max_dist)
			if not second_ignore:
				second_primer = detectpolya.detectSubSequence(seq21, primer_seq, min_len = primer_min_len, max_l_dist = primer_max_dist)

		# format results
		r = []
		if not first_ignore:
			r += detectpolya.formatResults(first_polya, first_primer, first_seqinfo, gene_id, transcript_features)
		if not second_ignore:
			r += detectpolya.formatResults(second_polya, second_primer, second_seqinfo, gene_id, transcript_features)
		results[gene_id] += r

	bar.finish()

	# add total count information to results
	for gene_id in total_counts.keys():
		for i in xrange(len(results[gene_id])):
			results[gene_id][i]["read_count"] = total_counts[gene_id]

	return results

def writeResults(results, outf):

	# write header
	print >> outf, "gene_id", "transcript_start", "transcript_end", "transcript_length", \
	"read_name", \
	"chrom", "read_start", "read_end", "read_length",\
	"read_seq", "read_clipped_seq", "read_qual", "read_cigar", "read_reversed_complemented", "read_count", \
	"polya_start_in_genome", "polya_end_in_genome", \
	"polya_start_in_read", "polya_end_in_read", \
	"polya_length", "polya_score", \
	"primer_start_in_genome", "primer_end_in_genome", \
	"primer_start_in_read", "primer_end_in_read", \
	"primer_length", "primer_score"

	# write entries for each gene
	for gene_id in sorted(results.keys()):

		# for a specific gene, retrieve information
		entry = results[gene_id]

		# print one line for every read
		for p in entry:
			print >> outf, gene_id, p["transcript_start"], p["transcript_end"], p["transcript_length"], \
			p["read_name"], \
			p["read_chrom"], p["read_start"], p["read_end"], p["read_length"], \
			p["read_seq"], p["read_clipped_seq"], p["read_qual"], p["read_cigar"], p["read_reversed_complemented"], p["read_count"],\
			p["polya_start_in_genome"], p["polya_end_in_genome"], \
			p["polya_start_in_read"], p["polya_end_in_read"], \
			p["polya_length"], p["polya_score"], \
			p["primer_start_in_genome"], p["primer_end_in_genome"], \
			p["primer_start_in_read"], p["primer_end_in_read"], \
			p["primer_length"], p["primer_score"]

	return 0
