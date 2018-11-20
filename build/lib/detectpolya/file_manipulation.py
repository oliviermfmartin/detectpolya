#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import collections
from itertools import izip
import HTSeq
import numpy as np
import progressbar
import sys
import warnings
import detectpolya
from detectpolya.internals import *

def _retrieveFeaturesFromGTF_(gtf_filename, nlines = None, verbose = True):

	"""
	Parses GTF file and outputs exon and transcript features.
	Exon features are used for gene counts and is a HTSeq.GenomicArrayOfSets
	Transcript features are used for annotations and is a dictionnary with the
	key being the gene id and the value a list of lists [[start, end, length]].
	"""

	if verbose:
		print "Reading GTF file"
		# setup progress bar
		widgets =  [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()] 
		if nlines == None:
			nlines = sum(1 for line in open(gtf_filename)) # get number of rows in file
		bar = progressbar.ProgressBar(widgets = widgets, maxval = nlines) # init progress bar
		bar.start()

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
		if verbose:
			bar.update(i)

		# get information from transcript
		if feature.type == "transcript":
			transcript_features[feature.attr["gene_id"]] += [[feature.iv.start, feature.iv.end, feature.iv.length]]

		# get information from exon
		if feature.type == "exon": 
			exon_features[feature.iv] += feature.attr["gene_id"]

	if verbose:
		bar.finish()

	# format information for transcripts
	transcript_features = {k: np.array(v) for k, v in transcript_features.items()} # convert to array
	transcript_features = {k: [min(v[:,0]), max(v[:,1]), np.mean(v[:,2])] for k, v in transcript_features.items()}

	return transcript_features, exon_features

def analyseFile(filename,
	filetype = "fa",
	paired_ends = False,
	gtf_filename = None,
	polya_method = "seed",
	polya_min_len = 5,
	polya_max_prop_non_a = 0.2,
	polya_seed_len = 3,
	primer_seq = None,
	primer_min_len = 10,
	primer_max_dist = 1,
	nlines = None,
	gtf_nlines = None,
	verbose = True):

	"""
	Detects poly-adenylation in a file containing reads. Supports FASTA, FASTQ, SAM and BAM

	Args:
		filename (str): Path of file to be analyzed. For paired ended fasta and fastq, must be a list of filenames.
		qual (str): Type of file. Can be "fa" for fasta, "fq" for fastq, "sam" or "bam".
		paired_ends (bool): Was the experiment paired-ended?
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
		verbose (bool): Print progress

	Returns:
		dict
			Keys are Gene IDs as found in the GTF file, if specified. 
			If no GTF was specified, there is only one key: "_unmapped".
			Values is a list of dictionnaries. 
			Every dictionnary correspond to a read with a poly-A or primer detected.
			The keys are described in README.

	Examples:
		>>> results = detectpolya.analyseFile("./test/paired_ends.sam", "sam")
		>>> print(results)
		defaultdict(list,
			{'_unmapped': [{'gene_id': '_unmapped',
			   'polya_end_in_genome': '15063255',
			   'polya_end_in_read': '72',
			   'polya_length': '21',
			   'polya_score': '16.9969746427',
			   'polya_start_in_genome': '15063184',
			   'polya_start_in_read': '52',
			   'primer_end_in_genome': None,
			   'primer_end_in_read': None,
			   'primer_length': None,
			   'primer_score': None,
			   'primer_start_in_genome': None,
			   'primer_start_in_read': None,
			   'read_chrom': 'I',
			   'read_cigar': '51M25S',
			   'read_clipped_seq': '===================================================AAAAAAAAAAAAAAAAGTACTCTGC',
			   'read_count': 18,
			   'read_end': 15063184,
			   'read_length': 51,
			   'read_mate': None,
			   'read_name': 'D00224L:232:CCB68ANXX:4:1101:14325:85615',
			   'read_qual': '0000<<7<707070FBB<00000000000<0FB<0007BBFFFFFFFB<0FFFFFFFFFFFFFFFFFFFFFBB<<0',
			   'read_reversed_complemented': True,
			   'read_seq': 'AGTTTTTCGTTTCCGGGGGTAGTATGGTTGAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTACTCTGC',
			   'read_start': 15063133,
			   'transcript_end': None,
			   'transcript_length': None,
			   'transcript_start': None}]})
	"""

	if filetype == "fa":
		fasta = True

		if paired_ends:
			assert len(filename) == 2, "Filenames of pair-ended FASTA files must be a list of size 2"
			def reader(x): 
				return izip(SeqIO.parse(x[0], "fasta"), SeqIO.parse(x[1], "fasta"))

		else:
			def reader(x): 
				return SeqIO.parse(x, "fasta")

	elif filetype == "fq":
		fasta = True

		if paired_ends:
			assert len(filename) == 2, "Filenames of pair-ended FASTA files must be a list of size 2"
			def reader(x): 
				return izip(SeqIO.parse(x[0], "fastq"), SeqIO.parse(x[1], "fastq"))

		else:
			def reader(x): 
				return SeqIO.parse(x, "fastq")

	elif filetype == "sam":
		fasta = False
		if not paired_ends:
			def reader(x):
				return HTSeq.bundle_multiple_alignments(HTSeq.SAM_Reader(x))
		if paired_ends:
			def reader(x):
				return HTSeq.pair_SAM_alignments(HTSeq.SAM_Reader(x), bundle=True)

	elif filetype == "bam":
		fasta = False
		if not paired_ends:
			def reader(x):
				return HTSeq.bundle_multiple_alignments(HTSeq.BAM_Reader(x))
		if paired_ends:
			def reader(x):
				return HTSeq.pair_SAM_alignments(HTSeq.BAM_Reader(x), bundle=True)
	else:
		raise NotImplementedError("Choose fasta, fastq, sam or bam")


	# parse GTF file
	if gtf_filename != None:
		if fasta:
			warnings.warn("GTF file is not used for FASTA and FASTQ files")
			transcript_features = None
			exon_features = None
		else:
			if verbose:
				transcript_features, exon_features = _retrieveFeaturesFromGTF_(gtf_filename, nlines = gtf_nlines)
	else:
		transcript_features = None
		exon_features = None

	# inialize counters and results dictionnary
	total_counts = collections.Counter()
	results      = collections.defaultdict(list)

	# if fasta file or no exon feature is provided, then gene id is always unmapped
	if fasta or exon_features == None:
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

	# iterate reads in bam files
	if not nlines:
		nlines = sum(1 for line in reader(filename)) # get number of rows in file

	# set up progress bar
	if verbose:
		print "Iterating file"
		widgets =  [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()] 
		bar = progressbar.ProgressBar(widgets = widgets, maxval = nlines) 
		bar.start()

	# iterate through file
	i = 0
	min_len = min(polya_min_len, primer_min_len)
	for bundles in reader(filename):

		# update progress bar
		i += 1
		if nlines: 
			if i > nlines:
				break
		if verbose:
			bar.update(i)

		# retrieve entry
		if not fasta:	 
			if len(bundles) == 0:
				continue
			# if multiple mapping, select first
			bundles = bundles[0]

		if paired_ends:
			first_read, second_read = bundles # extract pair
		else:
			first_read = bundles

		# get info from reads
		if fasta:
			first_seqinfo  = detectpolya.getSeqInfoSeqIO(first_read, filetype)
			if paired_ends:
				second_seqinfo = detectpolya.getSeqInfoSeqIO(second_read, filetype)

		else:
			first_seqinfo  = detectpolya.getSeqInfoHTSeq(first_read)
			if paired_ends:
				second_seqinfo = detectpolya.getSeqInfoHTSeq(second_read)

		if exon_features:
			gene_id = getGeneID(first_read, second_read, exon_features)

		if paired_ends:
			first_seqinfo["mate"] = "1"
			second_seqinfo["mate"] = "2"

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
			if not first_ignore:
				seq11 = first_seqinfo["seq"]
			if not second_ignore:
				seq21 = second_seqinfo["seq"]
		else:
			if not first_ignore:
				seq11 = first_seqinfo["clipped_seq"]
			if not second_ignore:
				seq21 = second_seqinfo["clipped_seq"]

		# detect poly-adenlynation (we only look at the 3')
		if not first_ignore:
			first_polya  = detectpolya.detectPolyA(seq11, 
				qual = first_seqinfo.get("qual"), \
				min_len = polya_min_len, \
				max_prop_non_a = polya_max_prop_non_a, \
				seed_len = polya_seed_len, \
				method = polya_method)

		if not second_ignore:
			second_polya = detectpolya.detectPolyA(seq21, 
				qual = second_seqinfo.get("qual"), \
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

	if verbose:
		bar.finish()

	# add total count information to results
	for gene_id in total_counts.keys():
		for i in xrange(len(results[gene_id])):
			results[gene_id][i]["read_count"] = total_counts[gene_id]

	return results

def printResults(results, outf = None, header = True):

	"""
	Write results from analyseFile to a file.

	Args:
		results (dict): Output of analyseFile
		outf (file): Filed open for writing. If not specified, results are printed.
		header (bool): Should header be printed?

	Returns:
		None

	Examples:
		>>> results = detectpolya.analyseFile("./test/paired_ends.sam", "sam")
		>>> print detectpolya.printResults(results)
		gene_id,transcript_start,transcript_end,transcript_length,read_name,read_mate,chrom,read_start,read_end,read_length,read_seq,read_clipped_seq,read_qual,read_cigar,read_reversed_complemented,read_count,polya_start_in_genome,polya_end_in_genome,polya_start_in_read,polya_end_in_read,polya_length,polya_score,primer_start_in_genome,primer_end_in_genome,primer_start_in_read,primer_end_in_read,primer_length,primer_score
		_unmapped,NA,NA,NA,D00224L:232:CCB68ANXX:4:1101:14325:85615,NA,I,15063133,15063184,51,AGTTTTTCGTTTCCGGGGGTAGTATGGTTGAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTACTCTGC,===================================================AAAAAAAAAAAAAAAAGTACTCTGC,0000<<7<707070FBB<00000000000<0FB<0007BBFFFFFFFB<0FFFFFFFFFFFFFFFFFFFFFBB<<0,51M25S,True,18,15063184,15063255,52,72,21,16.9969746427,NA,NA,NA,NA,NA,NA
	"""

	def g(d, v):
		r = d.get(v)
		if r == None: r = "NA"
		else: r = str(r)
		return r

	# write header
	if header:
		row = ["gene_id", "transcript_start", "transcript_end", "transcript_length", \
			"read_name", "read_mate", \
			"chrom", "read_start", "read_end", "read_length",\
			"read_seq", "read_clipped_seq", "read_qual", \
			"read_cigar", "read_reversed_complemented", "read_count", \
			"polya_start_in_genome", "polya_end_in_genome", \
			"polya_start_in_read", "polya_end_in_read", \
			"polya_length", "polya_score", \
			"primer_start_in_genome", "primer_end_in_genome", \
			"primer_start_in_read", "primer_end_in_read", \
			"primer_length", "primer_score"]
		row = ','.join(row)

		if outf != None:
			print >> outf, row
		else:
			print row

	# write entries for each gene
	for gene_id in sorted(results.keys()):

		# for a specific gene, retrieve information
		entry = results[gene_id]

		# print one line for every read
		for p in entry:

			row = [gene_id, g(p, "transcript_start"), g(p, "transcript_end"), g(p, "transcript_length"), \
				g(p, "read_name"), g(p, "read_mate"), \
				g(p, "read_chrom"), g(p, "read_start"), g(p, "read_end"), g(p, "read_length"), \
				g(p, "read_seq"), g(p, "read_clipped_seq"), g(p, "read_qual"), \
				g(p, "read_cigar"), g(p, "read_reversed_complemented"), g(p, "read_count"),\
				g(p, "polya_start_in_genome"), g(p, "polya_end_in_genome"), \
				g(p, "polya_start_in_read"), g(p, "polya_end_in_read"), \
				g(p, "polya_length"), g(p, "polya_score"), \
				g(p, "primer_start_in_genome"), g(p, "primer_end_in_genome"), \
				g(p, "primer_start_in_read"), g(p, "primer_end_in_read"), \
				g(p, "primer_length"), g(p, "primer_score")]
			row = ','.join(row)

			if outf != None:
				print >> outf, row
			else:
				print row

	return None
