#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
import detectpolya
from collections import namedtuple

# match objects
Match = namedtuple('Match', ['start', 'end', 'score', 'strand'])
WindowMatch = namedtuple('WindowMatch', ['start', 'length', 'score'])

def revComp(seq):
	"""Returns reverse complement of sequence."""
	return str(Seq(seq).reverse_complement())

def comp(seq):
	"""Returns complement of sequence."""
	return str(Seq(seq).complement())

def getGeneID(first_read, second_read, exon_features):

	"""
	Retrieve of identification of gene to which read maps to.
	"""
	# if first_read.aligned and not second_read.aligned:
	# 	import pdb; pdb.set_trace()

	gene_ids = set()
	if first_read.aligned:
		for iv, val in exon_features[first_read.iv].steps():
			gene_ids |= val

	if second_read:
		if second_read.aligned:
			for iv, val in exon_features[second_read.iv].steps():
				gene_ids |= val

	if len(gene_ids) == 1:
		gene_id = list(gene_ids)[0]
	elif len(gene_ids) == 0:
		gene_id = "_no_feature"
	else:
		gene_id = "_ambiguous"

	return gene_id

def getSeqInfoHTSeq(read):

	"""
	Retrieves information from HTSeq object alignement and returns it as
	a dictionnary. This into account the fact that sequences in BAM file 
	may be reversed complemented to correspond to the reference sequence strand.
	"""

	def _makeCigarString_(cigar):
		return ''.join((''.join((str(i.size), i.type)) for i in cigar))

	def _makeCigarOperations_(cigar):
		return ''.join((i.size * i.type for i in cigar))

	# to be able to align to genome, fasta sequences are reversed complemented
	# read.read.seq is the original fasta sequence
	flag = "{0:b}".format(read.flag + 4096)
	reversed_complemented = bool(int(flag[-5]))

	if reversed_complemented:
		seq = revComp(read.read.seq.decode('UTF-8')) # aligned to reference, sequence in the BAM file
		qual = read.read._qualstr.decode('UTF-8')[::-1] # we reverse quality string so it matches BAM sequence
	else: 
		seq = read.read.seq.decode('UTF-8') # aligned to reference, sequence in the BAM and FASTA file
		qual = read.read._qualstr.decode('UTF-8') # quality string already matches

	# we only keep clipped nucleotides, matches are not of immediate interest
	if read.aligned:
		chrom = read.iv.chrom
		start = read.iv.start
		end = read.iv.end
		cigar_string = _makeCigarString_(read.cigar)
		cigar_operations = _makeCigarOperations_(read.cigar)
		clipped_seq = detectpolya.removeMatches(seq, cigar = cigar_string)
	else:
		chrom = None
		start = None
		end = None
		cigar_string = str(len(seq)) + "S"
		cigar_operations = "S" * len(seq)
		clipped_seq = seq

	return {"name":   read.read.name,
			"chrom":  chrom, 
			"start":  start, 
			"end":    end,
			"length": len(seq), # read.iv.length,
			"seq":    seq, 
			"clipped_seq": clipped_seq, 
			"reversed_complemented": reversed_complemented,
			"qual": qual,
			"cigar_operations": cigar_operations,
			"cigar_string": cigar_string,
			"aligned": read.aligned}

def getSeqInfoSeqIO(read, filetype):

	"""
	Retrieves information from Bio.SeqIO object alignement and returns it as
	a dictionnary.
	"""

	# def _mask5Prime_(seq): 
	# 	return "".join([seq[i] if i >= float(len(seq))/1.5 else "=" for i in range(len(seq))])

	seqinfo = {"name":   read.id,
			   "seq":    str(read.seq), 
			   # "clipped_3p_seq": _mask5Prime_(read.seq),
			   "length": str(len(read.seq))}

	if filetype == "fq":
		seqinfo["qual"] =''.join(chr(x + 33) for x in read.letter_annotations["phred_quality"])

	return seqinfo

def formatResults(polya, primer, seqinfo, gene_id, transcript_features):

	'''
	Formats results into a dictionnary
	'''

	if polya == None and primer == None:
		return []

	# handle read information
	read_name                  = seqinfo.get("name")
	read_mate                  = seqinfo.get("mate")
	read_chrom                 = seqinfo.get("chrom")
	read_start                 = seqinfo.get("start")
	read_end                   = seqinfo.get("end")
	read_length                = seqinfo.get("length")
	read_seq                   = seqinfo.get("seq")
	read_clipped_seq           = seqinfo.get("clipped_seq")
	read_clipped_3p_seq        = seqinfo.get("clipped_3p_seq")
	read_qual                  = seqinfo.get("qual")
	read_cigar                 = seqinfo.get("cigar_string")
	read_reversed_complemented = seqinfo.get("reversed_complemented")

	# handle transcript information
	if transcript_features == None:
		transcript_start  = None
		transcript_end    = None
		transcript_length = None
		transcript_strand = None

	else:
		tf = transcript_features.get(gene_id)
		if tf == None:
			transcript_start  = None
			transcript_end    = None
			transcript_length = None
			transcript_strand = None

		else:
			transcript_start  = str(tf[0]) # transcript start position in reference genome
			transcript_end    = str(tf[1]) # transcript end position
			transcript_length = str(tf[2]) # transcript mean length
			transcript_strand = "/".join(tf[3])

	# handle poly-a information
	if polya == None:
		polya_start_in_genome = None
		polya_end_in_genome   = None
		polya_start_in_read   = None
		polya_end_in_read     = None
		polya_length          = None
		polya_score           = None
		polya_strand          = None

	else:
		if seqinfo.get("start"):
			if polya.strand == "+":
				polya_start_in_genome = str(seqinfo.get("start") + polya.start)
				polya_end_in_genome   = str(seqinfo.get("start") + polya.end - 1)
			elif polya.strand == "-":
				polya_start_in_genome = str(seqinfo.get("end")   - polya.end - 1)
				polya_end_in_genome   = str(seqinfo.get("end")   - polya.start)
		else:
			polya_start_in_genome = None
			polya_end_in_genome   = None

		polya_start_in_read = str(polya.start + 1)
		polya_end_in_read   = str(polya.end)
		polya_length        = str(polya.end - polya.start)
		polya_score         = str(polya.score)
		polya_strand        = str(polya.strand) 

	# handle primer information
	if primer == None:
		primer_start_in_genome = None
		primer_end_in_genome   = None
		primer_start_in_read   = None
		primer_end_in_read     = None
		primer_length          = None
		primer_score           = None
		primer_strand          = None

	else:
		if seqinfo.get("start"):
			if primer.strand == "+":
				primer_start_in_genome = str(seqinfo.get("start") + primer.start)
				primer_end_in_genome   = str(seqinfo.get("start") + primer.end - 1)
			elif primer.strand == "-":
				primer_start_in_genome = str(seqinfo.get("end")   - primer.end - 1)
				primer_end_in_genome   = str(seqinfo.get("end")   - primer.start)
		else:
			primer_start_in_genome = None
			primer_end_in_genome   = None

		primer_start_in_read = str(primer.start + 1)
		primer_end_in_read   = str(primer.end)
		primer_length        = str(primer.end - primer.start)
		primer_score         = str(primer.score)
		primer_strand        = str(primer.strand) 

	return [{"gene_id": gene_id,
			"transcript_start": transcript_start,
			"transcript_end": transcript_end,
			"transcript_length": transcript_length,
			"transcript_strand": transcript_strand,
			"read_name":  read_name,
			"read_mate":  read_mate,
			"read_chrom": read_chrom,
			"read_start": read_start,
			"read_end": read_end,
			"read_length": read_length,
			"read_seq": read_seq,
			"read_clipped_seq": read_clipped_seq,
			"read_clipped_3p_seq": read_clipped_3p_seq,
			"read_qual": read_qual,
			"read_cigar": read_cigar,
			"read_reversed_complemented": read_reversed_complemented,
			"polya_start_in_genome": polya_start_in_genome,
			"polya_end_in_genome": polya_end_in_genome,
			"polya_start_in_read": polya_start_in_read,
			"polya_end_in_read": polya_end_in_read,
			"polya_length": polya_length,
			"polya_score": polya_score,
			"polya_strand": polya_strand,
			"primer_start_in_genome": primer_start_in_genome,
			"primer_end_in_genome": primer_end_in_genome,
			"primer_start_in_read": primer_start_in_read,
			"primer_end_in_read": primer_end_in_read,
			"primer_length": primer_length,
			"primer_score": primer_score,
			"primer_strand": primer_strand}]
