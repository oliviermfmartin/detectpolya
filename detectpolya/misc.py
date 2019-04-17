#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
import re
from detectpolya.internals import *

def estimateProbabilityNucleotide(seq, qual = None, nuc = "A"):

	"""
	Estimate probability of a specific nucleotide to be at a
	position of a read sequence given call errors specified by
	quality string.
	Returns zeros and ones if quality string is absent.

	Args:
		seq (str): Read sequence string
		qual (str): Quality string. This has to be in the same order as `seq`.
			This is not always the case in BAM file where the sequence can be 
			reversed complemented to be aligned to reference. Quality string 
			needs to be reversed in this case.
		nuc (str): nucleotide of intestest. By default, this is "A" so we can
			look for poly-adelynation.
	Returns:
		list
			List of floats between zero and one giving the probability that `nuc`
			is called at positions given by `seq`.

	Examples:
		>>> p = detectpolya.estimateProbabilityNucleotide("CGTTAAATA", "BBBCFFF!B")
		>>> print(p)
		[0.00016706241120909083,
		 0.00016706241120909083,
		 0.00016706241120909083,
		 0.0001327023901844991,
		 0.9998004737685031,
		 0.9998004737685031,
		 0.9998004737685031,
		 0.3333333333333333,
		 0.9994988127663728]
	"""

	# make sure we only use upper cases
	seq = seq.upper()
	nuc = nuc.upper()

	# no quality string is given, return ones if match, zero if other nucleotide
	if qual == None:
		padn = [float(x == nuc) for x in seq]

	else:
		assert len(seq) == len(qual), "Sequence length does not match quality string length"

		# compute base calling error 
		perr = [10 ** -(float(ord(x) - 33) / 10.0) for x in qual]

		# compute probability that base is an adenosine
		padn = []
		for i in range(len(seq)):
			if seq[i] == nuc:
				padn.append(1 - perr[i])
			elif seq[i] not in ["A", "C", "G", "T"]:
				padn.append(0)
			else:
				padn.append(perr[i] * 1.0/3.0)

	return padn

def removeMatches(seq, cigar, remove_five_prime = False, remove_three_prime = False):

	"""
	Replaces nucleotide matching reference by equal sign in read sequence.
	This allows these nucleotides to be ignored by detection algorithms
	while still outputting correct position in read.

	Args:
		seq (str): Read sequence string
		cigar (str): CIGAR string
		remove_five_prime (bool): Remove 5' soft clipped nucleotides
		remove_three_prime (bool): Remove 3' soft clipped nucleotides

	Returns:
		str
			Sequence with nucleotide matching reference replaced by equal sign.

	Examples:
		>>> seq = detectpolya.removeMatches("ACTG", "3M1S")
		>>> print(seq)
		"===G"

		>>> seq = detectpolya.removeMatches("ACTG", "3M1S", True)
		>>> print(seq)
		"A==="
	"""

	def _format_cigar_(cigar):
		result = ""
		n = ""
		operations = None
		for x in cigar:
			if x.isdigit():
				n += x
			else:
				if not x in ["N", "D", "H"]:
					result += x * int(n)
				n = ""

		return result

	# format cigar into string of same length as sequence
	cigar2 = _format_cigar_(cigar)
	assert len(seq) == len(cigar2), "Sequence length does not match CIGAR string length"

	# remove one or two ends (replace S by M in cigar string if at one end)
	if remove_five_prime:
		first_match = re.search("M", cigar2)
		if first_match:
			cigar2 = ''.join([cigar2[i] if i > first_match.start() else "M" for i in range(len(seq))])

	if remove_three_prime:
		first_match = re.search("M", cigar2[::-1])
		if first_match:
			cigar2 = ''.join([cigar2[i] if i < (len(cigar2) - first_match.start()) else "M" for i in range(len(seq))])

	# replace matches by equal sign
	seq2 = ''.join([seq[i] if cigar2[i] == "S" else "=" for i in range(len(seq))])
	return seq2
