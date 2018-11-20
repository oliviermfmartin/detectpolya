#!/usr/bin/python
# -*- coding: utf-8 -*-

import fuzzysearch
import numpy as np
import re
from detectpolya.internals import *
import detectpolya

def detectPolyA(seq, qual = None, method = "seed", min_len = 5, 
	max_prop_non_a = 0.2, seed_len = 4):

	"""
	Detects poly-adenylation in a read sequence. 

	If quality string is specified (`qual`) computes the probability that a 
	nucleotide is an adenosines and uses the expected number of adenosines in 
	a subsequence instead of the exact count of called adenosines. Only 
	returns longest match if it is more than a certain length.

	* `seed`: Heuristic similar to the BLAST algorithm. First matches a seed 
		subsequence of consecutive adenosines, then tries to extend it until until 
		a the subsequence reaches the max proportion of non-adenosines nucleotides 
		(`max_prop_non_a`);

	* `window`: Exact method that recursively computes the proportion of 
		adenosines for all possible subsequences. Returns longest match that is 
		not below the max proportion of non-adenosines nucleotides if any.

	Notes:
		The quality string is not always in the same order as the sequence. 
		This is true in the BAM file where the sequence can be reversed 
		complemented to be aligned to reference. Quality string needs to be 
		reversed in this case.

	Args:
		seq (str): Read sequence string.
		qual (str): Quality string.
		method (str): Detection algorithm can be `seed` or `window`.
		min_len (int): Minimum length of a poly-adelynated tail.
		max_prop_non_a (float): Maximum proportion of non-adenosines a 
			poly-adelynated tail may contain.
		seed_len (int): Length of seed for seed algorithm. 

	Returns:
		collection.namedtuple
			If a match is found, return a named tupple: start, end and score. 
			The score corresponds to the number of (expected) matched adenosines.

	Examples:
		>>> polya = detectpolya.detectPolyA("ACTGGTAAAAAA")
		>>> print(polya)
		Match(start=5, end=12, score=6.0)

		>>> polya = detectpolya.detectPolyA("ACTGGTGTACAT")
		>>> print(polya)
		None
	"""

	if method == "seed":
		return _detectPolyASeed_(seq = seq, qual = qual, min_len = min_len, max_prop_non_a = max_prop_non_a, seed_len = seed_len)

	elif method == "window":
		return _detectPolyAWindow_(seq = seq, qual = qual, min_len = min_len, max_prop_non_a = max_prop_non_a)

	else:
		raise NotImplementedError("Choose seed or window method")


def _detectPolyAWindow_(seq, qual, min_len, max_prop_non_a):

	"""
	Subroutine of detectPolyA. Implements window method.
	See detectPolyA documentation for more information.
	"""

	# make sure nucleotides are upper case
	seq = seq.upper()

	# shape
	shape = (max_polya+1, len(seq))

	# compute A count/expected value for window size of one
	count = np.zeros(shape, dtype=float)
	count[1,:] = detectpolya.estimateProbabilityNucleotide(seq = seq, qual = qual, nuc = "A")
	
	# compute A count for other window sizes
	for i in xrange(2, max_polya+1): # window size
		for j in xrange(0, len(seq)): # sequence position
			if i + j > len(seq): # these positions are not computable
				continue
			count[i, j] = count[i-1, j] + count[i-1, j+1] - count[i-2, j+1]

	# proportion of adenosines
	prop = np.zeros(shape, dtype=float)
	for i in xrange(2, shape[0]):
		prop[i,:] = [count[i, j] / float(i) for j in xrange(shape[1])]

	# match matrix
	match = np.zeros(shape, dtype=int)
	for i in xrange(min_len, shape[0]):
		match[i,:] = [int(i - a <= i * max_prop_non_a) for a in count[i,:]]

	# select left-low match (longest and furthest)
	if 1 in match:
		start = np.argmax(np.amax(match[:,], 0))
		len_seq = max_polya - np.argmax(match[:,start][::-1])
		end = start + len_seq
		expade = count[len_seq, start]
		return Match(start, end, expade)
	else:
		return None

def _detectPolyASeed_(seq, qual, min_len, max_prop_non_a, seed_len):

	"""
	Subroutine of detectPolyA. Implements seed method.
	See detectPolyA documentation for more information.
	"""

	def _matchPolya_(seq, pattern, padn, max_prop_non_a):

		# get pattern matches
		seed = pattern.finditer(seq)
		matches = []
		for match in pattern.finditer(seq):
			expade = sum(padn[match.start():match.end()]) # expected number of adenosines
			matches.append([match.start(), match.end() - 1, expade]) # append to matches (start, end, expected number of adenosines)

		# try to extend matches
		for i in xrange(len(matches)):

			# try to extend to 3'
			start = matches[i][0]
			end   = matches[i][1]
			expade    = matches[i][2]
			newexpade = matches[i][2]
			for j in xrange(end + 1, len(seq), 1):
				newexpade = expade + padn[j]

				# abort if low probability or finds another match
				len_seq = j - start + 1

				if seq[j] == "=" or len_seq - newexpade > len_seq * max_prop_non_a:
					matches[i][1] = j-1
					matches[i][2] = expade
					break
				elif j in [x[0] for x in matches]: # all start positions
					matches[i][1] = j-1
					matches[i][2] = expade
					break
				elif j == len(seq):
					matches[i][1] = j
					matches[i][2] = newexpade

				expade = newexpade

			# try to extend to 5'
			start = matches[i][0]
			end   = matches[i][1]
			expade 	 = matches[i][2]
			newexpade = matches[i][2]

			for j in xrange(start - 1, -1, -1):
				newexpade = expade + padn[j]

				# abort if low probability or finds another match
				len_seq = end - j + 1

				if seq[j] == "=" or len_seq - newexpade > len_seq * max_prop_non_a:
					matches[i][0] = j+1
					matches[i][2] = expade
					break
				elif j in [x[1] for x in matches]: # all end positions
					matches[i][0] = j+1
					matches[i][2] = expade
					break
				elif j == 0:
					matches[i][0] = j
					matches[i][2] = newexpade

				expade = newexpade


		matches = [[x[0], x[1]+1, x[2]] for x in matches]
		return matches

	def _mergePolyaMatches_(x):

		def __mergePolyaMatches__(x):

			if x[0][1] == x[1][0]: # we should merge
				return [[x[0][0], x[1][1], x[0][2] + x[1][2]]]
			else: # don't merge
				return x

		merged = __mergePolyaMatches__(x[0:2])

		for i in xrange(2, len(x)):
			merged = merged + [x[i]]
			merged = merged[:-2] + __mergePolyaMatches__(merged[-2:])

		return merged

	# make sure sequence is all upper cases
	seq  = seq.upper()

	# compute probability that base is an adenosine
	padn = detectpolya.estimateProbabilityNucleotide(seq, qual, "A")

	# match poly-a tails
	matches = _matchPolya_(seq, re.compile("A{" + str(seed_len) + ",}"), padn, max_prop_non_a)

	# if no match return nothing
	if len(matches) == 0:
		return None

	# merge hits if they are contingent
	elif len(matches) > 1:
		matches = _mergePolyaMatches_(matches)

	# keep only longest match
	result = matches[np.argmax([m[2] for m in matches])]
	result = Match(result[0], result[1], result[2])

	# return result if long enough
	if result.end - result.start >= min_len:
		return result
	else:
		return None

def detectSubSequence(seq, subseq, min_len = 7, max_l_dist = 1):

	"""
	Uses the `taleinat/fuzzysearch` package to look for a subsequence within a 
	sequence while allowing for mismatches. Looks for subsequence itself and 
	its reverse complement. Only returns longest match if it is  more than a 
	certain length.

	Args:
		seq (str): Read sequence string.
		subseq (str): Subsequence to be detected.
		min_len (int): Minimum length of subsequence match.
		max_l_dist (float): Maximum Levenshtein distance.
		
	Returns:
		collection.namedtuple
			If a match is found, return a named tupple: start, end and score. 
			The score corresponds to the number of matched nucleotides.

	Examples:
		>>> match = detectpolya.detectSubSequence("AAATATAAATACCC", "TATATATA");
		>>> print(match)
		Match(start=3, end=10, score=6)
	"""

	def _returnLongestMatch_(near_matches):
		lengths      = [x.end - x.start for x in near_matches]
		max_index    = lengths.index(max(lengths))
		near_match = near_matches[max_index]
		return Match(near_match.start, near_match.end, near_match.end - near_match.start - near_match.dist)

	# subsequence could be in two different orientations
	subseqs = [subseq, revComp(subseq)]
	revcomp = revComp(seq)

	# try to get near matches
	near_matches = []
	for i in xrange(len(subseqs)):

		p = subseqs[i]

		# actual sequence - match longer and longer sequences
		for length in xrange(min_len, len(subseq)):
			nm = fuzzysearch.find_near_matches(p[0:length], seq, max_l_dist = max_l_dist)
			near_matches += nm

			if len(nm) == 0:
				break

		# reverse complement of sequence
		for length in xrange(min_len, len(subseq)):
			nm = fuzzysearch.find_near_matches(p[0:length], revcomp, max_l_dist = max_l_dist)
			nm = [fuzzysearch.Match(start = len(seq) - x.end, end = len(seq) - x.start, dist = x.dist) for x in nm]
			near_matches += nm

			if len(nm) == 0:
				break

	# only return longest match
	if len(near_matches) == 0:
		return None
	else:
		return _returnLongestMatch_(near_matches)
