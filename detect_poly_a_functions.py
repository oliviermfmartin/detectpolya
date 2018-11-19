import sys
import progressbar
import re
import collections
import numpy as np
import HTSeq
import fuzzysearch
import pdb
from Bio.Seq import Seq

def mask5Prime(seq, mask_until = 30):
	return ''.join([seq[i] if i >= mask_until else "=" for i in xrange(len(seq))])

def revComp(seq):
	return str(Seq(seq).reverse_complement())


def comp(seq):
	return str(Seq(seq).complement())

def removeMatches(seq, cigar, reversed_complemented = False):

	'''
	replace call matches by equal sign so they are ignored 
	by poly-a and primer detection algorithms
	'''

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

	cigar2 = _format_cigar_(cigar)
	assert len(seq) == len(cigar2), "sequence length does not cigar length"
	if reversed_complemented:
		cigar2 = cigar2[::-1]

	seq2 = ''.join([seq[i] if cigar2[i] != "M" else "=" for i in xrange(len(seq))])
	return seq2

def _estimateProbAdenosine_(seq, qual = None):

	if qual == None:
		padn = [float(nuc == "A") for nuc in seq]

	else:
		assert len(seq) == len(qual), "sequence length does not match quality string length"

		# compute base calling error 
		perr = [10 ** -(float(ord(x) - 33) / 10.0) for x in qual]

		# compute probability that base is an adenosine
		padn = []
		for i in xrange(len(seq)):
			if seq[i] == "A":
				padn.append(1 - perr[i])
			elif seq[i] == "=":
				padn.append(0)
			else:
				padn.append(perr[i] * 1.0/3.0)

	return padn

def detectPolyA(seq, qual = None, method = "seed", min_len = 5, max_percent_non_adenosines = 0.2, seed_len = 4):
	if method == "seed":
		return detectPolyASeed(seq = seq, qual = qual, min_len = min_len, max_percent_non_adenosines = max_percent_non_adenosines, seed_len = seed_len)
	elif method == "window":
		return detectPolyAWindow(seq = seq, qual = qual, min_len = min_len, max_percent_non_adenosines = max_percent_non_adenosines)	
	else:
		raise NotImplementedError("choose seed or window")


def detectPolyAWindow(seq, qual, min_len, max_percent_non_adenosines):

	max_polya = len(seq)

	# make sure nucleotides are upper case
	seq = seq.upper()

	# shape
	shape = (max_polya+1, len(seq))

	# compute A count/expected value for window size of one
	count = np.zeros(shape, dtype=float)
	count[1,:] = _estimateProbAdenosine_(seq = seq, qual = qual)
	
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
		match[i,:] = [int(i - a <= i * max_percent_non_adenosines) for a in count[i,:]]

	# select left-low match (longest and furthest)
	if 1 in match:
		start = np.argmax(np.amax(match[:,], 0))
		len_seq = max_polya - np.argmax(match[:,start][::-1])
		end = start + len_seq
		expade = count[len_seq, start]
		return [[start, end, expade]]
	else:
		return []

def detectPolyASeed(seq, qual, min_len, max_percent_non_adenosines, seed_len):
	'''
	detects poly-adenlynation for a sequence while taking into account call errors
	first matches a "seed", a sequence of adenosines of length `seed_len`,
	then tries to extend it
	only returns poly-a tails longer than `min_len` nucleotides
	'''

	def _match_polya_(seq, pattern, padn, max_percent_non_adenosines):

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

				if seq[j] == "=" or len_seq - newexpade > len_seq * max_percent_non_adenosines:
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

				if seq[j] == "=" or len_seq - newexpade > len_seq * max_percent_non_adenosines:
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

	def _merge_polya_matches_(x):


		def __merge_polya_matches__(x):

			if x[0][1] == x[1][0]: # we should merge
				return [[x[0][0], x[1][1], x[0][2] + x[1][2]]]
			else: # don't merge
				return x

		merged = __merge_polya_matches__(x[0:2])

		for i in xrange(2, len(x)):
			merged = merged + [x[i]]
			merged = merged[:-2] + __merge_polya_matches__(merged[-2:])

		return merged

	# make sure sequence is all upper cases
	seq  = seq.upper()

	# compute probability that base is an adenosine
	padn = _estimateProbAdenosine_(seq, qual)

	# match poly-a tails
	matches = _match_polya_(seq, re.compile("A{" + str(seed_len) + ",}"), padn, max_percent_non_adenosines)

	# if no match return nothing
	if len(matches) == 0:
		return []

	# merge hits if they are contingent
	elif len(matches) > 1:
		matches = _merge_polya_matches_(matches)

	# keep only longest match
	result = [matches[np.argmax([m[2] for m in matches])]]

	# return result if long enough
	if result[0][1] - result[0][0] >= min_len:
		return result
	else:
		return []

def getSeqInfo(almnt):

	'''
	retrieves information from HTSeq object alignements
	'''

	def _makeCigarString_(cigar):
		return ''.join((''.join((str(i.size), i.type)) for i in almnt.cigar))

	def _makeCigarOperations_(cigar):
		return ''.join((i.size * i.type for i in almnt.cigar))

	# to be able to align to genome, fasta sequences are reversed complemented
	# almnt.read.seq is the original fasta sequence
	# the cigar string refers to this sequence
	# we want the aligned sequence
	# this can be gathered in the bam flag
	reversed_complemented = bool(int("{0:b}".format(almnt.flag)[-5]))

	if reversed_complemented:
		seq = revComp(almnt.read.seq) # aligned to reference, sequence in the BAM file
		qual = almnt.read._qualstr[::-1] # we reverse quality string so it matches BAM sequence
		cigar_string = _makeCigarString_(almnt.cigar[::-1])
		cigar_operations = _makeCigarOperations_(almnt.cigar[::-1])
	else: 
		seq = almnt.read.seq # aligned to reference, sequence in the BAM and FASTA file
		qual = almnt.read._qualstr # quality string already matches
		cigar_string = _makeCigarString_(almnt.cigar)
		cigar_operations = _makeCigarOperations_(almnt.cigar)

	# we only keep clipped nucleotides, matches are not of immediate interest
	clipped_seq = removeMatches(seq, cigar_string)

	return {"name":   almnt.read.name,
			"chrom":  almnt.iv.chrom, 
			"start":  almnt.iv.start, 
			"end":    almnt.iv.end,
			"length": almnt.iv.length,
			"seq":    seq, 
			"clipped_seq":    clipped_seq, 
			"reversed_complemented": reversed_complemented,
			"qual":                  qual,
			"cigar_operations":     cigar_operations,
			"cigar_string":         cigar_string}

def formatResults(polya, primer, seqinfo, identif):

	'''
	formats results into a dictionnary
	'''

	if len(primer) == 0 and len(polya) == 0:
		return []

	if len(primer) == 0:
		primer_start_in_genome = "NA"
		primer_end_in_genome = "NA"
		primer_start_in_read = "NA"
		primer_end_in_read = "NA"
		primer_length = "NA"
		primer_distance = "NA"

	else:
		primer_start_in_genome = str(seqinfo["start"] + primer.start)
		primer_end_in_genome = str(seqinfo["end"] + primer.end - 1)
		primer_start_in_read = str(primer.start + 1)
		primer_end_in_read = str(primer.end)
		primer_length = str(primer.end - primer.start)
		primer_distance = str(primer.dist)

	if len(polya) == 0:
		results = [{"alignment": identif,
					"read_name":  seqinfo["name"],
					"read_chrom": seqinfo["chrom"],
					"read_start": seqinfo["start"],
					"read_end": seqinfo["start"],
					"read_length": seqinfo["length"],
					"read_seq": seqinfo["seq"],
					"read_clipped_seq": seqinfo["clipped_seq"],
					"read_qual": seqinfo["qual"],
					"read_cigar": seqinfo["cigar_string"],
					"read_reversed_complemented": str(seqinfo["reversed_complemented"]),
					"polya_start_in_genome": "NA",
					"polya_end_in_genome": "NA",
					"polya_start_in_read": "NA",
					"polya_end_in_read": "NA",
					"polya_length": "NA",
					"polya_nadenosines": "NA",
					"primer_start_in_genome": primer_start_in_genome,
					"primer_end_in_genome": primer_end_in_genome,
					"primer_start_in_read": primer_start_in_read,
					"primer_end_in_read": primer_end_in_read,
					"primer_length": primer_length,
					"primer_distance": primer_distance}]

	else:

		results = [{"alignment": identif,
					"read_name":  seqinfo["name"],
					"read_chrom": seqinfo["chrom"],
					"read_start": seqinfo["start"],
					"read_end": seqinfo["start"],
					"read_length": seqinfo["length"],
					"read_seq": seqinfo["seq"],
					"read_clipped_seq": seqinfo["clipped_seq"],
					"read_qual": seqinfo["qual"],
					"read_cigar": seqinfo["cigar_string"],
					"read_reversed_complemented": str(seqinfo["reversed_complemented"]),
					"polya_start_in_genome": str(seqinfo["start"] + p[0]),
					"polya_end_in_genome": str(seqinfo["end"] + p[1] - 1),
					"polya_start_in_read": str(p[0] + 1),
					"polya_end_in_read": str(p[1]),
					"polya_length": str(p[1] - p[0]),
					"polya_nadenosines": str(p[2]),
					"primer_start_in_genome": str(primer_start_in_genome),
					"primer_end_in_genome": str(primer_end_in_genome),
					"primer_start_in_read": str(primer_start_in_read),
					"primer_end_in_read": str(primer_end_in_read),
					"primer_length": str(primer_length),
					"primer_distance": str(primer_distance)}\
					for p in polya]

	return results



def _returnLongestMatch_(near_matches):
	lengths      = [x.end - x.start for x in near_matches]
	max_index    = lengths.index(max(lengths))
	near_matches = near_matches[max_index]
	return near_matches

def detectSubSequence(seq, subseq, min_len = 7, max_l_dist = 1):

	# subsequence could be in two different orientations
	# subseqs = [subseq, revComp(subseq), subseq[::-1], comp(subseq)]
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
	if near_matches:
		near_matches = _returnLongestMatch_(near_matches)

	return near_matches

def retrieveFeaturesFromGTF(gtf_filename, nlines = None):

	'''
	read GTF files and returns exonic features for counting 
	and information about mRNA transcripts
	'''

	# setup progress bar
	widgets =  [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()] 
	gtf_nrow = sum(1 for line in open(gtf_filename)) # get number of rows in file
	bar = progressbar.ProgressBar(widgets = widgets, maxval = gtf_nrow) # init progress bar
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
		bar.update(i)

		# get information from transcript
		if feature.type == "transcript":
			transcript_features[feature.attr["gene_id"]] += [[feature.iv.start, feature.iv.end, feature.iv.length]]

		# get information from exon
		if feature.type == "exon": 
			exon_features[feature.iv] += feature.attr["gene_id"]

	bar.finish()

	# format information for transcripts
	transcript_features = {k: np.array(v) for k, v in transcript_features.items()} # convert to array
	transcript_features = {k: [min(v[:,0]), max(v[:,1]), np.mean(v[:,2])] for k, v in transcript_features.items()}

	return transcript_features, exon_features

def detectPrimerAndPolyA(almnt_filename, 
	exon_features, 
	polya_min_len, 
	polya_max_percent_non_adenosines, 
	polya_seed_len, 
	polya_method,
	primer_seq, 
	primer_min_len,
	primer_max_dist, 
	nlines = None, 
	bam = True):
	
	mask_until = 40

	if bam:
		reader = HTSeq.BAM_Reader
	else:
		reader = HTSeq.SAM_Reader

	# inialize counters and results dictionnary
	total_counts = collections.Counter()
	results      = collections.defaultdict(list)

	# iterate reads in bam files
	if not nlines:
		almnt_nrow = sum(1 for line in HTSeq.pair_SAM_alignments(reader(almnt_filename), bundle=True)) # get number of rows in file
	else:
		almnt_nrow = nlines

	# progress bar
	widgets =  [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()] 
	bar = progressbar.ProgressBar(widgets = widgets, maxval = almnt_nrow) 
	bar.start()

	# iterate through BAM file
	i = 0
	min_len = min(polya_min_len, primer_min_len)
	for bundles in HTSeq.pair_SAM_alignments(reader(almnt_filename), bundle=True):

		# update progress bar
		i += 1
		if nlines: 
			if i > nlines:
				break
		bar.update(i)

		if len(bundles) == 0:
			continue

		bundle = bundles[0] # only use first mapping

		# extract pair
		first_almnt, second_almnt = bundle

		if first_almnt == None or second_almnt == None: 
			print "ERROR: missing alignment pair"
			break

		# get info from reads
		first_seqinfo  = getSeqInfo(first_almnt)
		second_seqinfo = getSeqInfo(second_almnt)

		# retrieve gene id of alignments
		if not first_almnt.aligned and not second_almnt.aligned:
			gene_id = "_unmapped"
		else:
			gene_ids = set()
			if first_almnt.aligned:
				for iv, val in exon_features[first_almnt.iv].steps():
					gene_ids |= val
			if second_almnt.aligned:
				for iv, val in exon_features[second_almnt.iv].steps():
					gene_ids |= val

			if len(gene_ids) == 1:
				gene_id = list(gene_ids)[0]
			elif len(gene_ids) == 0:
				gene_id = "_no_feature"
			else:
				gene_id = "_ambiguous"

		# add read to count 
		total_counts[gene_id] += 1

		# check if perfect math
		first_ignore = len(first_seqinfo["cigar_operations"]) - first_seqinfo["cigar_operations"].count("M") < min_len
		second_ignore = len(second_seqinfo["cigar_operations"]) - second_seqinfo["cigar_operations"].count("M") < min_len

		# remove nucleotides that correspond to matches
		if not first_ignore:
			seq11 = first_seqinfo["clipped_seq"]
		if not second_ignore:
			seq21 = second_seqinfo["clipped_seq"]

		# detect poly-adenlynation (we only look at the 3')
		if not first_ignore:
			first_polya  = detectPolyA(mask5Prime(seq11, mask_until), first_seqinfo["qual"], \
				min_len = polya_min_len, \
				max_percent_non_adenosines = polya_max_percent_non_adenosines, \
				seed_len = polya_seed_len, \
				method = polya_method)

		if not second_ignore:
			second_polya = detectPolyA(mask5Prime(seq21, mask_until), second_seqinfo["qual"], \
				min_len = polya_min_len, \
				max_percent_non_adenosines = polya_max_percent_non_adenosines, \
				seed_len = polya_seed_len, \
				method = polya_method)

		# detect primer in sequence (primer is already reversed complement in function)
		if not first_ignore:
			first_primer = detectSubSequence(seq11, primer_seq, min_len = primer_min_len, max_l_dist = primer_max_dist)
		if not second_ignore:
			second_primer = detectSubSequence(seq21, primer_seq, min_len = primer_min_len, max_l_dist = primer_max_dist)

		# format results
		r = []
		if not first_ignore:
			r += formatResults(first_polya, first_primer, first_seqinfo, "first")
		if not second_ignore:
			r += formatResults(second_polya, second_primer, second_seqinfo, "second")
		results[gene_id] += r

	bar.finish()

	# add total count information to results
	for gene_id in total_counts.keys():
		for i in xrange(len(results[gene_id])):
			results[gene_id][i]["read_count"] = total_counts[gene_id]

	return results

def writeResults(results, transcript_features, outf):

	# write header
	print >> outf, "gene_id", "transcript_start", "transcript_end", "transcript_length", \
	"read_name", "alignment", \
	"chrom", "read_start", "read_end", "read_length",\
	"read_seq", "read_clipped_seq", "read_qual", "read_cigar", "read_reversed_complemented", "read_count", \
	"polya_start_in_genome", "polya_end_in_genome", \
	"polya_start_in_read", "polya_end_in_read", \
	"polya_length", "polya_nadenosines", \
	"primer_start_in_genome", "primer_end_in_genome", \
	"primer_start_in_read", "primer_end_in_read", \
	"primer_length", "primer_dist"

	# write entries for each gene
	for gene_id in sorted(results.keys()):

		# for a specific gene, retrieve information
		entry = results[gene_id]

		# retrieve transcript info
		try:
			transcript_start = transcript_features[gene_id][0] # transcript start position in reference genome
			transcript_end = transcript_features[gene_id][1] # transcript end position
			transcript_length = transcript_features[gene_id][2] # transcript mean length
		except:
			transcript_start = "NA"
			transcript_end = "NA"
			transcript_length = "NA"

		# print one line for every read
		for p in entry:
			print >> outf, gene_id, transcript_start, transcript_end, transcript_length, \
			p["read_name"], p["alignment"], \
			p["read_chrom"], p["read_start"], p["read_end"], p["read_length"], \
			p["read_seq"], p["read_clipped_seq"], p["read_qual"], p["read_cigar"], p["read_reversed_complemented"], p["read_count"],\
			p["polya_start_in_genome"], p["polya_end_in_genome"], \
			p["polya_start_in_read"], p["polya_end_in_read"], \
			p["polya_length"], p["polya_nadenosines"], \
			p["primer_start_in_genome"], p["primer_end_in_genome"], \
			p["primer_start_in_read"], p["primer_end_in_read"], \
			p["primer_length"], p["primer_distance"]

	return 0
