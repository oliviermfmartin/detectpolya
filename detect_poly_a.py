#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import progressbar
import re
import collections
import numpy as np
import HTSeq
import fuzzysearch
import pdb
from detect_polya_functions import *

almnt_filename = "./data/alignments/c_elegans/CCB68ANXX_4_N703-S503-NX-xt.sorted.bam" # name of BAM file
# almnt_filename = "ribopolya.sam"
gtf_filename = "../../data/c_elegans.PRJNA13758.WS265.canonical_geneset.gtf" # name of GTF file
output_filename = "./tables/polya_window.csv"

polya_min_len = 5
polya_max_percent_non_adenosines = 0.20
polya_seed_len = 3
polya_method = "window"

primer_seq = "AAGCAGTGGTATCAACGCAGAGTAC"
primer_min_len = 10
primer_max_dist = 1

nlines_gtf =  None
# nlines_gtf = 1
nlines_bam = 4816016
# nlines_bam = 30000

# seq  = "=ATGTAATGTCTAACAAAAA"
# qual = "KKKKKKKKKKKKKK!KKKKK"
# # print detectPolyA(seq, qual =  qual, method = "window")
# print detectPolyA(seq, qual =  qual, method = "seed")
# sys.exit(1)

# open file to write
outf = open(output_filename, "w")

# retrieve features from GTF file
print "reading GTF file"
transcript_features, exon_features = retrieveFeaturesFromGTF(gtf_filename, nlines = nlines_gtf)

# read counts, poly-adenlynation and primer
print "reading BAM file"
results = detectPrimerAndPolyA(almnt_filename, \
	exon_features   = exon_features,
	polya_min_len   = polya_min_len, \
	polya_max_percent_non_adenosines = polya_max_percent_non_adenosines, \
	polya_seed_len  = polya_seed_len, \
	primer_seq      = primer_seq, \
	primer_min_len  = primer_min_len, \
	primer_max_dist = primer_max_dist, \
	polya_method    = polya_method, 
	nlines          = nlines_bam, \
	bam             = ".bam" in almnt_filename)

# write results to file
print "writing results to file"
writeResults(results, transcript_features, outf)

outf.close()
print "done!"
