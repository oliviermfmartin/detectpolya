#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import detectpolya

if __name__ == "__main__":

	# handle arguments
	parser = argparse.ArgumentParser(description='Detecting poly-adenylation sites in FASTA and BAM files from next-generation sequencing')
	
	parser.add_argument('input', type=str, help='Input file path (FASTA, SAM or BAM)')
	parser.add_argument('output', type=str, help='Output file')
	parser.add_argument('-f', '--filetype', type=str, default="bam", help='Filetype: "fa", "fq", "sam", "bam" (default: bam)')
	parser.add_argument('-d', '--paired_ends', action='store_false', help='For SAM and BAM files specify paired-ends')
	parser.add_argument('-g', '--gtf', default=None, type=str, help='GTF filename')
	parser.add_argument('-m', '--method', type=str, default="seed", help='')
	parser.add_argument('-l', '--min_len', type=int, default=5, help='Minimum length of poly-A tail')
	parser.add_argument('-n', '--max_prop_non_a', type=float, default=0.2, help='Maximum proportion of non-adenosines')
	parser.add_argument('-s', '--seed_len', type=int, default=3, help='Seed length of seed-based poly-A detection')
	parser.add_argument('-p', '--primer', action='store_false', help='Should script also look for primer')
	parser.add_argument('--primer_seq', type=str, help='Primer sequence')
	parser.add_argument('--primer_min_len', type=int, default=10, help='Minimum length of primer match')
	parser.add_argument('--primer_max_dist', type=float, default=1, help='Maximum Levenshtein distance between primer and read subsequence')
	parser.add_argument('--max_nlines', type=int, default=None, help='Maximum number of lines to be read from input file')
	parser.add_argument('--max_nlines_gtf', type=int, default=None, help='Maximum number of lines to be read from GTF file')

	args = parser.parse_args()

	# open file to write
	outf = open(args.output, "w")

	# read counts, poly-adenlynation and primer
	print "Detection algorithm running"
	results = detectpolya.analyseFile(args.input,
		gtf_filename         = args.gtf,
		polya_min_len        = args.min_len,
		polya_max_prop_non_a = args.max_prop_non_a,
		polya_seed_len       = args.seed_len,
		primer_seq           = args.primer_seq,
		primer_min_len       = args.primer_min_len,
		primer_max_dist      = args.primer_max_dist,
		polya_method         = args.method, 
		nlines               = args.max_nlines,
		gtf_nlines           = args.max_nlines_gtf,
		filetype             = args.filetype)

	# write results to file
	print "Writing results"
	detectpolya.writeResults(results, outf)

	outf.close()
	print "Done!"
