#!/usr/bin/python
# -*- coding: utf-8 -*-

if __name__ == "__main__":

	import argparse
	import detectpolya

	# handle arguments
	parser = argparse.ArgumentParser(description='Detecting poly-adenylation sites in FASTA and BAM files from next-generation sequencing')
	
	parser.add_argument('-i', '--input',  required=True, type=str, help='Input file path (FASTA, FASTQ, SAM or BAM)')
	parser.add_argument('-i2', '--input2', type=str, default=None, help='Second input file path for pair-end FASTA, FASTQ')
	parser.add_argument('-o', '--output', type=str, default=None, help='Output file (default is console).')
	parser.add_argument('-t', '--filetype', type=str, default="bam", help='Filetype: "fa", "fq", "sam", "bam" (default: bam)')
	parser.add_argument('-d', '--paired_ends', action='store_true', help='Specifies if experiment used paired ends')
	parser.add_argument('-g', '--gtf', default=None, type=str, help='GTF filename')
	parser.add_argument('-m', '--method', type=str, default="seed", help='Poly-A tail detection method: "seed" or "window" (default: seed)')
	parser.add_argument('-l', '--min_len', type=int, default=5, help='Minimum length of poly-A tail')
	parser.add_argument('-n', '--max_prop_non_a', type=float, default=0.2, help='Maximum proportion of non-adenosines')
	parser.add_argument('-s', '--seed_len', type=int, default=3, help='Seed length for seed-based poly-A detection')
	parser.add_argument('-p', '--primer', action='store_true', help='Should script also look for primer?')
	parser.add_argument('--primer_seq', type=str, help='Primer sequence')
	parser.add_argument('--primer_min_len', type=int, default=10, help='Minimum length of primer match')
	parser.add_argument('--primer_max_dist', type=int, default=1, help='Maximum Levenshtein distance between primer and read subsequence')
	parser.add_argument('--max_nlines', type=int, default=None, help='Maximum number of lines to be read from input file')
	parser.add_argument('--max_nlines_gtf', type=int, default=None, help='Maximum number of lines to be read from GTF file')
	parser.add_argument('--noheader', action='store_true', default=None, help='Don\'t print header')
 	parser.add_argument('--silent', action='store_true', default=None, help='Don\'t print progress')

	args = parser.parse_args()

	if (args.filetype == "fa" or args.filetype == "fq") and args.paired_ends:
		assert args.input2 != None, "Must specify second mates with --input2 for pair-end FASTA and FASTQ files."
		infiles = [args.input, args.input2]
	else:
		infiles = args.input

	# open file to write
	if args.output != None:
		outf = open(args.output, "w")
	else:
		outf = None

	# read counts, poly-adenlynation and primer
	results = detectpolya.analyseFile(infiles,
		filetype             = args.filetype,
		paired_ends          = args.paired_ends,
		gtf_filename         = args.gtf,
		polya_method         = args.method, 
		polya_min_len        = args.min_len,
		polya_max_prop_non_a = args.max_prop_non_a,
		polya_seed_len       = args.seed_len,
		primer_seq           = args.primer_seq,
		primer_min_len       = args.primer_min_len,
		primer_max_dist      = args.primer_max_dist,
		nlines               = args.max_nlines,
		gtf_nlines           = args.max_nlines_gtf,
		verbose              = not args.silent)

	# write results to file
	if not args.silent:
		print "Writing results"
	detectpolya.printResults(results, outf = outf, header = not args.noheader)

	if args.output != None:
		outf.close()
