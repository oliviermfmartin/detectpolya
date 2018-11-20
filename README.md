# detectpolya

## Description

This folder contains a python package and a python script used to detect poly-adenylatoin sites and primer sequences in files from next-generation sequencing. Currently FASTA, FASTQ, SAM and BAM files are supported.

## Dependencies

The package was developped for Python 2.7 and requires the following packages: `Bio`, `collections`, `fuzzysearch`,`HTSeq`, `itertools`, `numpy`, `warnings`. These can be easily installed using `pip` or `conda`.

## Installation

You can install the Python by typing in the command line.

```
python setup.py install
```

## Usage

After installation, you can run a Python script `./detectpolya.py`. For information on how to use it, type: `python ./detectpolya.py --help` As an example, you may use one the pair end SAM test file as follows.

```
python ./detectpolya.py --input ./test/paired_ends.sam -t "sam" --paired_ends
```

You can also directly use the Python package after importing it within a Python script: `import detectpolya`. Documentation can be found in `./doc/build/index.html`

The documentation was made using Sphinx and can be rebuild by typing: `sphinx-build -b html doc/source/ doc/build/`

## Column specification

gene_id: Gene ID (SAM or BAM and GTF)
transcript_start: Start of position of most upstream transcript (SAM or BAM and GTF)
transcript_end: End position of most downstream transcript (SAM or BAM and GTF)
transcript_length: Mean length of transcripts (SAM or BAM and GTF)
read_name: Name of read 
read_mate: Number of mate if paired end
read_chrom: Chromosome read maps to (SAM or BAM)
read_start: Start position read maps to (SAM or BAM)
read_end: End position read maps to (SAM or BAM)
read_length: Lenght of read 
read_seq: Sequence of sequenced read (not reversed complemented to match reference)
read_clipped_seq: Sequence of sequenced read with matched nucleotides masked by equal sign (SAM or BAM)
read_qual": Quality string of read in the same order as the sequence (FASTQ, SAM or BAM)
read_cigar: CIGAR string in the same order as the sequence (SAM or BAM)
read_reversed_complemented": Was the read reversed complemented to match sequence (SAM or BAM)
polya_start_in_genome: Start position of poly-A tail in genome (SAM or BAM)
polya_end_in_genome: End position of poly-A tail in genome (SAM or BAM)
polya_start_in_read: Start position of poly-A tail in read
polya_end_in_read: End position of poly-A tail in read
polya_length: Length of poly-A tail
polya_score: Number of expected adenosines
primer_start_in_genome: Start position of primer in genome (SAM or BAM)
primer_end_in_genome: End position of primer in genome (SAM or BAM)
primer_start_in_read: Start position of primer in read
primer_end_in_read: Start position of primer in read
primer_length: Length of primer
primer_score: Number of matched nucleotides

## Author Info

Olivier M. F. Martin, 2018
oliviermfmartin@gmail.com
