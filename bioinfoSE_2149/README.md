This file is based off of [this question](https://bioinformatics.stackexchange.com/questions/2149/how-to-safely-and-efficiently-convert-subset-of-bam-to-fastq) on the bioinformatics stackexchange. In short, it's needed to convert a paired-end BAM file to two fastq files. The caveats are:

 1. The BAM file may or may not be sorted in any particular way.
 2. There may be duplicate entries (e.g., secondary alignments).
 3. Pairs should be excluded if **BOTH** mates align to a given list of contigs.

This solution uses pysam and creates a buffer with read entries (much like HTSeq-count). As mates are found, they're written to gzipped fastq files and removed from the buffer. This solution does not require any particular BAM file sorting, anything will work.

Note that having duplicate read names that do not arise from secondary or supplemental alignments will break this script!

Usage: `./convert.py alignments.bam filterList.txt /some/output/path/prefix`
