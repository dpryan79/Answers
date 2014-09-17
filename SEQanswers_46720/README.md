People often ask about recreating the metrics output by bowtie2 at the end of its run. The `bam2bowtie_stats.py` program is a simple python script that uses pysam to parse and reproduce these metrics. Note that it currently only supports single-end reads, though I could trivially change that. You will obviously need python and pysam installed to use this.

Usage is quite simple. Simply type `bam2bowtie_stats.py file.bam` after making the file executable and moving it into your PATH.
