Problem: A reported bug in samtools sort (multiple versions) resulted in loss of alignments after query-sorting a BAM file.

The BAM file exhibiting this problem was privately made available and I query sorted it with samtools 0.1.19 and 1.0. I manually checked for what was reported to have been a lost alignment but was able to find it in both sorted BAMs. Consequently, I wrote the following two (slow) C programs. The first, `sort_names.c`, query sorts the read names in a given BAM file. The original BAM file had ~1.2 billion alignments, so this was quite slow. The second program, `compare_bam_txt.c`, compares the resulting TXT file with a query sorted BAM file and reports any missing alignments. Compilation is as follows:

`gcc -g -o sort_names sort_names.c -I/home/ryand/include -L/home/ryand/lib -lhts`
`gcc -g -o compare_bam_txt compare_bam_txt.c -I/home/ryand/include -L/home/ryand/lib -lhts`
