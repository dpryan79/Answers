Problem: A reported bug in samtools sort (multiple versions) resulted in loss of alignments after query-sorting a BAM file.

The BAM file exhibiting this problem was privately made available and I query sorted it with samtools 0.1.19 and 1.0. I manually checked for what was reported to have been a lost alignment but was able to find it in both sorted BAMs. Consequently, I wrote the following three (slow) C programs. The first, `sort_names.c`, query sorts the read names in a given BAM file. The original BAM file had ~1.2 billion alignments, so this was quite slow. The second program, `compare_bam_txt.c`, compares the resulting TXT file with a query sorted BAM file and reports any missing alignments. The third program, `compare_bam_bam` compares two name-sorted BAM files. The two compare programs are quick enough, the sort program is very slow, but since this is a one-off sort of thing there's not much reason to optimize things.

Compilation is as follows:

`gcc -g -o sort_names sort_names.c -I/home/ryand/include -L/home/ryand/lib -lhts`
`gcc -g -o compare_bam_txt compare_bam_txt.c -I/home/ryand/include -L/home/ryand/lib -lhts`
`gcc -g -o compare_bam_bam compare_bam_bam.c -I/home/ryand/include -L/home/ryand/lib -lhts`

Usage:
`./sort_names unsorted.bam > sorted_names`
`./compare_bam_txt name_sorted.bam sorted_names > differences`
`./compare_bam_bam name_sorted1.bam name_sorted2.bam > differences`

The last program is useful for comparing name sorted files produced by different versions of samtools (namely, 0.1.19 and 1.0, since they should produce the same output).
