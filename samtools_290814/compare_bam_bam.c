//gcc -o compare_bam_bam compare_bam_bam.c -I/home/ryand/include -L/home/ryand/lib -lhts
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <string.h>

int strnum_cmp(char *a, char *b) {
    char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

int main(int argc, char *argv[]) {
    htsFile *fp1 = NULL;
    htsFile *fp2 = NULL;
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    bam_hdr_t *header1 = bam_hdr_init();
    bam_hdr_t *header2 = bam_hdr_init();
    int rv;
    unsigned long long total_comparisons = 0;

    //Open the files
    if(argc != 3) return -1;
    fp1 = hts_open(argv[1], "rb");
    fp2 = hts_open(argv[2], "rb");

    header1 = bam_hdr_read(fp1->fp.bgzf);
    header2 = bam_hdr_read(fp2->fp.bgzf);
    bam_read1(fp1->fp.bgzf, read1);
    bam_read1(fp2->fp.bgzf, read2);
    while(1) {
        if(++total_comparisons % 10000000 == 0) {fprintf(stderr,"%llu\n", total_comparisons);fflush(stderr);}
        rv = strnum_cmp(bam_get_qname(read1), bam_get_qname(read2));
        if(rv == 0) {
            if(bam_read1(fp1->fp.bgzf, read1) <= 4) break;
            if(bam_read1(fp2->fp.bgzf, read2) <= 4) break;
        } else if(rv>0) {
            fprintf(stderr, "A %s has %s but %s has %s\n", argv[1], bam_get_qname(read1), argv[2], bam_get_qname(read2)); fflush(stderr);
            printf("A %s has %s but %s has %s\n", argv[1], bam_get_qname(read1), argv[2], bam_get_qname(read2)); fflush(stdout);
            if(bam_read1(fp2->fp.bgzf, read2) <= 4) break;
        } else { //rv<0
            fprintf(stderr, "B %s has %s but %s has %s\n", argv[1], bam_get_qname(read1), argv[2], bam_get_qname(read2)); fflush(stderr);
            printf("B %s has %s but %s has %s\n", argv[1], bam_get_qname(read1), argv[2], bam_get_qname(read2)); fflush(stdout);
            if(bam_read1(fp1->fp.bgzf, read1) <= 4) break;
        }
    }

    while(bam_read1(fp1->fp.bgzf, read1) > 4) {
        fprintf(stderr, "C %s has %s but %s has finished\n", argv[1], bam_get_qname(read1), argv[2]); fflush(stderr);
        printf("C %s has %s but %s has finished\n", argv[1], bam_get_qname(read1), argv[2]); fflush(stdout);
    }
    while(bam_read1(fp2->fp.bgzf, read2) > 4) {
        fprintf(stderr, "D %s has finished but %s has %s\n", argv[1], argv[2], bam_get_qname(read2)); fflush(stderr);
        printf("D %s has finished but %s has %s\n", argv[1], argv[2], bam_get_qname(read2)); fflush(stdout);
    }

    bam_hdr_destroy(header1);
    bam_hdr_destroy(header2);
    bam_destroy1(read1);
    bam_destroy1(read2);
    hts_close(fp1);
    hts_close(fp2);

    return 0;
}
