//gcc -o compare_bam_txt compare_bam_txt.c -I/home/ryand/include -L/home/ryand/lib -lhts
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
    htsFile *fp_bam = NULL;
    FILE *fp_txt = NULL;
    bam1_t *read = bam_init1();
    bam_hdr_t *header = bam_hdr_init();
    char *qname, line[1024];
    int rv;
    unsigned long long total_comparisons = 0;

    //Open the files
    if(argc != 3) return -1;
    fp_bam = hts_open(argv[1], "rb");
    fp_txt = fopen(argv[2], "r");

    header = bam_hdr_read(fp_bam->fp.bgzf);
    bam_read1(fp_bam->fp.bgzf, read);
    qname = bam_get_qname(read);
    if(fgets(line, 1024, fp_txt) == NULL) return -1;
    *(line+strlen(line)-1) = '\0';
    while(1) {
        if(++total_comparisons % 10000000 == 0) {fprintf(stderr,"%llu\n", total_comparisons);fflush(stderr);}
        rv = strnum_cmp(line, qname);
        if(rv == 0) {
            if(bam_read1(fp_bam->fp.bgzf, read) <= 4) break;
            qname = bam_get_qname(read);
            if(fgets(line, 1024, fp_txt) == NULL) break;
            *(line+strlen(line)-1) = '\0';
        } else if(rv>0) {
            printf("BAM file is missing %s\n", line); fflush(stdout);
            if(fgets(line, 1024, fp_txt) == NULL) break;
            *(line+strlen(line)-1) = '\0';
        } else {
            printf("TXT file is missing %s\n", qname); fflush(stdout);
            if(bam_read1(fp_bam->fp.bgzf, read)  <= 4) break;
            qname = bam_get_qname(read);
        }
    }

    while(bam_read1(fp_bam->fp.bgzf, read) != 0) {
        qname = bam_get_qname(read);
        printf("TXT file is missing %s\n", qname); fflush(stdout);
    }
    while(fgets(line, 1024, fp_txt) != NULL) {
        *(line+strlen(line)-1) = '\0';
        printf("BAM file is missing %s\n", line); fflush(stdout);
    }

    bam_hdr_destroy(header);
    bam_destroy1(read);
    hts_close(fp_bam);
    fclose(fp_txt);

    return 0;
}
