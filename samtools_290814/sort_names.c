#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <unistd.h>

static int strnum_cmp(const void *_a, const void *_b) {
    const unsigned char *a = *(const unsigned char**)_a, *b = *(const unsigned char**)_b;
    const unsigned char *pa = a, *pb = b;
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
    htsFile *fp = NULL;
    bam1_t *read = bam_init1();
    bam_hdr_t *header = bam_hdr_init();
    FILE *tmp_file = NULL, **tmp_files = NULL;
    char *tmp_file_name = malloc(1024 * sizeof(char));
    char **rnames, **lines;
    size_t j = 0, i = 0, tmp_file_num = 0, maxsize = 10000000;
    size_t nfinished = 0;
    int to_print = -1;
    unsigned long long total = 0;

    //Open the files
    if(argc != 2) return -1;
    fp = hts_open(argv[1], "rb");
    rnames = malloc(sizeof(char *) * maxsize);

    fprintf(stderr,"Reading in %s\n", argv[1]); fflush(stderr);
    header = bam_hdr_read(fp->fp.bgzf);
    sprintf(tmp_file_name, "tmp_file.%04i", (int) tmp_file_num);
    tmp_file = fopen(tmp_file_name, "w");
    while(bam_read1(fp->fp.bgzf, read) >= 0) {
        if((total+1) % 10000000 == 0) { fprintf(stderr, "%llu\n", total+1); fflush(stderr);}
        *(rnames+i) = strdup(bam_get_qname(read));
        assert(*(rnames+i) != NULL);
        i++;
        if(i == maxsize) {
            fprintf(stderr,"\tsorting..."); fflush(stderr);
            qsort(rnames, i, sizeof(char *), strnum_cmp);
            fprintf(stderr, "done\n"); fflush(stderr);
            for(j=0; j<i; j++) {
                fprintf(tmp_file, "%s\n", rnames[j]);
                free(rnames[j]);
            }
            fclose(tmp_file);
            sprintf(tmp_file_name, "tmp_file.%04i", (int) ++tmp_file_num);
            tmp_file = fopen(tmp_file_name, "w");
            i=0;
        }
    }
    if(i != 0) {
        qsort(rnames, (size_t) i, sizeof(char *), strnum_cmp);
        for(j=0; j<i; j++) {
            fprintf(tmp_file, "%s\n", rnames[j]);
            free(rnames[j]);
        }
    }
    fprintf(stderr, "Finished reading in read names\n"); fflush(stderr);
    free(rnames);
    bam_hdr_destroy(header);
    bam_destroy1(read);
    hts_close(fp);
    fclose(tmp_file);

    tmp_file_num = 122;
    //Open tmp files
    tmp_files = malloc(sizeof(FILE*) * tmp_file_num);
    lines = malloc(sizeof(char *) * tmp_file_num);
    for(i=0; i<tmp_file_num; i++) {
        lines[i] = malloc(sizeof(char) * 1024);
        sprintf(tmp_file_name, "tmp_file.%04i", (int) i);
        tmp_files[i] = fopen(tmp_file_name, "r");
        fgets(lines[i], 1024, tmp_files[i]);
        *(lines[i] + strlen(lines[i])-1) = '\0'; //Trim \n
    }

    //Merge
    while(nfinished < tmp_file_num) {
        nfinished = 0;
        to_print = 0;
        for(i=1; i<tmp_file_num; i++) {
            if(lines[i] == NULL) { //File i is finished
                nfinished++;
                continue;
            }
            if(lines[to_print] == NULL) { //File 0 is finished
                to_print = i;
                nfinished++;
                continue;
            }
            if(strnum_cmp((void *) (lines+to_print), (void *) (lines+i)) > 0) to_print = i;
        }
        printf("%s\n", lines[to_print]);
        if(fgets(lines[to_print], 1024, tmp_files[to_print]) == NULL) {
            free(lines[to_print]);
            lines[to_print] = NULL;
            nfinished++;
        } else {
            *(lines[to_print] + strlen(lines[to_print])-1) = '\0'; //Trim \n
        }
    }

    //Free things
    for(i=0; i<=tmp_file_num; i++) {
        fclose(tmp_files[i]);
//        sprintf(tmp_file_name, "tmp_file.%04i", (int) i);
//        unlink(tmp_file_name);
    }
    free(lines);
    free(tmp_file_name);
    free(tmp_files);

    return 0;
}
