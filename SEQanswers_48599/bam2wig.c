#include "htslib/sam.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>

typedef struct {
    samFile *fp;
    bam_hdr_t *hdr;
    int minMapq;
} mplp_data;

int filter_func(void *data, bam1_t *b) {
    int rv;
    mplp_data *ldata = (mplp_data *) data;

    while(1) {
        rv = sam_read1(ldata->fp, ldata->hdr, b);
        if(rv<0) return rv;
        if(b->core.tid == -1 || b->core.flag & BAM_FUNMAP) continue; //Unmapped
        if(b->core.qual < ldata->minMapq) continue; //-q
        if(b->core.flag & 0x300) continue; //Ignore secondary alignments and those with QC failed
        if(b->core.flag & BAM_FDUP) continue;
        if((b->core.flag & 0x9) == 0x1) b->core.flag |= 0x2; //Discordant pairs can cause double counts
        break;
    }
    return rv;
}

void processDepths(samFile *ifile, bam_hdr_t *hdr, int minMapq, int minPhred, FILE *of1, FILE *of2) {
    bam_mplp_t iter;
    int tid, pos, i, n_plp, ret;
    const bam_pileup1_t **plp = NULL;
    mplp_data *data = NULL;
    uint32_t plus_depth, minus_depth;
    int32_t ctid = -1;

    data = calloc(1,sizeof(mplp_data));
    if(data == NULL) {
        fprintf(stderr, "Couldn't allocate space for the data structure in processDepths()!\n");
        return;
    }
    data->fp = ifile;
    data->hdr = hdr;
    plp = calloc(1, sizeof(bam_pileup1_t *));
    if(plp == NULL) {
        fprintf(stderr, "Couldn't allocate space for the plp structure in extractCalls()!\n");
        return;
    }

    //Start the pileup
    iter = bam_mplp_init(1, filter_func, (void **) &data);
    bam_mplp_init_overlaps(iter);
    bam_mplp_set_maxcnt(iter, 1e6);
    while((ret = bam_mplp_auto(iter, &tid, &pos, &n_plp, plp)) > 0) {
        plus_depth = 0; minus_depth = 0;
        if(tid != ctid) {
            ctid = tid;
            fprintf(of1, "variableStep chrom=%s span=1\n", hdr->target_name[tid]);
            fprintf(of2, "variableStep chrom=%s span=1\n", hdr->target_name[tid]);
        }

        for(i=0; i<n_plp; i++) {
            if(bam_get_qual(plp[0][i].b)[plp[0][i].qpos] < minPhred) continue;
            if(plp[0][i].b->core.flag & 0x80) { //Read #2
                if(plp[0][i].b->core.flag & 0x10) { //- strand
                    minus_depth++;
                } else { //+ strand
                    plus_depth++;
                }
            } else {
                if(plp[0][i].b->core.flag & 0x10) { //+ strand
                    plus_depth++;
                } else { //- strand
                    minus_depth++;
                }
            }
        }
        fprintf(of1, "%"PRId32"\t%"PRIu32"\n", pos+1, plus_depth);
        fprintf(of2, "%"PRId32"\t%"PRIu32"\n", pos+1, minus_depth);
    }

    bam_mplp_destroy(iter);
    free(data);
    free(plp);
}

void usage(char *prog) {
    fprintf(stderr, "Usage: %s [OPTIONS] file.bam prefix\n", prog);
    fprintf(stderr, "\nThis program accepts a sorted and indexed BAM file containing RNAseq alignments\n"
"and produces two WIG files (prefix.plus.wig and prefix.minus.wig) with the\n"
"per-base depth of the appropriate strand.\n");
    fprintf(stderr, "\n"
"OPTIONS:\n"
"-s INT Indicate which of the reads in a pair denote the strand. The default is\n"
"       2, meaning that read #2 denotes the strand of a pair. For single-end\n"
"       data, this means that if an alignment is reverse complemented then it\n"
"       will increase the depth of the + strand in the appropriate region. This\n"
"       is appropriate for dUTP-based libraries, which are the most common.\n"
"-q INT Minimum MAPQ value to process an alignment. The default is 0, meaning\n""       that all alignments will be included.\n"
"-p INT Minimum Phred score to include a base in an alignment. The default is 1\n"
"       and that's also the lowest possible value (if a value of 0 were used,\n"
"       then depths would be artificially inflated when reads in a pair overlap)\n");
}

int main(int argc, char *argv[]) {
    int strand = 2, c;
    int minMapq = 0, minPhred = 1;
    samFile *ifile;
    bam_hdr_t *hdr;
    char *oname;
    FILE *of1, *of2;

    static struct option lopts[] = {
        {"help", 0, NULL, 'h'}
    };

    while((c = getopt_long(argc, argv, "s:q:p:h", lopts, NULL)) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
            break;
        case 's' :
            strand = atoi(optarg);
            if(strand != 1 && strand != 2) {
                fprintf(stderr, "-s must be either 1 or 2. We'll use the default of 2.\n");
                fflush(stderr);
                strand = 2;
             }
             break;
        case 'q' :
            minMapq = atoi(optarg);
            if(minMapq < 0) minMapq = 0;
            break;
        case 'p' :
             minPhred = atoi(optarg);
             if(minPhred < 1) minPhred = 1;
             break;
        default :
             fprintf(stderr, "Invalid option '%c'\n", c);
            usage(argv[0]);
            return 1;
        }
    }

    if(argc-optind != 2) {
        usage(argv[0]);
        return 1;
    }

    //Open the files
    ifile = sam_open(argv[optind], "rb");
    hdr = sam_hdr_read(ifile);
    oname = malloc(sizeof(char)*(strlen(argv[optind+1]) + strlen(".minus.wig ")));
    assert(oname);
    sprintf(oname, "%s.plus.wig", argv[optind+1]);
    if(strand == 2) of1 = fopen(oname, "w");
    else of2 = fopen(oname, "w");
    sprintf(oname, "%s.minus.wig", argv[optind+1]);
    if(strand == 2) of2 = fopen(oname, "w");
    else of1 = fopen(oname, "w");
    free(oname);

    processDepths(ifile, hdr, minMapq, minPhred, of1, of2);

    sam_close(ifile);
    fclose(of1);
    fclose(of2);
    return 0;
}
