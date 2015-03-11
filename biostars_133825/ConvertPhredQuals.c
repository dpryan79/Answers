#include "../htslib/htslib/sam.h"
#include "../htslib/htslib/hts.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

int convert(htsFile *of, htsFile *fp, bam_hdr_t *hdr) {
    bam1_t *b = bam_init1();
    uint8_t *QUAL;
    int i;

    while(sam_read1(fp, hdr, b) >= 0) {
        QUAL = bam_get_qual(b);
        if(*QUAL != 0xFF) {
            for(i=0; i<b->core.l_qseq; i++) {
                if(QUAL[i] < 64) {
                    fprintf(stderr, "This file seems to have proper QUAL encoding! Exiting...\n");
                    return 1;
                }
                QUAL[i] -= 31;
            }
        }
        sam_write1(of, hdr, b);
    }

    bam_destroy1(b);

    return 0;
}

void usage(char *prog) {
    fprintf(stderr, "Usage: %s [OPTIONS] <file.bam>\n", prog);
    fprintf(stderr, "\n"
"This program accepts a SAM, BAM, or CRAM file as input and will convert its QUAL\n"
"field to be phred+33 encoded if it's incorrectly phred+64 encoded. Note that if\n"
"you input a phred+33 encoded file that the results are likely to be gibberish\n"
"(though the program will likely catch this). By default, the output file is in\n"
"the same format as the input, except for SAM input, which will produce BAM\n"
"output. See the options to change this. Like samtools, the resulting file is\n"
"written to standard output by default (but see the -o option).\n"
"\n"
"<file.bam>   The input file in SAM, BAM or CRAM format. For CRAM format, you\n"
"             must also specify the reference fasta file (see -f, below).\n" 
"\nOPTIONS\n"
"-f  FILE     A reference sequence in fasta format. This is only needed with\n"
"             CRAM input/output.\n"
"-@  INT      The number of output compression threads (default: 1).\n"
"-C           Produce a CRAM file, regardless of the input format (requires -f).\n"
"-B           Produce a BAM file, regardless of the input format.\n"
"-o  STR      Output file name (default: stdout).\n");
}

int main(int argc, char *argv[]) {
    htsFile *fp, *of;
    char *oname = "-", *fa = NULL;
    int nt = 1, err;
    bam_hdr_t *hdr;
    char c;
    int outType = 0; //1: BAM, 2: CRAM

    opterr = 0;
    while((c = getopt(argc, argv, "hf:@:CBo:")) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
        case 'f' :
            fa = optarg;
            break;
        case '@' :
            nt = atoi(optarg);
            if(nt<1) nt=1;
            break;
        case 'C' :
            outType = 2;
            break;
        case 'B' :
            outType = 1;
            break;
        case 'o' :
            oname = optarg;
            break;
        case '?' :
        default :
            fprintf(stderr, "Invalid option '%s'\n", argv[optind-1]);
            usage(argv[0]);
            return 1;
            break;
        }
    }

    if(argc == 1) {
        usage(argv[0]);
        return 0;
    }
    if(argc-optind != 1) {
        fprintf(stderr, "You must specify an input file!\n");
        usage(argv[0]);
        return 1;
    }

    //Open input and output files
    fp = sam_open(argv[optind], "r");
    if((outType==2 || fp->is_cram) && !fa) {
        fprintf(stderr, "When either the input or output is in CRAM format, you must specify a reference sequence with -f!\n");
        return 1;
    }
    hdr = sam_hdr_read(fp);
    if(outType) {
        if(outType==2) {
            of = sam_open(oname, "wc");
        } else {
            of = sam_open(oname, "wb");
        }
    } else {
        if(fp->is_cram) {
            of = sam_open(oname, "wc");
        } else {
            of = sam_open(oname, "wb");
        }
    }
    if(fa) {
        hts_set_fai_filename(fp, fa);
        hts_set_fai_filename(of, fa);
    }
    if(nt>1) hts_set_threads(of, nt);
    sam_hdr_write(of, hdr);

    err = convert(of, fp, hdr);

    return err;
}
