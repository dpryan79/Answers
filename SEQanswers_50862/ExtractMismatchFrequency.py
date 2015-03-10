#!/usr/bin/env python
import pysam
import argparse
import sys

#Get the number of mismatches
def getNM(b) :
    if(b.has_tag("NM")) :
        return b.get_tag("NM")
    elif(b.has_tag("nM")) :
        return b.get_tag("nM")
    else :
        return 0

#parse the CIGAR string, summing the M=X elements
def getNumAligned(b) :
    l = 0
    for t in b.cigartuples :
        if(t[0] == 0 or t[0] == 7 or t[0] == 8) :
            l += t[1]
    return l

parser = argparse.ArgumentParser(description="Extract the number of NM bases and the number of aligned bases from a BAM file")
parser.add_argument("fname", metavar="alignments.bam", help="Input BAM file")
args = parser.parse_args()

if(args.fname == None) :
    parser.print_help()
    sys.exit()

samfile = pysam.AlignmentFile(args.fname, "rb")
NM = 0
bpAligned = 0
while 1 :
    try :
        b = samfile.next()
    except :
        break

    #Ignore secondard and supplemental alignments
    if(b.flag & 0x800 or b.flag & 0x100) :
        continue

    NM += getNM(b)
    bpAligned += getNumAligned(b)

samfile.close()

print("%i mismatches in %i aligned bases (%5.2f%% mismatch rate)" % (NM, bpAligned, float(100*NM)/float(bpAligned)))
