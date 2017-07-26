#!/usr/bin/env python
import pysam
import argparse
import sys
import gzip


def parseArgs():
    parser = argparse.ArgumentParser(description="Convert a paired-end BAM file to fastq files after filtering out pairs where both align to a list of filtered contigs.")
    parser.add_argument("BAMfile", help="The BAM file to parse")
    parser.add_argument("filterList", help="A file containing a list of contigs/scaffolds/chromosomes to exclude (one per line). Pairs of reads will only be excluded from output if BOTH mates align to one of these.")
    parser.add_argument("outputPrefix", help="The prefix for the output files. They will be names prefix_R1.fastq.gz and prefix_R2.fastq.gz")
    return parser


def loadFilterList(fname):
    s = set()
    for line in open(fname):
        s.add(line.strip())
    return s


def filterAlignment(b, filterList):
    """
    Exclude only cases where BOTH mates align to a filtered scaffold/contig
    Additionally exclude secondary/supplemental alignments
    """
    if b.flag & 2304:
        return True
    if b.is_unmapped or b.mate_is_unmapped:
        return False
    if b.reference_name in filterList or b.next_reference_name in filterList:
        return True
    return False


def inBuffer(b, buf):
    if b.query_name in buf:
        return True
    return False


def revComp(s):
    d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    s = [d[c] for c in s]
    return ''.join(s[::-1])


def storeInBuffer(b, buf):
    s = b.query_sequence
    q = ''.join([chr(c+33) for c in b.query_qualities])
    if b.is_read1:
        rNum = 1
    else:
        rNum = 2
    if b.is_reverse:
        s = revComp(s)
        q = q[::-1]
    buf[b.query_name] = (s, q, rNum)


def dumpAlignments(b, buf, of1, of2):
    b2 = buf[b.query_name]

    # Write the mate
    of = of1
    if b2[2] == 2:
        of = of2
    of.write("@{}\n{}\n+\n{}\n".format(b.query_name, b2[0], b2[1]))
    
    # Write the read
    s = b.query_sequence
    q = ''.join([chr(c+33) for c in b.query_qualities])
    if b.is_reverse:
        s = revComp(s)
        q = q[::-1]
    of = of1
    if b.is_read2:
        of = of2
    of.write("@{}\n{}\n+\n{}\n".format(b.query_name, s, q))

    # Remove from the buffer
    del buf[b.query_name]


def main(args):
    args = parseArgs().parse_args(args)

    bam = pysam.AlignmentFile(args.BAMfile)
    filterList = loadFilterList(args.filterList)
    of1 = gzip.GzipFile("{}_R1.fastq.gz".format(args.outputPrefix), "w")
    of2 = gzip.GzipFile("{}_R2.fastq.gz".format(args.outputPrefix), "w")

    buf = dict()
    for b in bam:
        if filterAlignment(b, filterList):
            continue
        if inBuffer(b, buf):
            dumpAlignments(b, buf, of1, of2)
        else:
            storeInBuffer(b, buf)

    bam.close()
    of1.close()
    of2.close()


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)
