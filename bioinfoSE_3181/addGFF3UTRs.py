#!/usr/bin/env python
import argparse


def parseAttributes(kvps):
    d = dict()
    for kvp in kvps:
        k, v = kvp.split("=")
        d[k] = v
    return d


def parseGFF3(fname, topLevel="gene", secondLevel="mRNA", CDS="CDS", exon="exon"):
    f = open(fname)
    genes = dict()
    transcripts = dict()
    g2t = dict()  # A dictionary with gene IDs associated to their transcript IDs

    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        attribs = parseAttributes(cols[8].split(";"))
        ID = attribs["ID"]

        if cols[2] == topLevel:
            genes[ID] = line
            g2t[ID] = []
        elif cols[2] == secondLevel:
            parent = attribs["Parent"]
            g2t[parent].append(ID)

            # Generate a list of [line, [exon entries], [CDS entries]]
            if ID not in transcripts:
                transcripts[ID] = [line, [], []]
            else:
                transcripts[ID][0] = line
        else:
            parent = attribs["Parent"]
            if parent not in transcripts:
                transcripts[parent] = [None, [], []]
            if cols[2] == exon:
                transcripts[parent][1].append((int(cols[3]), int(cols[4]), line))
            elif cols[2] == CDS:
                transcripts[parent][2].append((int(cols[3]), int(cols[4]), line))
    f.close()

    return genes, transcripts, g2t


def findUTRs(transcripts, fivePrime = True):
    for k, t in transcripts.items():
        # If there are no CDS entries then skip
        if len(t[2]) == 0:
            continue
        # Get the strand ("." will be treated as "+")
        strand = t[0].split("\t")[6]

        # For 3' UTR, just swap the strand
        if not fivePrime:
            strand = "+" if strand == "-" else "-"

        exons = [(s, e) for s, e, _ in t[1]]
        CDSs = [(s, e) for s, e, _ in t[2]]
        exons.sort()
        CDSs.sort()
        UTRs = []
        if strand != "-":
            final = CDSs[0][0]
            for s, e in exons:
                if e < final:
                    UTRs.append((s,e))
                elif s < final:
                    UTRs.append((s, final - 1))
                else:
                    break
        else:
            final = CDSs[-1][-1]
            for s, e in exons:
                if e < final:
                    continue
                elif s > final:
                    UTRs.append((s, e))
                else:
                    UTRs.append((final + 1, e))
        t.append(UTRs)


def saveGFF(oname, genes, transcripts, g2t, fivePrime, threePrime):
    o = open(oname, "w")
    o.write("##gff-version 3\n")
    for geneID, geneLine in genes.items():
        o.write(geneLine)
        # Go through each transcript
        for transID in g2t[geneID]:
            t = transcripts[transID]
            if t[0] is None:
                continue
            o.write(t[0])
            cols = t[0].strip().split("\t")

            # Write the exons, then CDS, then 5'UTR, then 3'UTR
            for exon in t[1]:
                o.write(exon[2])
            for CDS in t[2]:
                o.write(CDS[2])
            for idx, UTR in enumerate(t[3]):
                o.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t".format(cols[0], cols[1], fivePrime, UTR[0], UTR[1], cols[6]))
                o.write("ID={}.fivePrimeUTR{};Parent={}\n".format(transID, idx, transID))
            for idx, UTR in enumerate(t[4]):
                o.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t".format(cols[0], cols[1], threePrime, UTR[0], UTR[1], cols[6]))
                o.write("ID={}.threePrimeUTR{};Parent={}\n".format(transID, idx, transID))
    o.close()


parser = argparse.ArgumentParser(description="Parse a GFF3 file lacking UTR entries and add them in. Note that this program assumes a reasonably formatted GFF3 file")
parser.add_argument("--fiveUTRname", default="five_prime_UTR", help="The label for 5' UTR entries (default: %(default)s)")
parser.add_argument("--threeUTRname", default="three_prime_UTR", help="The label for 3' UTR entries (default: %(default)s)")
parser.add_argument("--topLevelID", default="gene", help="The 'type' designating the top-level entry (default: %(default)s)")
parser.add_argument("--secondLevelID", default="mRNA", help="The 'type' designating the second-level entry, typically something like 'mRNA' or 'transcript' (default: %(default)s)")
parser.add_argument("--exonID", default="exon", help="The 'type' designating exons (default: %(default)s)")
parser.add_argument("--CDSID", default="CDS", help="The 'type' designating CDS (default: %(default)s)")
parser.add_argument("input", metavar="input.gff3", help="Input file name")
parser.add_argument("output", metavar="output.gff3", help="Output file name")
args = parser.parse_args()

genes, transcripts, g2t = parseGFF3(args.input, topLevel=args.topLevelID, secondLevel=args.secondLevelID, exon=args.exonID, CDS=args.CDSID)

# Add the UTRs
findUTRs(transcripts)
findUTRs(transcripts, fivePrime=False)

# Write the output
saveGFF(args.output, genes, transcripts, g2t, args.fiveUTRname, args.threeUTRname)
