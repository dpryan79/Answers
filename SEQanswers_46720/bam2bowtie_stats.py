#!/usr/bin/env python
import pysam
import sys

if(len(sys.argv) != 2) :
    print("Usage: %s file.bam" % sys.argv[0]);
    print("This program reproduces the metrics printed by bowtie2 after it completes its alignments.")
    print("Only single end metrics are supported at the moment.")
    sys.exit()

f = pysam.Samfile(sys.argv[1], 'rb')

total = 0
unpaired = 0
mapped = 0
multimap = 0
for read in f.fetch():
    total += 1
    if(read.flag & 1 == 0) :
        unpaired += 1
    if(read.flag & 4 == 0) :
        mapped += 1
        if(read.mapq < 2) :
            AS = None
            XS = None
            for tag in read.tags :
                if(tag[0] == 'AS') :
                    AS = int(tag[1])
                elif(tag[0] == 'XS') :
                    XS = int(tag[1])
                if(AS != None and XS != None) :
                    break
            if(XS != None and XS == AS) :
                multimap += 1

print("%i reads; of these:" % total)
print("%i (%6.2f%%) were unpaired; of these:" % (unpaired, 100*unpaired/total))
print("%i (%6.2f%%) aligned 0 times" % (total-mapped, 100*(1-mapped/total)))
print("%i (%6.2f%%) aligned exactly 1 time" % (mapped-multimap, 100*(1-multimap/mapped)))
print("%i (%6.2f%%) aligned > 1 times" % (multimap, 100*multimap/mapped))
print("%6.2f%% overall alignment rate" % (100*mapped/total))
f.close()
