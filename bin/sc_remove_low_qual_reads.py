#!/usr/bin/env python

# Remove reads with phred less than 20
# Have to make this as there doesn't seem to be anything around that allows me to do this!
# Write out as BAM file to save space
import sys
import pysam
import numpy as np

pysam.index(sys.argv[1])

phred_threshold = int(sys.argv[2])
bamfile = pysam.Samfile(sys.argv[1], "rb")
out_bam = pysam.Samfile(sys.argv[3], "wb", template=bamfile)

for row in bamfile:
    phred_mean = np.mean(row.query_qualities)
    if phred_mean >= phred_threshold:
        out_bam.write(row)

bamfile.close()
out_bam.close()
