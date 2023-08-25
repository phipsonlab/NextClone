#!/usr/bin/env python

# Remove reads with phred less than 20
# Have to make this as there doesn't seem to be anything around that allows me to do this!
# Write out as BAM file to save space
import sys
import pysam
import numpy as np

def remove_low_qual_reads(input_bam_filename, out_bam_filename, phred_threshold):
    pysam.index(input_bam_filename)

    bamfile = pysam.Samfile(input_bam_filename, "rb")
    out_bam = pysam.Samfile(out_bam_filename, "wb", template=bamfile)

    for row in bamfile:
        phred_mean = np.mean(row.query_qualities)
        if phred_mean >= phred_threshold:
            out_bam.write(row)

    bamfile.close()
    out_bam.close()

if __name__ == '__main__':
    input_bam_filename = sys.argv[1]
    phred_threshold = int(sys.argv[2])
    out_bam_filename = sys.argv[3]

    remove_low_qual_reads(input_bam_filename, out_bam_filename, phred_threshold)


