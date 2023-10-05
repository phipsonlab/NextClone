#!/usr/bin/env python

import sys
import argparse
import pysam
import os

from Bio import SeqIO, Seq, SeqRecord
from math import ceil

def split_reads(input_bam_filename, outdir, n_chunks):
    """
    Splits unmapped reads into chunks.

    Parameters
    ----------
    input_bam_filename : str
        Name of the BAM file containing the unmapped reads
    outdir : str
        Directory where the generated chunk files will be saved.
    n_chunks : int
        Specifies the desired number of chunks to divide the reads into.
    """

    pysam.index(input_bam_filename)

    # See: https://www.biostars.org/p/6970/
    input_bamfile = pysam.Samfile(input_bam_filename, "rb")

    bamfilename = os.path.splitext(os.path.basename(input_bam_filename))[0]

    # mapped and unmapped functions will return the total number of reads
    n_reads_per_chunk = ceil(input_bamfile.unmapped / n_chunks)
    
    next_chunk_id = 0
    fout = None
    n_row_written = 0

    for row in input_bamfile:
        # open a new fasta file to write
        if n_row_written % n_reads_per_chunk == 0:
            if fout is not None:
                fout.close()
            fout = open(f"{outdir}/{bamfilename}_unmapped_chunk_{next_chunk_id}.fasta", "w")
            next_chunk_id += 1

        # Attach the query name as well 
        cell_id = f"Cell_{row.get_tag('CB')}|{row.query_name}"
        seq = Seq.Seq(row.query_sequence)

        rec = SeqRecord.SeqRecord(seq, cell_id, "", "")

        SeqIO.write(rec, fout, "fasta")

        n_row_written += 1
            
    if fout is not None:
        fout.close()
    
    input_bamfile.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_bam_filename", help="BAM file containing the unmapped reads")
    parser.add_argument("--n_chunks", help="Number chunks to divide the reads into", type=int)
    parser.add_argument("--outdir", help="Output directory")

    args = parser.parse_args()

    # TODO remove me
    # input_bam_filename = sys.argv[1]
    # outdir = sys.argv[2]
    # n_reads_per_chunk = int(sys.argv[3])
    # reads_missing_cb_file = sys.argv[4]

    split_reads(
        input_bam_filename=args.input_bam_filename,
        outdir=args.outdir,
        n_chunks=args.n_chunks
    )