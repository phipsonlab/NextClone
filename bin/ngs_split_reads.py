#!/usr/bin/env python

# coding: utf-8

import argparse
import pandas as pd

from math import ceil

# Use BioPython - yes even though we can just write it out using normal fwrite.
# Just so it's more conforming to standards.
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def split_reads(infile_name, sample_name, n_chunks, outdir):
    """
    Splits Flexiplex read counts into chunks.

    This function divides the read counts produced by Flexiplex (in Discovery mode) 
    into chunks. 
    The input read counts format is a tab-separated value containing the 
    read sequence followed by its occurrence count.

    Parameters
    ----------
    infile_name : str
        Path to the file containing read counts produced by Flexiplex.
    sample_name : str
        Identifier for the sample associated with the read counts. Used as a prefix for the chunk filenames.
    n_chunks : int
        Specifies the desired number of chunks to divide the reads into.
    outdir : str
        Directory where the generated chunk files will be saved.
    """
    
    df = pd.read_csv(infile_name, sep="\t", header=None)
    df.columns = ['barcode', 'count']

    chunk = 0
    outfile_handle = None

    # calculate how many reads to put in each chunk
    n_reads_per_chunk = ceil(df.shape[0] / n_chunks)

    for i, row in df.iterrows():
        if i % n_reads_per_chunk == 0:
            if outfile_handle is not None:
                outfile_handle.close()
            outfile_handle = open(f"{outdir}/{sample_name}_chunk{chunk}.fasta", "w")
            chunk += 1
            
        
        # Add the word START as the dummy adapter for the barcode for Flexiplex
        rec = SeqRecord(
            seq=Seq(f"START{row['barcode']}"), 
            id=f"READ{i}", 
            description=f"count.{row['count']},sample.{sample_name}")

        SeqIO.write(rec, outfile_handle, "fasta")
        
    outfile_handle.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--barcode_file", help="Barcode count output by Flexiplex discovery")
    parser.add_argument("--sample_name", help="Name of the sample")
    parser.add_argument("--n_chunks", help="Number chunks to divide the reads into", type=int)
    parser.add_argument("--outdir", help="Output directory")

    args = parser.parse_args()

    split_reads(
        infile_name=args.barcode_file, 
        sample_name=args.sample_name, 
        n_chunks=args.n_chunks,
        outdir=args.outdir
    )


