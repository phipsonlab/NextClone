#!/usr/bin/env python

# coding: utf-8
# Break the barcode into chunk of 100 reads

import argparse
import pandas as pd

# Use BioPython - yes even though we can just write it out using normal fwrite.
# Just so it's more conforming to standards.
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("--barcode_file", help="Barcode count output by Flexiplex discovery")
parser.add_argument("--sample_name", help="Name of the sample")
parser.add_argument("--nreads_per_chunk", help="Number of reads per chunk/Fasta file", type=int)

args = parser.parse_args()
infile_name = args.barcode_file
sample_name = args.sample_name
n_reads_per_chunk = args.nreads_per_chunk

df = pd.read_csv(infile_name, sep="\t", header=None)
df.columns = ['barcode', 'count']

chunk = 0
outfile_handle = None

for i, row in df.iterrows():
    if i % n_reads_per_chunk == 0:
        if outfile_handle is not None:
            outfile_handle.close()
        chunk += 1
        outfile_handle = open(f"{sample_name}_chunk{chunk}.fasta", "w")
    
    rec = SeqRecord(
        seq=Seq(f"START{row['barcode']}"), 
        id=f"READ{i}", 
        description=f"count.{row['count']},sample.{sample_name}")

    SeqIO.write(rec, outfile_handle, "fasta")
    
outfile_handle.close()
