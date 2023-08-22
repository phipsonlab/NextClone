#!/usr/bin/env python

# coding: utf-8
# Combine the barcodes we found from flexiplex discovery and mapping, per fasta chunk.

import argparse
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

parser = argparse.ArgumentParser()

parser.add_argument("--unmapped_chunk", help="Barcode count output by Flexiplex discovery")
parser.add_argument("--mapped_chunk", help="Barcode count output by Flexiplex map")
parser.add_argument("--out_file", help="Name of the output file")

args = parser.parse_args()
mapped_barcodes_file = args.mapped_chunk
chunk_file = args.unmapped_chunk
outfile = args.out_file

# Read in the barcode chunk that is mapped by Flexiplex
# It looks something like the following:
# Read    CellBarcode     FlankEditDist   BarcodeEditDist UMI
# READ0   AAAAAGTGATTCTAGACAGT    0       0
# READ1   CCGGCGTACTGGCAGCCGCA    0       0
barcode_mapping = pd.read_csv(mapped_barcodes_file, sep='\t', index_col=False)

# Read in the barcode chunk found by Flexiplex discovery
# It looks like the following:
# >READ0 count.17791,sample.vexGFP-1k-Palbociclib-R5
# STARTAAAAAGTGATTCTAGACAGT
# >READ1 count.15311,sample.vexGFP-1k-Palbociclib-R5
# STARTCCGGCGTACTGGCAGCCGCA
fasta_chunks = []
with open(chunk_file) as in_fasta:
    for title, seq in SimpleFastaParser(in_fasta):
        read_title = title.split(" ")
        read_id = read_title[0]

        # Get the count and sample name
        description = read_title[1].split(",")
        read_count = int(description[0].split(".")[1])
        sample_name = description[1].split(".")[1]
        
        # Get the actual read's sequence
        read = str(seq)[5:]
        
        a_read = [read_id, read_count, sample_name, read]
        fasta_chunks.append(a_read)
        
fasta_chunks_df = pd.DataFrame(fasta_chunks, columns=['read_id','read_count','sample_name', 'read_sequence'])

# Left join so if anything is not mapped, we will capture it.
barcode_counts_mapped = fasta_chunks_df.merge(barcode_mapping, how='left', left_on='read_id', right_on='Read')
barcode_counts_mapped.rename(columns={"CellBarcode": "known_barcode"}, inplace=True)

barcode_counts_mapped.to_csv(outfile, index=False, na_rep=pd.NA)