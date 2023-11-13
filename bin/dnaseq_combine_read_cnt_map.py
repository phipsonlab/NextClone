#!/usr/bin/env python

# coding: utf-8
# Combine the barcodes we found from flexiplex discovery and mapping, per fasta chunk.

import argparse
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

def combine_read_cnt_and_map(barcode_mapping_file, reads_fasta_file, output_filename):
    """
    Merges read counts and mapping information produced by flexiplex.

    Parameters
    ----------
    barcode_mapping_file : str
        File containing the mapping of reads to known barcodes as produced by flexiplex.
        Example format:
        Read    CellBarcode     FlankEditDist   BarcodeEditDist UMI
        READ0   AAATTGTGATTTTAGATAGT    0       0
        READ1   CCGTCGTATTGTCTGCCGCT    0       0

    reads_fasta_file : str
        FASTA file comprising the reads that correspond to the mapping in `barcode_mapping_file`.
        Example format:
        >READ0 count.23,sample.test
        STARTAAATTGTGATTTTAGATAGT
        >READ1 count.55,sample.test
        STARTCCGTCGTATTGTCTGCCGCT

    output_filename : str
        The name of the output file where the combined data will be saved.
    """

    
    barcode_mapping = pd.read_csv(barcode_mapping_file, sep='\t', index_col=False)

    
    fasta_chunks = []
    with open(reads_fasta_file) as in_fasta:
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

    barcode_counts_mapped.to_csv(output_filename, index=False, na_rep=pd.NA)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--unmapped_chunk", help="Barcode count output by Flexiplex discovery")
    parser.add_argument("--mapped_chunk", help="Barcode count output by Flexiplex map")
    parser.add_argument("--out_file", help="Name of the output file")

    args = parser.parse_args()

    combine_read_cnt_and_map(
        barcode_mapping_file=args.mapped_chunk,
        reads_fasta_file=args.unmapped_chunk,
        output_filename=args.out_file
    )

