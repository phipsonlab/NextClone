#!/usr/bin/env python

import pandas as pd
import sys

# Uncomment me for actual workflow
chunk_files = sys.argv[1: len(sys.argv)]

# For local testing only
# files_loc = '/stornext/Bioinf/data/lab_phipson/givanna/senescence_cells_project/pilot_dataset/nextflow_pipeline/nextflow_int_files/mapped_chunks'
# chunk_files = [
#     f'{files_loc}/unmapped_chunk_1_reads_barcodes.txt',
#     f'{files_loc}/unmapped_chunk_2_reads_barcodes.txt']

chunks_to_merge = []

for chunk in chunk_files:
    df = pd.read_csv(chunk, sep='\t', index_col=False)
    df.rename(columns={'CellBarcode': 'CloneBarcode'}, inplace=True)
    df['CellBarcode'] = [x.split('|')[0].split('_')[1] for x in df['Read'].to_numpy()]
    df['ReadId'] = [x.split('|')[1] for x in df['Read'].to_numpy()]
    chunks_to_merge.append(df)


chunks_df = pd.concat(chunks_to_merge)

chunks_df.to_csv("clone_barcodes.csv", index=False)