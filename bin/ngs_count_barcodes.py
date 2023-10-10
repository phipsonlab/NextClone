#!/usr/bin/env python
# coding: utf-8

# Collapse barcode count by sample and known barcode
import sys
import pandas as pd

def count_barcodes(outdir, chunk_files):
    
    if len(chunk_files) == 1:
        df = pd.read_csv(chunk_files[0])
    else:
        df = pd.concat((pd.read_csv(x) for x in chunk_files), ignore_index=True)
    barcode_count = df[['known_barcode','read_count', 'sample_name']].groupby(['sample_name','known_barcode']).sum()
    barcode_count.to_csv(f"{outdir}/clone_barcode_counts.csv")

if __name__ == '__main__':
    chunk_files = sys.argv[2: len(sys.argv)]
    count_barcodes(sys.argv[1], chunk_files)