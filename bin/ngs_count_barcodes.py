#!/usr/bin/env python
# coding: utf-8

# Collapse barcode count by sample and known barcode
import sys
import pandas as pd

chunk_files = sys.argv[1: len(sys.argv)-1]

df = pd.concat((pd.read_csv(x) for x in chunk_files), ignore_index=True)
barcode_count = df[['known_barcode','read_count', 'sample_name']].groupby(['sample_name','known_barcode']).sum()
barcode_count.to_csv("clone_barcode_counts.csv")
