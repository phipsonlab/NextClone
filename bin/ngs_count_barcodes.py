#!/usr/bin/env python
# coding: utf-8

# Collapse barcode count by sample and known barcode
import sys
import pandas as pd

from os import listdir
from os.path import isfile, join

basedir = sys.argv[1]
chunks = [join(basedir, f) for f in listdir(basedir) if isfile(join(basedir, f))]

df = pd.concat((pd.read_csv(x) for x in chunks), ignore_index=True)
barcode_count = df[['known_barcode','read_count', 'sample_name']].groupby(['sample_name','known_barcode']).sum()
barcode_count.to_csv("clone_barcode_counts.csv")
