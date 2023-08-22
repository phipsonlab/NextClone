#!/usr/bin/env python

# Extract cell barcode, umi, and read name from sam file
import pysam
import pandas as pd
import sys

pysam.index(sys.argv[1])

bamfile = pysam.AlignmentFile(sys.argv[1], "rb", check_sq=False)
reads = []
for row in bamfile.fetch(until_eof=True):
    if not row.has_tag("CB") or not row.has_tag("UB"):
        continue
    cell_barcode = row.get_tag("CB")
    umi = row.get_tag("UB")
    seq_name = row.query_name
    read = [cell_barcode, umi, seq_name]
    reads.append(read)
    

df = pd.DataFrame(reads, columns=['CellBarcode', 'UMI', 'ReadName'])
df.to_csv("unmapped_reads_cb_umi_readname.csv", index=False)
