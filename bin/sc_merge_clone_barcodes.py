#!/usr/bin/env python

import pandas as pd
import sys

# Uncomment me for actual workflow
chunk_files = sys.argv[1: len(sys.argv)-1]

chunks_to_merge = []

for chunk in chunk_files:
    df = pd.read_csv(chunk, sep='\t', index_col=False)

    df.rename(columns={'CellBarcode': 'CloneBarcode'}, inplace=True)
    df['CellBarcode'] = [x.split('|')[0].split('_')[1] for x in df['Read'].to_numpy()]
    df['ReadId'] = [x.split('|')[1] for x in df['Read'].to_numpy()]

    # Don't think i need this column as it's a concatenation of cell barcode
    # and read id
    df = df.drop(columns = 'Read')

    # Have to identify which "sample" (essentially the input bam file) the barcode come from.
    # This is because the method is merging the clone barcodes from EVERY SINGLE INPUT BAM FILE.
    # We can do this by looking at the filename of the chunk.
    source_bam_file = chunk.split("_reads_unmapped")[0]
    df['SourceBAMFile'] = source_bam_file
    chunks_to_merge.append(df)


chunks_df = pd.concat(chunks_to_merge)

# Reorder columns
chunks_df = chunks_df[['CellBarcode', 'CloneBarcode', 'SourceBAMFile',
                      'ReadId', 'FlankEditDist', 'BarcodeEditDist','UMI']]

chunks_df.to_csv(sys.argv[len(sys.argv)-1], index=False)