#!/usr/bin/env python
import sys
import pysam
from Bio import SeqIO, Seq, SeqRecord
import os

pysam.index(sys.argv[1])

# See: https://www.biostars.org/p/6970/
bamfile = pysam.Samfile(sys.argv[1], "rb")

bamfilename = os.path.splitext(os.path.basename(sys.argv[1]))[0]

outdir = sys.argv[2]
n_reads_per_chunk = int(sys.argv[3])

# for storing reads that have no CB
reads_missing_cb_file = sys.argv[4]

next_chunk_id = 0
fout = None
n_row_written = 0

for row in bamfile:

    # have to check if the tag exists first as some reads appear to be missing CB
    if not row.has_tag("CB"):
        # If the conde gets here, n_row_written, won't be updated
        with open(reads_missing_cb_file, "a") as mcb:
            mcb.write(row.query_name)
            mcb.write("\n")
            continue

    # open a new fasta file to write
    if n_row_written % n_reads_per_chunk == 0:
        if fout is not None:
            fout.close()
        fout = open(f"{outdir}/{bamfilename}_unmapped_chunk_{next_chunk_id}.fasta", "w")
        next_chunk_id += 1

    # Attach the query name as well 
    cell_id = f"Cell_{row.get_tag('CB')}|{row.query_name}"
    seq = Seq.Seq(row.query_sequence)

    rec = SeqRecord.SeqRecord(seq, cell_id, "", "")

    SeqIO.write(rec, fout, "fasta")

    n_row_written += 1
        
if fout is not None:
    fout.close()
bamfile.close()