#!/usr/bin/env python
import sys
import pysam
from Bio import SeqIO, Seq, SeqRecord

pysam.index(sys.argv[1])

# See: https://www.biostars.org/p/6970/
bamfile = pysam.AlignmentFile(sys.argv[1], "rb", check_sq=False)

chunk = 0
fout = None
n_row = -1

outdir = sys.argv[2]

n_reads_per_chunk = int(sys.argv[3])

for row in bamfile.fetch(until_eof=True):
    n_row += 1

    # have to check if the tag exists first as some reads appear to be missing CB
    if not row.has_tag("CB"):
        with open("reads_missing_cb.txt", "a") as mcb:
            mcb.write(row.query_name)
            mcb.write("\n")
            continue

    # open a new fasta file to write
    if n_row % n_reads_per_chunk == 0:
        if fout is not None:
            fout.close()
        chunk += 1
        fout = open(f"{outdir}/unmapped_chunk_{chunk}.fasta", "w")

    # Attach the query name as well 
    cell_id = f"Cell_{row.get_tag('CB')}|{row.query_name}"
    seq = Seq.Seq(row.query_sequence)

    rec = SeqRecord.SeqRecord(seq, cell_id, "", "")

    SeqIO.write(rec, fout, "fasta")
        
if fout is not None:
    fout.close()
bamfile.close()