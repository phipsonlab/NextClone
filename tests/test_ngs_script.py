import os
from Bio import SeqIO
import pandas as pd
import sys
sys.path.append('../bin')

from ngs_split_reads import split_reads
from ngs_combine_read_cnt_map import combine_read_cnt_and_map
from ngs_count_barcodes import count_barcodes

def test_split_reads_uneven_chunks(tmp_path):

    infile = "data/ngs/barcodes_counts.txt"
    sample_name = "test"
    n_chunks = 3
    # each chunk should only have at most 2 reads as we only have 5 reads in the data.
    exp_n_reads = [2,2,1]
    
    split_reads(
        infile_name=infile, 
        sample_name=sample_name, 
        n_chunks=n_chunks, 
        outdir=tmp_path
    )

    # read the barcode count and store the sequence in a list and the count
    # in another list
    # Initialize empty lists for the two columns
    sequences = []
    counts = []

    with open(infile, 'r') as f:
        for line in f:
            # Split each line by the tab character
            seq_value, cnt_value = line.strip().split('\t')
            
            sequences.append(seq_value)
            counts.append(int(cnt_value))

    # to keep track of the code that check the sequences and counts of each reads
    n_line = 0
    for i in range(n_chunks):
        # first check the chunk exists
        fasta_file = f"{tmp_path}/{sample_name}_chunk{i}.fasta"
        assert os.path.exists(fasta_file)

        with open(fasta_file, 'r') as f:

            # to check the number of reads in each chunk
            n_reads_in_file = 0

            # Check that the reads in each file is correct,
            # i.e., the same as those in the infile
            for record in SeqIO.parse(f, "fasta"):
                exp_desc = f"READ{n_line} count.{counts[n_line]},sample.{sample_name}"
                exp_sequence = f"START{sequences[n_line]}"
                assert record.description == exp_desc, f"Expected: {exp_desc}, but got {record.description}"
                assert str(record.seq) == exp_sequence, f"Expected: {exp_sequence}, but got {str(record.seq)}"
                n_line += 1
                n_reads_in_file += 1

            assert n_reads_in_file == exp_n_reads[i]

def test_split_reads_more_chunks_than_reads(tmp_path):

    infile = "data/ngs/barcodes_counts.txt"
    sample_name = "test"
    n_chunks = 6
    # each chunk should only have at most 2 reads as we only have 5 reads in the data.
    exp_n_reads = [1,1,1,1,1]
    
    split_reads(
        infile_name=infile, 
        sample_name=sample_name, 
        n_chunks=n_chunks, 
        outdir=tmp_path
    )

    # read the barcode count and store the sequence in a list and the count
    # in another list
    # Initialize empty lists for the two columns
    sequences = []
    counts = []

    with open(infile, 'r') as f:
        for line in f:
            # Split each line by the tab character
            seq_value, cnt_value = line.strip().split('\t')
            
            sequences.append(seq_value)
            counts.append(int(cnt_value))

    # to keep track of the code that check the sequences and counts of each reads
    n_line = 0
    # -1 because we have more chunks than reads
    for i in range(n_chunks - 1):
        # first check the chunk exists
        fasta_file = f"{tmp_path}/{sample_name}_chunk{i}.fasta"
        assert os.path.exists(fasta_file)

        with open(fasta_file, 'r') as f:

            # to check the number of reads in each chunk
            n_reads_in_file = 0

            # Check that the reads in each file is correct,
            # i.e., the same as those in the infile
            for record in SeqIO.parse(f, "fasta"):
                exp_desc = f"READ{n_line} count.{counts[n_line]},sample.{sample_name}"
                exp_sequence = f"START{sequences[n_line]}"
                assert record.description == exp_desc, f"Expected: {exp_desc}, but got {record.description}"
                assert str(record.seq) == exp_sequence, f"Expected: {exp_sequence}, but got {str(record.seq)}"
                n_line += 1
                n_reads_in_file += 1

            assert n_reads_in_file == exp_n_reads[i]

    assert not os.path.exists(f"{tmp_path}/{sample_name}_chunk5.fasta")

def test_split_reads_onechunk(tmp_path):

    infile = "data/ngs/barcodes_counts.txt"
    sample_name = "test"
    n_chunks = 1
    # the chunk should have all the 5 reads
    exp_n_reads = [5]
    
    split_reads(
        infile_name=infile, 
        sample_name=sample_name, 
        n_chunks=n_chunks, 
        outdir=tmp_path
    )

    # read the barcode count and store the sequence in a list and the count
    # in another list
    # Initialize empty lists for the two columns
    sequences = []
    counts = []

    with open(infile, 'r') as f:
        for line in f:
            # Split each line by the tab character
            seq_value, cnt_value = line.strip().split('\t')
            
            sequences.append(seq_value)
            counts.append(int(cnt_value))

    # to keep track of the code that check the sequences and counts of each reads
    n_line = 0
    for i in range(n_chunks):
        # first check the chunk exists
        fasta_file = f"{tmp_path}/{sample_name}_chunk{i}.fasta"
        assert os.path.exists(fasta_file)

        with open(fasta_file, 'r') as f:

            # to check the number of reads in each chunk
            n_reads_in_file = 0

            # Check that the reads in each file is correct,
            # i.e., the same as those in the infile
            for record in SeqIO.parse(f, "fasta"):
                exp_desc = f"READ{n_line} count.{counts[n_line]},sample.{sample_name}"
                exp_sequence = f"START{sequences[n_line]}"
                assert record.description == exp_desc, f"Expected: {exp_desc}, but got {record.description}"
                assert str(record.seq) == exp_sequence, f"Expected: {exp_sequence}, but got {str(record.seq)}"
                n_line += 1
                n_reads_in_file += 1

            assert n_reads_in_file == exp_n_reads[i]

def test_combine_read_cnt_map(tmp_path):
    barcode_mapping_file = "data/ngs/unmapped_chunk_1_reads_barcodes.txt"
    reads_fasta_file = "data/ngs/unmapped_chunk_1.fasta"
    output_file_path = f"{tmp_path}/test_out.csv"

    combine_read_cnt_and_map(barcode_mapping_file, reads_fasta_file, output_file_path)

    df = pd.read_csv(output_file_path)
    assert df.shape[0] == 10
    assert df.shape[1] == 9

    exp_read_ids = [f"READ{x}" for x in [0,1,2,3,4,495,496,497,498,499]]
    exp_read_cnts = [24186,6989,6083,5549,4832,1325,1324,1323,1323,1322]
    exp_read_barcode = [
        'AGGGGAGTCGCGTGGTAGGC',
        'TGTCTAATGGGGGTGTCACT',
        'GCGGAAACGTACTCAATAAG',
        'GAGGATATTTTGGTTTCTAG',
        'AGAGGTCTTTAGTTCAACTT',
        'GTTGACAGTATTTGCCGGCA',
        'TGGACGCGTACTCTCCAGTC',
        'ATGGCTCCCTTGGGGAATGA',
        'CCTGTTGATCGTGCGCGCAG',
        'GATTGTATCGGATTGGGGGA'
    ]
    exp_known_barcode = [
        'AGGGGAGTCGCGTGGTAGGC',
        'TGTCTAATGGGGGTGTCACT',
        pd.NA,
        'GAGGATATTTTGGTTTCTAG',
        'AGAGGTCTTTAGTTCAACTT',
        'GTTGACAGTATTTGCCGGCA',
        'TGGACGCGTACTCTCCAGTC',
        'ATGGCTCCCTTGGGGAATGA',
        'CCTGTTGATCGTGCGCGCAG',
        'GATTGTATCGGATTGGGGGA'
    ]
    for i in range(len(exp_read_ids)):
        assert df[df['read_id'] == exp_read_ids[i]]['read_count'].iloc[0] == exp_read_cnts[i]
        assert df[df['read_id'] == exp_read_ids[i]]['read_sequence'].iloc[0] == exp_read_barcode[i]

        # No known barcode
        if i == 2:
            assert pd.isna(df[df['read_id'] == exp_read_ids[i]]['known_barcode'].iloc[0])
        else:
            assert df[df['read_id'] == exp_read_ids[i]]['known_barcode'].iloc[0] == exp_known_barcode[i]

def test_count_barcodes_2chunks(tmp_path):
    count_barcodes(tmp_path, ["data/ngs/barcodes_counts_mapped_chunk1.csv", "data/ngs/barcodes_counts_mapped_chunk2.csv"])
    
    df = pd.read_csv(f"{tmp_path}/clone_barcode_counts.csv")
    
    exp_samples = ['test1', 'test2']
    exp_known_barcode = ['AGGGGAGTCGCGTGGTAGGC', 'AGGGGAGTCGCGTGGTAGGC']
    exp_count = [2, 10]

    for i in range(len(exp_samples)):
        assert df[df['sample_name'] == exp_samples[i]]['known_barcode'].iloc[0] == exp_known_barcode[i]
        assert df[df['sample_name'] == exp_samples[i]]['read_count'].iloc[0] == exp_count[i]

def test_count_barcodes_1chunk(tmp_path):
    count_barcodes(tmp_path, ["data/ngs/barcodes_counts_mapped_chunk1.csv"])
    
    df = pd.read_csv(f"{tmp_path}/clone_barcode_counts.csv")

    assert df[df['sample_name'] == 'test1']['known_barcode'].iloc[0] == 'AGGGGAGTCGCGTGGTAGGC'
    assert df[df['sample_name'] == 'test1']['read_count'].iloc[0] == 2
