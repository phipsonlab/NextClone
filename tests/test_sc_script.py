import os
import pysam
import numpy as np
import sys
sys.path.append('../bin')

from sc_remove_low_qual_reads import remove_low_qual_reads
from sc_split_reads import split_reads

def test_low_qual_reads_removed(tmp_path):

    input_bam = "data/sc/mix_qual_reads.bam"
    out_bam = f"{tmp_path}/high_qual_reads.bam"

    remove_low_qual_reads(input_bam, out_bam, 30)

    pysam.index(out_bam)

    actual_read_ids = []
    expected_read_ids = [
        "A00121:674:HVMHWDSX3:2:2201:30056:3756",
        "A00121:674:HVMHWDSX3:1:2327:9516:3443",
        "A00121:674:HVMHWDSX3:1:2130:30996:21042",
        "A00121:674:HVMHWDSX3:2:1246:24035:12587",
        "A00121:674:HVMHWDSX3:2:2625:29279:9267",
        "A00121:674:HVMHWDSX3:2:1343:10945:6480",
        "A00121:674:HVMHWDSX3:2:2653:30156:30232",
        "A00121:674:HVMHWDSX3:2:2568:13711:24330",
        "A00121:674:HVMHWDSX3:1:2165:30337:5556",
        "A00121:674:HVMHWDSX3:2:2308:30517:11350",
        "A00121:674:HVMHWDSX3:2:1572:7111:20415"
    ]
    bamfile = pysam.Samfile(out_bam, "rb")
    for row in bamfile:
        actual_read_ids.append(row.query_name)
    bamfile.close()

    assert len(np.setdiff1d(actual_read_ids, expected_read_ids)) == 0
    assert len(np.setdiff1d(expected_read_ids, actual_read_ids)) == 0

def test_split_reads_start(tmp_path):
    
    input_bam = "data/sc/unmapped_reads.bam"
    n_chunks = 4
    # 25 reads, 4 chunks
    exp_n_reads = [7,7,7,4]

    split_reads(
        input_bam_filename=input_bam,
        outdir=tmp_path,
        n_chunks=n_chunks
    )

    expected_splitted_read_ids = [
        "A00121:674:HVMHWDSX3:2:2462:6705:33207", 
        "A00121:674:HVMHWDSX3:1:1550:9091:25974", 
        "A00121:674:HVMHWDSX3:2:2667:10827:26287", 
        "A00121:674:HVMHWDSX3:2:1529:6063:21637", 
        "A00121:674:HVMHWDSX3:2:1643:29261:30592", 
        "A00121:674:HVMHWDSX3:1:2546:16857:10426", 
        "A00121:674:HVMHWDSX3:1:1277:22815:33082", 
        "A00121:674:HVMHWDSX3:1:1577:3956:6872", 
        "A00121:674:HVMHWDSX3:1:2666:5312:29543", 
        "A00121:674:HVMHWDSX3:1:1403:2935:14935", 
        "A00121:674:HVMHWDSX3:2:2161:27923:26929", 
        "A00121:674:HVMHWDSX3:2:1606:27371:32800", 
        "A00121:674:HVMHWDSX3:2:2259:28320:13369", 
        "A00121:674:HVMHWDSX3:1:2270:4607:11569", 
        "A00121:674:HVMHWDSX3:2:2354:16152:8672", 
        "A00121:674:HVMHWDSX3:2:2377:9390:16376", 
        "A00121:674:HVMHWDSX3:2:2516:11089:10394", 
        "A00121:674:HVMHWDSX3:1:2531:18439:16172", 
        "A00121:674:HVMHWDSX3:1:1265:15203:3615", 
        "A00121:674:HVMHWDSX3:2:2566:27407:4147", 
        "A00121:674:HVMHWDSX3:1:2152:18231:1908", 
        "A00121:674:HVMHWDSX3:1:1248:16089:33270", 
        "A00121:674:HVMHWDSX3:1:1355:23104:8500", 
        "A00121:674:HVMHWDSX3:1:1134:25536:14152", 
        "A00121:674:HVMHWDSX3:1:1540:28628:31093"
    ]

    # For checking the content of the fasta file
    actual_splitted_read_ids = []
    
    for i in range(n_chunks):
        # Check first that the file exists
        assert os.path.exists(f"{tmp_path}/unmapped_reads_unmapped_chunk_{i}.fasta")

        with pysam.FastxFile(f"{tmp_path}/unmapped_reads_unmapped_chunk_{i}.fasta") as fasta:
            
            n_reads_in_file = 0

            for read in fasta:
                actual_splitted_read_ids.append(read.name.split("|")[1])
                n_reads_in_file += 1
            
            assert n_reads_in_file == exp_n_reads[i]

    assert len(np.setdiff1d(actual_splitted_read_ids, expected_splitted_read_ids)) == 0
    assert len(np.setdiff1d(expected_splitted_read_ids, actual_splitted_read_ids)) == 0

def test_split_reads_more_chunks_than_reads(tmp_path):
    
    input_bam = "data/sc/unmapped_reads.bam"
    n_chunks = 27
    n_reads_in_bam = 25
    exp_n_reads = [1] * n_reads_in_bam

    split_reads(
        input_bam_filename=input_bam,
        outdir=tmp_path,
        n_chunks=n_chunks
    )

    expected_splitted_read_ids = [
        "A00121:674:HVMHWDSX3:2:2462:6705:33207", 
        "A00121:674:HVMHWDSX3:1:1550:9091:25974", 
        "A00121:674:HVMHWDSX3:2:2667:10827:26287", 
        "A00121:674:HVMHWDSX3:2:1529:6063:21637", 
        "A00121:674:HVMHWDSX3:2:1643:29261:30592", 
        "A00121:674:HVMHWDSX3:1:2546:16857:10426", 
        "A00121:674:HVMHWDSX3:1:1277:22815:33082", 
        "A00121:674:HVMHWDSX3:1:1577:3956:6872", 
        "A00121:674:HVMHWDSX3:1:2666:5312:29543", 
        "A00121:674:HVMHWDSX3:1:1403:2935:14935", 
        "A00121:674:HVMHWDSX3:2:2161:27923:26929", 
        "A00121:674:HVMHWDSX3:2:1606:27371:32800", 
        "A00121:674:HVMHWDSX3:2:2259:28320:13369", 
        "A00121:674:HVMHWDSX3:1:2270:4607:11569", 
        "A00121:674:HVMHWDSX3:2:2354:16152:8672", 
        "A00121:674:HVMHWDSX3:2:2377:9390:16376", 
        "A00121:674:HVMHWDSX3:2:2516:11089:10394", 
        "A00121:674:HVMHWDSX3:1:2531:18439:16172", 
        "A00121:674:HVMHWDSX3:1:1265:15203:3615", 
        "A00121:674:HVMHWDSX3:2:2566:27407:4147", 
        "A00121:674:HVMHWDSX3:1:2152:18231:1908", 
        "A00121:674:HVMHWDSX3:1:1248:16089:33270", 
        "A00121:674:HVMHWDSX3:1:1355:23104:8500", 
        "A00121:674:HVMHWDSX3:1:1134:25536:14152", 
        "A00121:674:HVMHWDSX3:1:1540:28628:31093"
    ]

    # For checking the content of the fasta file
    actual_splitted_read_ids = []
    
    for i in range(n_reads_in_bam):
        # Check first that the file exists
        assert os.path.exists(f"{tmp_path}/unmapped_reads_unmapped_chunk_{i}.fasta")

        with pysam.FastxFile(f"{tmp_path}/unmapped_reads_unmapped_chunk_{i}.fasta") as fasta:
            
            n_reads_in_file = 0

            for read in fasta:
                actual_splitted_read_ids.append(read.name.split("|")[1])
                n_reads_in_file += 1
            
            assert n_reads_in_file == exp_n_reads[i]

    assert not os.path.exists(f"{tmp_path}/unmapped_reads_unmapped_chunk_25.fasta")
    assert not os.path.exists(f"{tmp_path}/unmapped_reads_unmapped_chunk_26.fasta")

    assert len(np.setdiff1d(actual_splitted_read_ids, expected_splitted_read_ids)) == 0
    assert len(np.setdiff1d(expected_splitted_read_ids, actual_splitted_read_ids)) == 0
