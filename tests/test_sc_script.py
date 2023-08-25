import os
import pysam
import numpy as np
import sys
sys.path.append('../bin')

from sc_remove_low_qual_reads import remove_low_qual_reads
from sc_split_reads import split_reads

def test_low_qual_reads_removed(tmp_path):

    input_bam = "data/mix_qual_reads.bam"
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

def test_split_reads_start_withCB(tmp_path):
    input_bam = "data/mix_with_without_cb_start_withCB.bam"
    reads_missing_cb = f"{tmp_path}/mix_with_without_cb_start_withCB_missing_cb.txt"

    split_reads(
        input_bam_filename=input_bam,
        outdir=tmp_path,
        n_reads_per_chunk=3,
        reads_missing_cb_file=reads_missing_cb
    )

    expected_splitted_read_ids = [
        "A00121:674:HVMHWDSX3:1:2537:28664:28808",
        "A00121:674:HVMHWDSX3:2:2440:12861:33254",
        "A00121:674:HVMHWDSX3:1:2532:27389:27289",
        "A00121:674:HVMHWDSX3:2:2367:30689:29528",
        "A00121:674:HVMHWDSX3:2:2563:22815:36714",
        "A00121:674:HVMHWDSX3:2:2642:31575:32221",
        "A00121:674:HVMHWDSX3:2:2175:7500:13072",
        "A00121:674:HVMHWDSX3:1:2449:16034:8155",
        "A00121:674:HVMHWDSX3:1:1675:16721:4366",
        "A00121:674:HVMHWDSX3:2:2638:15130:6057",
        "A00121:674:HVMHWDSX3:1:1247:17291:26428"
    ]

    expected_missingCB_read_ids = [
        "A00121:674:HVMHWDSX3:2:2335:19334:1031",
        "A00121:674:HVMHWDSX3:1:2320:18322:36980",
        "A00121:674:HVMHWDSX3:2:1325:14036:5259",
        "A00121:674:HVMHWDSX3:2:2557:27190:7874",
        "A00121:674:HVMHWDSX3:2:2556:29053:22811",
        "A00121:674:HVMHWDSX3:2:2201:7319:16172",
        "A00121:674:HVMHWDSX3:2:1376:24985:28040",
        "A00121:674:HVMHWDSX3:2:2375:4327:36667",
        "A00121:674:HVMHWDSX3:1:1348:20076:26303",
        "A00121:674:HVMHWDSX3:1:1348:20048:26318",
        "A00121:674:HVMHWDSX3:1:2247:1750:19460"
    ]

    n_expected_chunks = 4

    for i in range(n_expected_chunks):
        assert os.path.exists(f"{tmp_path}/mix_with_without_cb_start_withCB_unmapped_chunk_{i}.fasta")
    
    # Check the content of the fasta file
    actual_splitted_read_ids = []
    for i in range(n_expected_chunks):
        with pysam.FastxFile(f"{tmp_path}/mix_with_without_cb_start_withCB_unmapped_chunk_{i}.fasta") as fasta:
            for read in fasta:
                actual_splitted_read_ids.append(
                    read.name.split("|")[1]
                )
    assert len(np.setdiff1d(actual_splitted_read_ids, expected_splitted_read_ids)) == 0
    assert len(np.setdiff1d(expected_splitted_read_ids, actual_splitted_read_ids)) == 0

    # Check the content of reads missing cb
    with open(reads_missing_cb, 'r') as file:
        lines = [line.strip() for line in file]
        assert len(np.setdiff1d(lines, expected_missingCB_read_ids)) == 0
        assert len(np.setdiff1d(expected_missingCB_read_ids, lines)) == 0

def test_split_reads_start_withoutCB(tmp_path):
    input_bam = "data/mix_with_without_cb_start_withoutCB.bam"
    reads_missing_cb = f"{tmp_path}/mix_with_without_cb_start_withoutCB_missing_cb.txt"

    split_reads(
        input_bam_filename=input_bam,
        outdir=tmp_path,
        n_reads_per_chunk=3,
        reads_missing_cb_file=reads_missing_cb
    )

    expected_splitted_read_ids = [
        "A00121:674:HVMHWDSX3:2:1417:1678:36620",
        "A00121:674:HVMHWDSX3:2:1417:1651:36135",
        "A00121:674:HVMHWDSX3:2:2128:20057:9298",
        "A00121:674:HVMHWDSX3:1:1631:30355:27352",
        "A00121:674:HVMHWDSX3:1:1616:18340:32471",
        "A00121:674:HVMHWDSX3:2:2239:6361:9251",
        "A00121:674:HVMHWDSX3:2:2239:5719:18756",
        "A00121:674:HVMHWDSX3:1:1227:28094:9001",
        "A00121:674:HVMHWDSX3:1:2126:28700:36573",
        "A00121:674:HVMHWDSX3:2:2402:5755:19664",
        "A00121:674:HVMHWDSX3:2:1152:19804:13902"
    ]

    expected_missingCB_read_ids = [
        "A00121:674:HVMHWDSX3:2:2335:19334:1031",
        "A00121:674:HVMHWDSX3:1:2320:18322:36980"
    ]

    n_expected_chunks = 4

    for i in range(n_expected_chunks):
        assert os.path.exists(f"{tmp_path}/mix_with_without_cb_start_withoutCB_unmapped_chunk_{i}.fasta")
    
    # Check the content of the fasta file
    actual_splitted_read_ids = []
    for i in range(n_expected_chunks):
        with pysam.FastxFile(f"{tmp_path}/mix_with_without_cb_start_withoutCB_unmapped_chunk_{i}.fasta") as fasta:
            for read in fasta:
                actual_splitted_read_ids.append(
                    read.name.split("|")[1]
                )
    assert len(np.setdiff1d(actual_splitted_read_ids, expected_splitted_read_ids)) == 0
    assert len(np.setdiff1d(expected_splitted_read_ids, actual_splitted_read_ids)) == 0

    # Check the content of reads missing cb
    with open(reads_missing_cb, 'r') as file:
        lines = [line.strip() for line in file]
        assert len(np.setdiff1d(lines, expected_missingCB_read_ids)) == 0
        assert len(np.setdiff1d(expected_missingCB_read_ids, lines)) == 0

    
    
     
    