#!/bin/bash nextflow

include { 
    dnaseq_trim_reads;
    dnaseq_filter_reads;
    dnaseq_count_reads;
    dnaseq_split_reads_to_chunks;
    dnaseq_map_barcodes;
    dnaseq_collapse_barcodes
} from "./modules/extract_dnaseq_barcodes"

include { 
    sc_get_unmapped_reads;
    sc_remove_low_qual_reads;
    sc_retain_reads_with_CB_tag;
    sc_split_unmapped_reads;
    sc_map_unmapped_reads;
    sc_merge_barcodes 
} from "./modules/extract_sc_clone_barcodes"

workflow {

    if (params.mode == 'DNAseq') {
        ch_barcode_chunks = Channel.fromPath("${params.dnaseq_fastq_files}/*.fastq.gz") | 
            dnaseq_trim_reads |
            dnaseq_filter_reads |
            dnaseq_count_reads |
            dnaseq_split_reads_to_chunks
        
        ch_barcode_mappings = dnaseq_map_barcodes(ch_barcode_chunks.flatten())
        dnaseq_collapse_barcodes(ch_barcode_mappings.collect())

    } 
    
    if (params.mode == 'scRNAseq') {
        ch_unmapped_fastas = Channel.fromPath("${params.scrnaseq_bam_files}/*.bam") | 
            sc_get_unmapped_reads |
            sc_remove_low_qual_reads |
            sc_retain_reads_with_CB_tag |
            sc_split_unmapped_reads
        
        ch_mapped_fastas = sc_map_unmapped_reads(ch_unmapped_fastas[0].flatten())
        sc_merge_barcodes(ch_mapped_fastas.collect())
    }
}