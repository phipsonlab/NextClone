#!/usr/bin/env nextflow

process sc_get_unmapped_reads {
    module 'samtools'
    publishDir "$params.publish_dir",  mode: 'copy'

    input:
        path possorted_aligned_reads_bam_file

    output:
        path 'unmapped_reads.bam'

    """
    samtools view -b -f 4 ${possorted_aligned_reads_bam_file} > unmapped_reads.bam
    """
}

process sc_remove_low_qual_reads {
    cpus 4
    memory '24 GB'
    time '1 hours'
    conda "${projectDir}/conda_env/extract_sc_env.yaml"
    publishDir "$params.publish_dir",  mode: 'copy'

    input:
        path unmapped_bam

    output:
        path "unmapped_filtered.bam"

    """
    sc_remove_low_qual_reads.py $unmapped_bam $params.phred_thres
    """
}

process sc_split_unmapped_reads {
    cpus 4
    memory '24 GB'
    time '3 hours'
    conda "${projectDir}/conda_env/extract_sc_env.yaml"
    publishDir "$params.publish_dir",  mode: 'copy'

    input:
    path unmapped_sam 

    output:
        path "unmapped_chunks/unmapped_chunk_*.fasta"
        path "reads_missing_cb.txt"

    """
    mkdir unmapped_chunks
    sc_split_reads.py $unmapped_sam unmapped_chunks $params.n_reads_per_chunk
    """
}

process sc_map_unmapped_reads {
    cpus 2
    memory '10 GB'
    time '24 hours'
    publishDir "${params.publish_dir}/mapped_chunks",  mode: 'copy'

    input:
        path unmapped_fasta

    output:
        path "${unmapped_fasta.baseName}_reads_barcodes.txt"

    """
    #!/usr/bin/bash
    flexiplex \
            -p ${params.adapter_5prime} \
            -T ${params.adapter_3prime} \
            -b 20 \
            -u 0 \
            -f ${params.adapter_edit_distance} \
            -e 3 \
            -n ${unmapped_fasta.baseName} \
            -k ${params.clone_barcodes_reference} \
            $unmapped_fasta
    
    """
}

process sc_merge_barcodes {
    cpus 4
    memory '24 GB'
    time '1 hours'
    conda "${projectDir}/conda_env/extract_sc_env.yaml"
    publishDir "${params.publish_dir}",  mode: 'copy'

    input:
        path mapped_reads 

    output:
        path "clone_barcodes.csv"

    """
    sc_merge_clone_barcodes.py $mapped_reads
    """
}