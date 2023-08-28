#!/usr/bin/env nextflow

process ngs_trim_reads {
    cpus 8
    memory '8 GB'
    time '2 hours'
    conda "${projectDir}/conda_env/trimgalore_env.yaml"
    publishDir "${params.publish_dir}"
    module "cutadapt"

    input:
        path fastq_file

    output:
        path "${out_file}"

    script:
    out_file = fastq_file.getSimpleName() + "_trimmed.fq.gz"
    orig_out_file = fastq_file.getSimpleName() + ".${params.trim_reads_to}bp_5prime.fq.gz"
    
    """
    trim_galore --hardtrim5 ${params.trim_reads_to} \
                --cores ${task.cpus} $fastq_file \
                --gzip
    
    mv ${orig_out_file} ${out_file}
    
    """
}

process ngs_filter_reads {
    cpus 8
    memory '8 GB'
    time '2 hours'
    conda "${projectDir}/conda_env/fastp_env.yaml"
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        path fastq_file

    output:
        path out_file

    script:
    out_file = fastq_file.getSimpleName() + "_filtered.fq.gz"
    
    """
    fastp -i $fastq_file \
            -o $out_file \
            -q ${params.fastp_phred_for_qualified_reads} \
            -u ${params.fastp_percent_bases_unqualified} \
            -w ${task.cpus}
    """
}

process ngs_count_reads {
    // Add dummy adapters, run flexiplex discovery, break up the barcodes into chunks
    cpus 1
    memory '2 GB'
    time '1 hours'
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"
    tag "${sample_name}"
    publishDir "${params.publish_dir}/${sample_name}/fasta_chunks", mode: 'copy'

    input:
        path fastq_file

    output:
        path '*_chunk*.fasta'

    script:
    sample_name = fastq_file.getSimpleName()
    fastq_w_adapter = sample_name + "_wDummyAdaptor.fq"
    barcode_counts = sample_name + "_barcodes_counts.txt"
    """
    zcat $fastq_file | sed 's/^/START/g' | sed 's/START@/@/g' > $fastq_w_adapter
    
    flexiplex \
        -p "START" \
        -T "" \
        -b 20 \
        -u 0 \
        -f 0 \
        -n $sample_name \
        ${fastq_w_adapter}

    ngs_break_up_barcodes.py --barcode_file ${barcode_counts} \
                                --sample_name ${sample_name} \
                                --nreads_per_chunk ${params.nreads_per_chunk}
    
    """
}

process ngs_map_barcodes {
    // Ran flexiplex per fasta chunk
    // Then combine the counting of read (flexiplex discovery)
    // and the mapped barcode
    cpus 20
    memory '4 GB'
    time '24 hours'
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"
    tag "${sample_name}"
    publishDir "${params.publish_dir}/${sample_name}/mapped_barcodes", mode: 'copy'

    input:
        path unmapped_fasta

    output:
        path out_file

    script:
    sample_name = unmapped_fasta.getSimpleName()
    mapped_chunk = sample_name + "_reads_barcodes.txt"
    out_file = sample_name + "_mapped.csv"
    """

    flexiplex \
        -p "START" \
        -T "" \
        -b 20 \
        -u 0 \
        -f 0 \
        -n ${sample_name} \
        -k ${params.clone_barcodes_reference} \
        -e ${params.barcode_edit_distance} \
        -p ${task.cpus} \
        $unmapped_fasta

    ngs_combine_mapped_and_discovery_barcodes.py --unmapped_chunk ${unmapped_fasta} \
                                                    --mapped_chunk ${mapped_chunk} \
                                                    --out_file ${out_file}
    """
}

process ngs_collapse_barcodes {
    cpus 1
    memory '4 GB'
    time '1 hours'
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        path mapped_reads

    output:
        path "clone_barcode_counts.csv"

    script:
    """
    ngs_count_barcodes.py ${mapped_reads}
    """
}