#!/usr/bin/env nextflow

process ngs_trim_reads {
    cpus 8
    memory '8 GB'
    time '2 hours'
    conda "${projectDir}/conda_env/trimgalore_env.yaml"
    publishDir "${params.publish_dir}", mode: 'copy'
    module "cutadapt"

    input:
        path fastq_file

    output:
        path "${out_file}"

    script:
    out_file = fastq_file.getSimpleName() + ".${params.trim_reads_to}bp_5prime.fq.gz"
    
    """
    trim_galore --hardtrim5 ${params.trim_reads_to} \
                --cores ${task.cpus} $fastq_file \
                --gzip
    
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
        path "${out_file}"

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
    // Add dummy adapters, run flexiplex discovery
    cpus 4
    memory '8 GB'
    time '1 hours'
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        path fastq_file

    output:
        path "${sample_name}_barcodes_counts.txt"

    script:
    sample_name = fastq_file.getSimpleName()
    fastq_w_adapter = sample_name + "_wDummyAdaptor.fq"
    """
    zcat $fastq_file | sed 's/^/START/g' | sed 's/START@/@/g' > ${fastq_w_adapter}
    
    flexiplex \
        -p "START" \
        -T "" \
        -b 20 \
        -u 0 \
        -f 0 \
        -n $sample_name \
        ${fastq_w_adapter}
    
    """
}

process ngs_split_reads_to_chunks {
    // Add dummy adapters, run flexiplex discovery, break up the barcodes into chunks
    cpus 4
    memory '8 GB'
    time '1 hours'
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        path barcode_counts

    output:
        path "${outdir}/${barcode_counts.baseName}_unmapped_chunk_*.fasta"

    script:
        outdir = "${barcode_counts.baseName}_unmapped_chunks"

    """
    mkdir ${outdir}
    ngs_break_up_barcodes.py --barcode_file ${barcode_counts} \
                                --sample_name ${barcode_counts.baseName} \
                                --nreads_per_chunk ${params.nreads_per_chunk} \
                                --outdir ${outdir}
    
    """
}

process ngs_map_barcodes {
    // Ran flexiplex per fasta chunk
    // Then combine the counting of read (flexiplex discovery)
    // and the mapped barcode
    cpus 2
    memory '4 GB'
    time '24 hours'
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
        path unmapped_fasta

    output:
        path "${out_file}"

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
        ${unmapped_fasta}

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