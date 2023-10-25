#!/usr/bin/env nextflow

process ngs_trim_reads {
    label 'medium'
    conda "${projectDir}/conda_env/trimgalore_env.yaml"
    module "cutadapt"

    input:
        path fastq_file

    output:
        path "${out_file}"

    script:
    out_file = fastq_file.getSimpleName() + ".${params.barcode_length}bp_5prime.fq.gz"
    
    """
    trim_galore --hardtrim5 ${params.barcode_length} \
                --cores ${task.cpus} $fastq_file \
                --gzip
    
    """
}

process ngs_filter_reads {
    label 'medium'
    conda "${projectDir}/conda_env/fastp_env.yaml"

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
    label 'small'

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
        -b ${params.barcode_length} \
        -u 0 \
        -f 0 \
        -n $sample_name \
        ${fastq_w_adapter}
    
    """
}

process ngs_split_reads_to_chunks {
    // break up the barcodes into chunks
    label 'small'
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"

    input:
        path barcode_counts

    output:
        path "${outdir}/${barcode_counts.baseName}_chunk*.fasta"

    script:
        outdir = "${barcode_counts.baseName}_unmapped_chunks"

    """
    mkdir ${outdir}
    ngs_split_reads.py --barcode_file ${barcode_counts} \
                                --sample_name ${barcode_counts.baseName} \
                                --n_chunks ${params.n_chunks} \
                                --outdir ${outdir}
    
    """
}

process ngs_map_barcodes {
    // Ran flexiplex per fasta chunk
    // Then combine the counting of read (flexiplex discovery)
    // and the mapped barcode
    label "${params.mapping_process_profile}"
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"

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
        -l "START" \
        -r "" \
        -b ${params.barcode_length} \
        -u 0 \
        -f 0 \
        -n ${sample_name} \
        -k ${params.clone_barcodes_reference} \
        -e ${params.barcode_edit_distance} \
        -p ${task.cpus} \
        ${unmapped_fasta}

    ngs_combine_read_cnt_map.py --unmapped_chunk ${unmapped_fasta} \
                                --mapped_chunk ${mapped_chunk} \
                                --out_file ${out_file}
    """
}

process ngs_collapse_barcodes {
    label 'small'
    conda "${projectDir}/conda_env/extract_ngs_env.yaml"

    input:
        path mapped_reads

    output:
        path "clone_barcode_counts.csv"

    script:
    """
    ngs_count_barcodes.py . ${mapped_reads}
    """
}