#!/usr/bin/env nextflow

process sc_get_unmapped_reads {
    // Using sambamba
    module 'sambamba'
    cpus 12
    memory '32 GB'
    time '24 hours'
    publishDir "$params.publish_dir",  mode: 'copy'

    input:
        path bam_file

    output:
        path "${out_bam_file}"

    script:
        out_bam_file = "${bam_file.baseName}_reads_unmapped.bam"
    """
    sambamba view -t ${task.cpus} -F "unmapped" -f bam -o ${out_bam_file} ${bam_file}
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
        path "${out_bam_file}"

    script:
        out_bam_file = "${unmapped_bam.baseName}_filtered.bam"

    """
    sc_remove_low_qual_reads.py ${unmapped_bam} ${params.phred_thres} ${out_bam_file}
    """
}

process sc_split_unmapped_reads {
    cpus 4
    memory '24 GB'
    time '3 hours'
    conda "${projectDir}/conda_env/extract_sc_env.yaml"
    publishDir "$params.publish_dir",  mode: 'copy'

    input:
        path unmapped_bam 

    output:
        path "${outdir}/${unmapped_bam.baseName}_unmapped_chunk_*.fasta"
        path "${reads_missing_cb_filename}"

    script:
        outdir = "${unmapped_bam.baseName}_unmapped_chunks"
        reads_missing_cb_filename = "${unmapped_bam.baseName}_reads_missing_cb.txt"
    """
    mkdir ${outdir}
    touch ${reads_missing_cb_filename}
    sc_split_reads.py ${unmapped_bam} ${outdir} ${params.nreads_per_chunk} ${reads_missing_cb_filename}
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
            -p ${params.adapter_5prime_clonmapper} \
            -T ${params.adapter_3prime_clonmapper} \
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
        path "${outfile}"

    script:
        outfile = "clone_barcodes.csv"

    """
    sc_merge_clone_barcodes.py ${mapped_reads} ${outfile}
    """
}