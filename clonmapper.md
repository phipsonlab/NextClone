---
layout: default
---

## Vignette for Clonmapper library

If you used the [Clonmapper](https://doi.org/10.1038/s43018-021-00222-8) library to barcode your cells, follow this vignette to use NextClone.

First thing first:

1. Copy the `nextflow.config` file on our github repository: https://github.com/phipsonlab/NextClone.
2. Set the directory/file path in the config file for the following parameters: 
   1. `publish_dir`
   2. `clone_barcodes_reference`
   3. `dnaseq_fastq_files` for DNA-seq data or `scrnaseq_bam_files` for scRNA-seq data.

Read the documentations on the homepage, make sure you understand what each of the parameter in the `nextclone.config` file is for.

Finally read on.

### For DNA-seq data

Make sure **at least** the following parameters are set in `nextflow.config` file.

```
mode = "DNAseq"
barcode_length = 20
```

If your clone barcode is not 20 bp, change `barcode_length` above to 20.

Then pull the github code in the main branch and run.

```
nextflow run phipsonlab/Nextclone -r main
```

### For scRNA-seq data

Before running NextClone, make sure you run cellranger first, and copy out the `possorted_genome_bam` file in the `outs` folder.

Then make sure **at least** the following parameters are set in `nextflow.config` file.

```
mode = "scRNAseq"
barcode_length = 20
adapter_5prime_clonmapper = "ATCTTGTGGAAAGGACGAAACACCG"
adapter_3prime_clonmapper = "GTTTCAGAGCTATGCTGGAAACAGC"
```

If your clone barcode is not 20 bp, change `barcode_length` above to 20.

You almost definitely will also need to change the 5' and 3' adapter sequence for your clone barcodes.
This is denoted by `adapter_5prime_clonmapper` and `adapter_3prime_clonmapper` parameter.
This sequence is determined by whoever made your library, and you should ask them for the sequence.
The default sequence in the `nextflow.config` file is unlikely to match the one that is in your ClonMapper library.

Then pull the github code in the main branch and run.

```
nextflow run phipsonlab/Nextclone -r main
```

[back](./)