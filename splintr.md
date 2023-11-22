---
layout: default
---

## Vignette for Splintr library

***WORK IN PROGRESS!***

If you used the [Splintr](https://doi.org/10.1038/s41586-021-04206-7) library to barcode your cells, follow this vignette to use NextClone.

First thing first:

1. Copy the `nextflow.config` file on our github repository: https://github.com/phipsonlab/NextClone.
2. Set the directory/file path in the config file for the following parameters: 
   1. `publish_dir`
   2. `clone_barcodes_reference`
   3. `dnaseq_fastq_files` for DNA-seq data or `scrnaseq_bam_files`

Then read on.

### For scRNAseq data

Before running NextClone, make sure you run cellranger first, and copy out the `possorted_genome_bam` file in the `outs` folder.

After you downloaded the FASTA file containing the reference clone barcodes, make sure you convert them to a normal text file where each line denote only the sequence for a clone barcode.
You can do this using a regex replacement. 
For example, using `vim` you can do something like this:
```
%s/>mCHERRY_Barcode_.*\n//n
```
What this will do is get rid of the name of the barcode and keep only the sequence. 

Then make sure **at least** the following parameters are set in `nextflow.config` file.

```
mode = "scRNAseq"
barcode_length = 60
adapter_5prime = "CGATTGACTA"
adapter_3prime = "TGCTAATGCG"
```

If your clone barcode is not 60 bp, change `barcode_length` above.

The sequence for the 5' and 3' adapter is from the original splintr protocol.
If you have asked for some customisation, you *probably* will have to modify the sequence above. 
This sequence is determined by whoever made your library.
Thus you should ask them for the sequence.

Then pull the github code in the `update_flexiplex_again` branch and run.

```
nextflow run phipsonlab/Nextclone -r update_flexiplex_again
```

[back](./)