---
title: BAM Operation
parent: Read Mapping
nav_order: 4
permalink: /bamoperation/
---

# BAM Operation

### Mark duplicate reads

#### biscuitBlaster Pipeline

The biscuitBlaster pipeline combines BISCUIT, samblaster, and samtools to mark duplicates and sort aligned reads.
The pipeline is:
```bash
$ biscuit align -M -R "my_rg" /path/to/idx read1.fq.gz read2.fq.gz | 
    samblaster -M | samtools sort -o my_output.bam -O BAM -
```
where `"my_rg"` is the read group (if applicable) to be used,
`/path/to/idx` is the FASTA file for your reference genome,
`read*.fq.gz` are the read1 and read2 FASTQ files from your sequencing run, and
`my_output.bam` is the name of your output BAM file.

Representation of the serial biscuitBlaster pipeline:

![serial biscuitBlaster](/biscuit/assets/serial_cookie_monster.gif)

Representation of the parallel biscuitBlaster pipeline:

![parallel biscuitBlaster](/biscuit/assets/parallel_cookie_monster.gif)

#### Deprecated functionality

**This functionality is DEPRECATED and will be removed in the next release!!!!
We suggest using the biscuitBlaster pipeline for duplicate marking.**
*This step is optional. The mark duplicate of BISCUIT is bisulfite strand aware.*
```bash
$ biscuit markdup input.bam output.bam
```

### Filter reads by bisulfite conversion

```bash
$ biscuit bsconv -g "chr1:1000000-1050000" GRCh37.fa input.bam
```
For some library preparation, incomplete conversion are enriched in a subset of reads that needs to be filtered. This command transforms bam into one that contains file `ZN` tag e.g., `ZN:Z:CA_R0C11,CC_R1C14,CG_R0C2,CT_R1C5`. This tag summarizes counts of retention and conversion for four different cytosine contexts `CpA`, `CpC`, `CpG` and `CpT`. It contains a minimum threshold of `CpA`, `CpC`, `CpT` or `CpH` in general. The `-b` option outputs the summary in tables instead of as tags in the BAM file.
