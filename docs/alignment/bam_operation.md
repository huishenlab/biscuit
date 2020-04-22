---
title: BAM Operation
parent: Read Mapping
nav_order: 2
permalink: /bamoperation/
---

# BAM Operation

## Mark Duplicate Reads

### biscuitBlaster Pipeline

The biscuitBlaster pipeline combines BISCUIT, samblaster, and samtools to mark
duplicates and sort aligned reads. The pipeline is:
```bash
$ biscuit align -M -R "my_rg" /path/to/my_reference.fa read1.fq.gz read2.fq.gz | 
    samblaster -M | samtools sort -o my_output.bam -O BAM -
```
where `"my_rg"` is the read group (if applicable) to be used,
`/path/to/my_reference.fa` is the FASTA file for your reference genome,
`read*.fq.gz` are the read1 and read2 FASTQ files from your sequencing run, and
`my_output.bam` is the name of your output BAM file.

Representation of the serial biscuitBlaster pipeline:

![serial biscuitBlaster](/biscuit/assets/serial_cookie_monster.gif)

Representation of the parallel biscuitBlaster pipeline:

![parallel biscuitBlaster](/biscuit/assets/parallel_cookie_monster.gif)

### `biscuit markdup` - Deprecated Functionality

**This functionality is DEPRECATED and will be removed in the next release!!!!
We suggest using the biscuitBlaster pipeline for duplicate marking.**
```bash
$ biscuit markdup input.bam output.bam
```

## Filter Reads by Bisulfite Conversion

For some library preparations, incomplete conversions are enriched in a subset
of reads which need to be filtered. The `bsconv` subcommand transforms an input
BAM into one that contains a `ZN` tag (like
`ZN:Z:CA_R0C11,CC_R1C14,CG_R0C2,CT_R1C5`). This tag summarizes counts of
retention and conversion for four different cytosine contexts `CpA`, `CpC`,
`CpG` and `CpT`. The `-b` flag outputs the counts in a table, instead of as a
tag in the BAM file.
```bash
$ biscuit bsconv /path/to/my_reference.fa input.bam
```

For help on available flags, enter `biscuit bsconv` on the command line.
