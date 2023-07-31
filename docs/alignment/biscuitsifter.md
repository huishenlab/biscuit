---
title: The biscuitSifter Pipeline
parent: Read Mapping
nav_order: 1
permalink: /biscuitsifter/
---

# The biscuitSifter Pipeline

The biscuitSifter pipeline combines BISCUIT, dupsifter, and samtools to align reads, mark duplicates, sort, and index
the aligned reads in an easy two-step process.

```bash
biscuit align -@ NTHREADS -R "my_rg" /path/to/my_reference.fa read1.fq.gz read2.fq.gz | \
    dupsifter /path/to/my_reference.fa | samtools sort -@ NTHREADS -o my_output.bam -O BAM -
samtools index my_output.bam
```
where `NTHREADS` is the number of threads, `"my_rg"` is the read group (if applicable) to be used,
`/path/to/my_reference.fa` is the FASTA file for the reference genome, `read*.fq.gz` are the read1 and read2 FASTQ files
from the sequencing run, and `my_output.bam` is the name of the output BAM file.
