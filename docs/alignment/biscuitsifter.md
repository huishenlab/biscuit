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
biscuit align \
    -@ NTHREADS \
    /path/to/my_reference.fa \
    read1.fq.gz \
    read2.fq.gz | \
dupsifter /path/to/my_reference.fa | \
samtools sort -@ NTHREADS -o my_output.bam -O BAM -

samtools index my_output.bam
```
where `NTHREADS` is the number of threads, `/path/to/my_reference.fa` is the FASTA file for the reference genome,
`read*.fq.gz` are the FASTQ files for reads 1 and 2 from the sequencing run, and `my_output.bam` is the name of the
output BAM file.

The `--write-index` option in `samtools sort` can alternatively be used to simplify the two-step process into a single
command:
```bash
biscuit align \
    -@ NTHREADS \
    /path/to/my_reference.fa \
    read1.fq.gz \
    read2.fq.gz | \
dupsifter /path/to/my_reference.fa | \
samtools sort \
    -@ NTHREADS \
    --write-index \
    -o my_output.bam \
    -O BAM \
    -
```
