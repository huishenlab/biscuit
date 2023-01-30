---
title: The biscuitBlaster Pipeline
parent: Read Mapping
nav_order: 1
permalink: /biscuitblaster/
---

# The biscuitBlaster Pipeline

The biscuitBlaster pipeline combines BISCUIT, samblaster, and samtools to align reads, mark duplicates, sort, and index
the aligned reads in an easy two-step process.

Depending on the downstream analysis goals of the user, there are two pipelines available.

## Version 1

Version 1 of the biscuitBlaster pipeline is the most basic version, producing a sorted, indexed BAM with duplicate
marked reads.

```bash
biscuit align -R "my_rg" /path/to/my_reference.fa read1.fq.gz read2.fq.gz | \
    samblaster | samtools sort -o my_output.bam -O BAM -
samtools index my_output.bam
```
where `"my_rg"` is the read group (if applicable) to be used, `/path/to/my_reference.fa` is the FASTA file for the
reference genome, `read*.fq.gz` are the read1 and read2 FASTQ files from the sequencing run, and `my_output.bam` is the
name of the output BAM file.

## Version 2

Version 2 of the biscuitBlaster pipeline builds off Version 1 by using samblaster generate SAM files of the split and
discordant alignments and a FASTQ file of the heavily clipped reads, in addition to generating a duplicate marked,
sorted, and indexed BAM file with mate tags added.

```bash
biscuit align -R "my_rg" \
    /path/to/my_reference.fa read1.fq.gz read2.fq.gz | \
    samblaster --addMateTags | \
    parallel --tmpdir temp_dir --pipe --tee {} ::: \
        "samblaster -a -e -u clipped.fastq -d disc.sam -s split.sam -o /dev/null" \
        "samtools view -hb | samtools sort -o all_my_reads.bam -O BAM -"
samtools index all_my_reads.bam
```
Version 2 leverages GNU parallel to produce each of these files in a fast, efficient manner. Note, each of the SAM files
produced will likely need to be sorted before using in downstream analyses.

It should also be pointed out that a relatively recent version of GNU parallel is needed to run Version 2 of
biscuitBlaster. In our experience, Version 2 has worked with GNU parallel versions >= 20191122.

## biscuitBlaster Representations

### Version 1 Representation

![serial biscuitBlaster](/biscuit/assets/serial_cookie_monster.gif)

### Version 2 Representation

![parallel biscuitBlaster](/biscuit/assets/parallel_cookie_monster.gif)
