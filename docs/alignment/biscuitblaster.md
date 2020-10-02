---
title: The biscuitBlaster Pipeline
parent: Read Mapping
nav_order: 1
permalink: /biscuitblaster/
---

# The biscuitBlaster Pipeline

The biscuitBlaster pipeline combines BISCUIT, samblaster, and samtools to align
reads, mark duplicates, and sort and index the aligned reads in an easy
one-liner.

Depending on your downstream analysis goals, there are a few pipelines available.

## Version 1

Version 1 of the biscuitBlaster pipeline is the most basic version, producing an
sorted, indexed BAM with duplicate marked reads.

```bash
$ biscuit align -R "my_rg" /path/to/my_reference.fa read1.fq.gz read2.fq.gz | 
    samblaster | samtools sort -o my_output.bam -O BAM -
$ samtools index my_output.bam
```
where `"my_rg"` is the read group (if applicable) to be used,
`/path/to/my_reference.fa` is the FASTA file for your reference genome,
`read*.fq.gz` are the read1 and read2 FASTQ files from your sequencing run, and
`my_output.bam` is the name of your output BAM file.

## Version 2

In addition to generating a duplicate marked, sorted, and indexed BAM file,
with mate tags added for use with structural variant finders, Version 2 builds
off Version 1 by using samblaster to generate SAM files of the split and
discordantly aligned reads and a FASTQ file of the heavily clipped reads.

```bash
$ biscuit align -R "my_rg" \
    /path/to/my_reference.fa read1.fq.gz read2.fq.gz | \
    samblaster --addMateTags | \
    parallel --tmpdir temp_dir --pipe --tee {} ::: \
        "samblaster -a -e -u clipped.fastq -d disc.sam -s split.sam -o /dev/null" \
        "samtools view -hb | samtools sort -o all_my_reads.bam -O BAM -"
$ samtools index all_my_reads.bam
```

Version 2 leverages GNU parallel to produce each of these files in a fast,
efficient manner. Note, each of the SAM files produced will likely need to be
sorted before using in downstream analyses.

## biscuitBlaster Representations

### Version 1 Representation

![serial biscuitBlaster](/biscuit/assets/serial_cookie_monster.gif)

### Version 2 Representation

![parallel biscuitBlaster](/biscuit/assets/parallel_cookie_monster.gif)
