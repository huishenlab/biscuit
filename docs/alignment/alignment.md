---
title: Read Mapping
nav_order: 3
has_children: true
permalink: docs/alignment
---

# Read Mapping

There are two subcommands used when mapping reads to a reference genome, `index` and `align`. `index` is used to index
the specified reference genome to allow for quicker mapping during the alignment process. `align` performs the actual
alignment of reads to the indexed reference genome.

## Creating the Reference Index

Before performing the read alignment, an index of the reference genome (referred to as `my_reference.fa` below) must be
created. This step only needs to be run **once** for each reference.

To create the index, run
```bash
biscuit index my_reference.fa
```

The index of BISCUIT is composed of the 2-bit packed reference (`*.bis.pac`, `*.bis.amb`, `*.bis.ann`). The FM-index and
suffix array of the parent strand (`*.par.bwt` and `*.par.sa`) and the daughter strand (`*.dau.bwt` and `*.dau.sa`).
Here, `*` refers to the input reference file, `my_reference.fa`.

For more help with `index`, run `biscuit index` in the terminal or check out the
[index help page]({{ site.baseurl }}{% link docs/subcommands/biscuit_index.md %}).

## Aligning Reads to the Reference

After creating an index of the reference genome, the sequenced reads can be mapped to the reference genome. In
conjunction with [samtools](https://github.com/samtools/samtools) and
[dupsifter](https://github.com/huishenlab/dupsifter), BISCUIT can be used to map, duplicate mark, sort, and index the
provided reads.

The suggested one-line command to generate an aligned, duplicate marked, and sorted BAM is referred to as the
[biscuitSifter]({{ site.baseurl }}{% link docs/alignment/biscuitsifter.md %}) pipeline:
```bash
# Uncompressed single end FASTQ example
biscuit align \
    -@ NTHREADS \
    /path/to/my_reference.fa \
    read1.fq | \
dupsifter /path/to/my_reference.fa | \
samtools sort -@ NTHREADS -o my_output.bam -O BAM -

# Uncompressed paired end FASTQ example
biscuit align \
    -@ NTHREADS \
    /path/to/my_reference.fa \
    read1.fq \
    read2.fq | \
dupsifter /path/to/my_reference.fa | \
samtools sort -@ NTHREADS -o my_output.bam -O BAM -

# Gzipped paired end FASTQ example
biscuit align \
    -@ NTHREADS \
    /path/to/my_reference.fa \
    read1.fq.gz \
    read2.fq.gz | \
dupsifter /path/to/my_reference.fa | \
samtools sort -@ NTHREADS -o my_output.bam -O BAM -
```
Indexing the aligned BAM is done with samtools:
```
samtools index my_output.bam
```

In its suggested form, the biscuitSifter pipeline creates a BAM file that retains the duplicate marked reads. `biscuit
pileup` will ignore marked duplicates by default, though these reads can be included in the analysis with the `-u` flag.
More details on `biscuit pileup` can be found on the [Read Pileup]({{ site.baseurl }}{% link docs/pileup.md %}) page.

For more help with `align`, run `biscuit align` in the terminal or check out the
[align help page]({{ site.baseurl }}{% link docs/subcommands/biscuit_align.md %}).

### Mapping a Single Read

While the above example highlights using FASTQ files for read mapping, BISCUIT can also map single reads on the fly
using the `-1` option for a single end read or `-1` and `-2` for a pair of reads from paired end sequencing.

```bash
# Single end example
biscuit align \
    /path/to/my_reference.fa \
    -1 AATTGGCC

# Paired end example
biscuit align \
    /path/to/my_reference.fa \
    -1 AATTGGCC \
    -2 AGCTAGCT
```

## Which Strand to Map?

Depending on the library construction strategy, bisulfite-converted (or bisulfite-like) reads can come from one of four
types of strands. These types are determined by whether the target is the top or bottom strand and whether the DNA is
the original strand that underwent conversion (**parent strand**, original top (OT) or original bottom (OB)), or the
complement strand created during PCR amplification (**daughter strand**, complement to original top (CTOT) or complement
to original bottom (CTOB)). In order to properly use a read aligner, it is critical to understand the strands from which
the sequenced reads can be generated. In BISCUIT, the strands the reads can be mapped to is controlled by the `-b`
option.

### Single-End Library

By default, BISCUIT maps reads to all four strands (OT/CTOT/OB/CTOB, `-b 0`).  This works for most cases, including
coventional libraries (such as TCGA WGBS and Roadmap Epigenome Project), the PBAT library, and single cell libraries
(where tagging happens after amplification). When `-b 1` is specified, reads are forced to map only to the parent
strands (OT/OB). This is more efficient and less error-prone for conventional libraries, such as TCGA and the Roadmap
Project, as their library preparations only generate reads from the parent strands. If, for some reason, reads must be
mapped to only the daughter strands (CTOT/CTOB), then include `-b 3` in the options for `biscuit align`.

### Paired-End Library

When `-b 0` is used (the default behavior), BISCUIT maps read 1 (from `read1.fq.gz`) to one strand (parent or daughter)
and read 2 (from `read2.fq.gz`) to the other strand. Alternatively, when using `-b 1`, read 1 is forced to map to
the parent strand and read 2 is mapped to the daughter strand. In order to map read 1 to the daughter and read 2 to the
parent, simply swap `read1.fq.gz` and `read2.fq.gz` in the command line prompt.

### An Example

An example of using the default functionality in BISCUIT is:
```bash
# Align one read to the mouse genome with default parameters
biscuit align \
    mm10.fa \
    -1 TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT

# The output below is a reduced representation of the full output
inputread    0    chr10    3386516    3    35M    *    0    0
    TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT    *    NM:i:0    MD:Z:35    ZC:i:1
    ZR:i:0    AS:i:35    XS:i:34    XL:i:35    XA:Z:chr10,+3386516,35M,0
    XB:Z:1,0    MC:Z:*    MQ:i:0    YD:A:f
```
In this example where the strand is not specified, BISCUIT maps the read to the original top (`YD:A:f`) strand without
any mismatches (`NM:i:0`).

Conversely, when using `-b 3`:
```bash
# Align one read to the mouse genome with alignment to the daughter strand
biscuit align \
    -b 3 \
    mm10.fa \
    -1 TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT

# Reduced representation output
inputread    0    chr10    3386516    60    35M    *    0    0
    TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT    *    NM:i:0    MD:Z:35    ZC:i:0
    ZR:i:17    AS:i:34    XS:i:0    XL:i:35    MC:Z:*    MQ:i:0    YD:A:r
```
In this case, BISCUIT maps the read to the original bottom (`YD:A:r`) strand with one mismatch (`NM:i:1`).

## Distinguishing Alternative/Decoy Chromosomes

As is done in BWA, BISCUIT makes a distinction between primary/unplaced/unlocalized contigs and alternative/decoy
contigs using the `my_reference.fa.alt` file. Any contig names found in this file will be marked as ALT contigs. BISCUIT
will preferentially align reads to contigs not marked as ALT. If no `my_reference.fa.alt` file is found, all contigs
will be treated as equally valid. [The FAQ page]({{ site.baseurl }}{% link docs/faq.md %}) contains information on how
to create a `my_reference.fa.alt` file if you do not have one.

As of version 1.10.0, BISCUIT provides the `biscuit postalt` subcommand to modify an output BAM from `biscuit align` to
prefer primary chromosomes over alternative/decoy chromosomes:
```bash
biscuit postalt \
    --unplaced \
    <in.bam> \
    [out.bam]
```
When a read's primary alignment is on a non-primary chromosome (as defined by an input regular expression) and its XA
tag contains a equal/near-equal score on a primary chromsome, BISCUIT will promote the alignment in the XA tag to the
primary alignment of the read and demote the original alignment to the XA tag. It will then recalculate the mapping
quality (MAPQ) score for the new primary alignment.

`biscuit postalt` provides an opt-in method for prioritizing primary chromosomes over unplaced/unlocalized contigs
without impacting the behavior of `biscuit align`.

## Other Alignment Features

  - Asymmetric scoring for C &rarr; T and G &rarr; A is used in local mapping
  - Generates a consistent mapping quality calculation with destination strand specification
  - Produces consistent `NM` and `MD` tags under asymmetric scoring
  - Produces `ZR` and `ZC` tags for retention count and conversion count
  - Separate seeding for parent and daughter strands for mapping efficiency
  - Indices make efficient use of disk-space and there is no need to store a bisulfite-converted reference
  - BWA-mem-like parameters that are visible to the users
  - Few dependencies (`zlib`, `ncurses`, `pthread`, and `curl`)
  - Robust to OS build and stable multi-threading
  - Optional parent and daughter strand restriction for both single- and paired-end reads
  - Optional OT/OB strand restriction, tightly integrated in mapping
  - Native inclusion of mate auxiliary tags (MC/MQ) in alignments
