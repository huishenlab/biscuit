---
title: Read Mapping
nav_order: 2
has_children: true
permalink: docs/alignment
---

# Read Mapping

There are two subcommands used when mapping reads to a reference genome, `index`
and `align`. `index` is used to index the chosen reference genome to allow
for quicker mapping during the alignment process. `align` performs the actual
mapping of reads to the indexed reference genome.

## Step 1: Index Reference for Mapping

Before performing the read alignment, you must create create an index of the
reference genome (called `my_reference.fa` below) to be used in your alignment:

```bash
$ biscuit index my_reference.fa
```
The index of BISCUIT is composed of the 2-bit packed reference (`*.bis.pac`,
`*.bis.amb`, `*.bis.ann`). The FM-index and suffix array of the parent strand
(`*.par.bwt` and `*.par.sa`) and the daughter strand (`*.dau.bwt` and `*.dau.sa`).
Here, `*` refers to the input reference file, `my_reference.fa`.

For more help with `index`, run `biscuit index` in the terminal.

## Step 2: Read Mapping

After creating an index of the reference genome, the sequenced reads can be
mapped. In conjunction with [samtools](https://github.com/samtools/samtools)
and [samblaster](https://github.com/GregoryFaust/samblaster), BISCUIT can be
used to map, mark duplicate, and sort the provided reads.  Additionally,
`samtools` can be used to create an indexed BAM file.

```bash
$ biscuit align -R "my_rg" /path/to/my_reference.fa read1.fq.gz read2.fq.gz |
    samblaster | samtools sort -o my_output.bam -O BAM -
$ samtools index my_output.bam
```
The first line is referred to as the *biscuitBlaster* pipeline and provides
the user with a sorted BAM file of aligned reads with duplicates marked, but not
removed (see
[The biscuitBlaster Pipeline[({{ site.baseurl }}{% link docs/alignment/biscuitblaster.md %})
for more information). During the
[Read Pileup]({{ site.baseurl }}{% link docs/pileup.md %})
state, BISCUIT will skip marked duplicates, by default. If desired, `samblaster`
has a flag for removing duplicates during the marking process. If duplicates
were retained and one desires to retain the marked duplicates during the pileup
stage, a flag has been included for `biscuit pileup` to retain marked duplicates.

For more help with `align`, run `biscuit align` in the terminal.

### Mapping a Single Read

While the above example highlights using FASTQ files for read mapping, BISCUIT
does not require a FASTQ file to map reads. Instead, you can use the `-1` option
in `biscuit align` to map single reads on the fly.

```bash
$ biscuit align /path/to/my_reference.fa -1 AATTGGCC
```
You can also map one *pair* of reads with an extra `-2` option.

## Which Strand to Map?

Depending on the library construction strategy, bisulfite-converted reads can
come from one of four types of strands. These types are determined by whether the
target is the Waston or Crick strand and whether the DNA is the original strand
that underwent bisulfite conversion (__parent strand__), or the synthesized
strand during PCR amplification (__daughter strand__). In order to properly use a
read mapper, it is critical to understand the strands from which the sequenced
reads can be generated. In BISCUIT, whether mapping occurs at the parent or
daughter strand is controlled using `-b` option.

### Single-End Library

By default, BISCUIT maps reads to both strands (`-b 0`). This works for most
cases, including coventional libraries (such as TCGA WGBS and Roadmap Epigenome
Project), the PBAT library, and single cell libraries (where tagging happens
after amplification). When `-b 1` is specified, reads are forced to map only to
the parent strands. This is more efficient and less error-prone for conventional
libraries, such as TCGA and the Roadmap Project, as their library preparations
only generate reads from the parent strands. If, for some reason, you need to
force reads to be mapped to only the daughter strands, then use `-b 3` in your
options for `biscuit align`.

### Paired-End Library

When `-b 0` is used (the default behavior), BISCUIT maps read 1 (from
`read1.fq.gz`) to one strand (parent or daughter) and read 2 (from `read2.fq.gz`)
to the other strand. Alternatively, when using `-b 1`, read 1 is forced to map to
the parent strand and read 2 is mapped to the daughter strand. In order to map
read 1 to the daughter and read 2 to the parent, simply swap `read1.fq.gz` and
`read2.fq.gz` in the command line prompt.

### An Example

An example of using the default funcionality in BISCUIT is:

```bash
$ biscuit align mm10.fa -1 TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT
```
```
# Note, this is a reduced representation of the output
inputread    0    chr10    3386516    3    35M    *    0    0
    TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT    *    NM:i:0    MD:Z:35    ZC:i:1
    ZR:i:0    AS:i:35    XS:i:34    XL:i:35    XA:Z:chr10,+3386516,35M,0
    XB:Z:1,0     YD:A:f
```

In this example, you can see that, without specifying the strand, BISCUIT maps
the read to the bisulfite Watson (`YD:A:f`) strand automatically, without any
mismatches (`NM:i:0`).

Conversely, when using `-b 3`:

```bash
$ biscuit align -b 3 mm10.fa -1 TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT
```
```
inputread    0    chr10    3386516    60    35M    *    0    0
    TTGGTGTGTGGGTTTTGATGTTGGGTGGAGGGTTT    *    NM:i:0    MD:Z:35    ZC:i:0
    ZR:i:17    AS:i:34    XS:i:0    XL:i:35    YD:A:r
```

In this case, BISCUIT maps the read to the bisulfite Crick (`YD:A:r`) strand
with one mismatch (`NM:i:1`).

## Distinguishing Decoy Chromosomes in Human and Mouse Genomes

As is done in BWA, BISCUIT makes a distinction of primary chromosomes from
alternative/decoy chromosomes. For human and mouse genomes, this inference is
turned on by default, which results in a preferential alignment of reads to the
primary chromosomes. This can be turned off through the `-i` option in
`biscuit align`. The inference of alternative/decoy chromosomes is (generically)
done in the following manner: If the chromosome names are listed as `chr1`,
`chr2`, and so on, then chromosomes with the name pattern `chrUn`, `_random`,
`_hap`, or `_alt` are set as ALT chromosomes.

## Other Useful Options

For more options available for indexing and aligning with BISCUIT, run
`biscuit index` or `biscuit align` in your terminal.

## Other Alignment Features

  - Asymmetric scoring for C to T and G to A is used in local mapping
  - Generates a consistent mapping quality calculation with destination strand
  specification
  - Produces consistent `NM` and `MD` tags under asymmetric scoring
  - Produces `ZR` and `ZC` tags for retention count and conversion count
  - Separate seeding for parent and daughter strands for mapping efficiency
  - Indices make efficient use of disk-space and there is no need to store a
  bisulfite-converted reference
  - BWA-mem-like parameters that are visible to the users
  - Rare dependencies
  - Robust to OS build and stable multi-threading
  - Optional parent and daughter strand restriction for both single- and
  paired-end reads
  - Optional BSW/top/BSC/bottom strand restriction, tightly integrated in mapping
