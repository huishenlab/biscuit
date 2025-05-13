---
title: Epireads and the epiBED Format
nav_order: 6
has_children: true
permalink: docs/epiread
---

# Epireads and the epiBED Format

Bisulfite sequencing data may be able to inform how methylation from neighboring positions are linked and how
methylation is linked to mutations (provided they can be unambiguously determined). We introduce the **epiBED** format,
which combines the utility of the [epiread format](http://smithlabresearch.org/downloads/methpipe-manual.pdf) with the
standard-compliant BED format. The epiBED format is a compact data format used to store CpG (and, in the case of
NOMe-seq, GpC) methylation, as well as SNP and indel information from the same read. It can be useful for estimating the
epiallele fraction and clonal structure/cell population, which is similar to
[MethylPurify](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0419-x).

## Generating Epireads in the epiBED Format

The original epiread format (proposed by the methpipe team), only includes information related to CpGs. BISCUIT builds
off this format to put it in a standard-compliant format, while also containing information related to filtered bases,
SNPs, indels, and, in NOMe-seq mode, GpC methylation. To generate an epiBED file, run
```bash
# For WGBS data
biscuit epiread -@ NTHREADS [-B snps.bed] /path/to/my_reference.fa my_output.bam | \
sort -k1,1 -k2,2n > my_epireads.epibed

# For NOMe-seq data
biscuit epiread -@ NTHREADS [-B snps.bed] -N /path/to/my_reference.fa my_output.bam | \
sort -k1,1 -k2,2n > my_epireads.epibed
```
The `snps.bed` file can be obtained by running `biscuit vcf2bed -t snp my_pileup.vcf.gz`. If no SNP file is given, the
output will not include information related to SNPs. Note, it is suggested that you sort your epiBED file before
running `bgzip` and `tabix` to compress and index your epiBED file. This is due to how softclipping is handled in the
BAM file versus how it is handled in the epiBED file. The difference can cause some reads to not be coordinate sorted,
which is required when bgzipping and tabixing.

As mentioned elsewhere, it is important to use the same reference to create the epiBED files as was used during the
alignment phase.

For more help with `epiread`, run `biscuit epiread` in the terminal or check out the
[epiread help page]({{ site.baseurl }}{% link docs/subcommands/biscuit_epiread.md %}).

## Generating epiBED Files from modBAMs

`biscuit epiread` can also be used to extract methylation from modification tags (MM/ML) in a BAM file via:
```bash
biscuit epiread -M -y 0.9 -L 10000 /path/to/my_reference.fa my_modbam.bam | \
sort -k1,1 -k2,2n > my_epireads.epibed
```
The `-M` flag turns on the modification tag reading and overrides reading methylation from C&rarr;T or G&rarr;A
conversions. The `-y 0.9` sets the probability of a successful modification call at 0.9 (which is the default value).
This probability can be adjusted anywhere from 0 to 1 (inclusive). Typically, modBAMs are associated with long read
sequencing, so the `-L` option increases the maximum read length to 10,000. This may be need to be adjusted depending
on the read length distribution of your dataset. If your data is short read sequencing, you can forego this option as
the default value is sufficient for short read lengths. While not required, some other options to consider modifying
when running `biscuit epiread` on long read data (particularly when comparing Nanopore data against
[modkit](https://nanoporetech.github.io/modkit/)) are `-b` (minimum base quality), `-m` (minimum mapping quality), `-a`
(minimum alignment score), `-5`/`-3` (exclude methylation within N bases of the 5' or 3' ends of the reads).

Currently, BISCUIT restricts the number of modifications to one per read and can be `C+m`, `C-m`, `G+m`, or `G-m`.
Coupled with the SAM FLAG, the position the methylation status is recorded at can either be at the C position or the G
position relative to the reference:

| Modification | SAM FLAG | Position |
|:------------:|:--------:|:--------:|
| `C+m`        | 0        | C        |
| `C+m`        | 16       | G        |
| `G-m`        | 0        | G        |
| `G-m`        | 16       | C        |

## Generating Legacy File Formats

In addition to the epiBED format, `biscuit epiread` continues to produce the BISCUIT epiread format and the pairwise
file format for `biscuit asm`. While these are available, it is suggested you use the epiBED format instead, as it
aligns with the output from `biscuit pileup`, whereas the epiread and pairwise file formats are not guaranteed to match
the results from `pileup`.
```bash
# Original epiread format with WGBS data
biscuit epiread -@ NTHREADS -o my_output.epiread -O [-B snps.bed] /path/to/my_reference.fa my_output.bam

# Original epiread format with NOMe-seq data
biscuit epiread -@ NTHREADS -o my_output.epiread -O -N [-B snps.bed] /path/to/my_reference.fa my_output.bam

# Pairwise format with WGBS data
biscuit epiread -@ NTHREADS -o my_output.epiread -P [-B snps.bed] /path/to/my_reference.fa my_output.bam

# Pairwise format with NOMe-seq data
biscuit epiread -@ NTHREADS -o my_output.epiread -P -N [-B snps.bed] /path/to/my_reference.fa my_output.bam
```
