---
title: Frequently Asked Questions
nav_order: 12
---

# Frequently Asked Questions

## General Questions

**I am new to DNA methylation research and do not understand a lot of the terminology. Are there sites that describe
relevant terms and phrases?**

> One place to look is the [Glossary]({{ site.baseurl }}{% link docs/glossary.md %}) page.

**Should I include the `-M` flag when running `biscuit align` or using the biscuitSifter pipeline?**

> Our suggestion is to avoid using the `-M` flag in BISCUIT unless working with legacy versions of tools or doing
structural variant calling. While current SV callers supposedly do not need the `-M` flag for compatibility (see
[this Twitter thread](https://twitter.com/biobenkj/status/1311660546232643584) for more details), our experience is that
the `-M` flag is still needed.

**I have a PBAT-style WGBS sample that I tried viewing in IGV's DNA methylation color scheme, but I see a lot of CpGs
colored brown, rather than the blue and red I expected. Why is that and is there anything I can do about it?**

> Generally, PBAT-style WGBS samples have reads 1 and 2 flipped relative to most other WGBS and WGBS-like protocols.
Because of this, the reads in BAMs from these samples have the forward and reverse strands flipped. IGV (at least up to
Version 2.8.2) can't handle these flipped strands and prints the color as brown, rather than blue or red.

> If you want to be able to view your samples in the blue/red color scheme, a script, `flip_pbat_strands.sh`, is
provided in the `scripts/` directory and on the release page. This script removes the "read reverse strand" flag if it
was included in the BAM FLAG field (or adds it if it was not included). Run
```bash
$ bash flip_pbat_strands.sh -h
```
for more details on usage.

## Installation Questions

**I ran `git clone git@github.com:huishenlab/biscuit.git`, but when trying to compile the executable, I get an error
message.**

> After version 1.4.0, make sure you `zlib`, `ncurses`, and `pthread` installed. If prior to version 1.4.0, make sure
you included `--recursive` when running `git clone`.

## Quality Control Related Questions

**I downloaded one of the pre-compiled binaries for BISCUIT, but want to run QC on my data. Where can I find the QC
script?**

> You can downloaded the `QC.sh` script from the
[BISCUIT release page](https://github.com/huishenlab/biscuit/releases/latest). Supplementary QC asset files will be
needed to run the script and can be found for the hg19, hg38, and mm10 genomes alongside the QC script. Note, if you
plan to always include `--no-cov-qc` when running `QC.sh`, you do not need to download or create the asset files.

**I have aligned my data to a reference genome that does not have QC asset files provided. What files do I need to run
`QC.sh`?**

> As of version 1.0.0, you will need three BED files. You will need a gzipped BED file with the locations of all
CpGs in your genome (called `cpg.bed.gz`). You will also need to break your genome into 100bp non-overlapping windows,
calculate the GC-content for each window, then find the windows with the top and bottom 10% of GC-content. These regions
will be placed into two separate gzipped BED files, `windows100bp.gc_content.bot10p.bed.gz` for the bottom 10% and
`windows100bp.gc_content.top10p.bed.gz` for the top 10%.  Make sure to sort each BED file using
(`sort -k1,1 -k2,2n bedfile.bed`).

> There is now a Perl script that can generate the asset files needed to run QC for BISCUIT. The
`build_biscuit_QC_assets.pl` script requires the name of the output directory to save the assets files to and the path
to the reference genome.

**I need to generate my own QC asset files and really want to make them myself without running the provided Perl script.
What is a more prescriptive way for creating these files?**

> **cpg.bed.gz**
  1. Find positions of all CpGs in genome reference FASTA.
  2. Put positions in
  [BED format](https://en.wikipedia.org/wiki/BED_(file_format)). You only need
  the chromosome, start, and end positions. Note, the BED format uses 0-based
  indexing, so make sure to number the first base as 0.
  3. Sort the BED file by chromosome position (`sort -k1,1 -k2,2n cpg.bed`).
  4. Gzip your sorted BED file and you now have your `cpg.bed.gz` QC file.

> **windows100bp.gc_content.top10p.bed.gz**
  1. Group genome reference FASTA into 100 bp windows.
  2. Calculate GC-content fraction for each window.
  3. Create BED file with window positions (chromosome, start, and end) and GC-content fraction as a fourth column.
  4. Sort BED file by GC-content fraction (`sort -k4,4n gc_content.bed`).
  5. Find top 10% of GC-content windows.
  6. Copy the four columns (chromosome, start, end, and GC-content fraction) for these windows into
  `windows100bp.gc_content.top10p.bed`.
  7. Sort by position (`sort -k1,1 -k2,2n`) and gzip `windows100bp.gc_content.top10p.bed` to create your top 10%
  GC-content QC file.

> **windows100bp.gc_content.bot10p.bed.gz**
  1. Follow Steps 1-4 for creating `windows100bp.gc_content.top10p.bed.gz` QC file.
  2. Instead of finding the top 10% of GC-content windows, find the bottom 10% of GC-content windows.
  3. Copy the four columns (chromosome, start, end, and GC-content fraction) for these windows into
  `windows100bp.gc_content.bot10p.bed`.
  4. Sort by position (`sort -k1,1 -k2,2n`) and gzip `windows100bp.gc_content.bot10p.bed` to create your bottom 10%
  GC-content QC file.
