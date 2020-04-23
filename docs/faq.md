---
title: Frequently Asked Questions
nav_order: 9
---

# Frequently Asked Questions

## General Questions

*I am new to DNA methylation research and do not understand a lot of the
terminology. Are there sites that describe relevant terms and phrases?*

> One place to look is the
[Glossary]({{ site.baseurl }}{% link docs/glossary.md %}) page.

*I downloaded one of the pre-compiled binaries for BISCUIT, but want to run QC
on my data. Where can I find the QC script?*

> You can downloaded the `QC.sh` script from the
[BISCUIT release page](https://github.com/huishenlab/biscuit/releases/latest).
Supplementary QC asset files will be needed to run the script and can be found
for the hg19, hg39, and mm10 genomes alongside the QC script.

*I have aligned my data to a reference genome that does not have QC asset files
provided. What do I files do I need to run `QC.sh`?

> As of the latest release, you will need three BED files. You will need a
gzipped BED file with the locations of all CpGs in your genome (called
`cpg.bed.gz`). You will also need to break your genome into 100bp
non-overlapping windows, calculate the GC-content for each window, then find the
windows with the top and bottom 10% of GC-content. These regions will be placed
into two separate gzipped BED files, `windows100bp.gc_content.bot10p.bed.gz` for
the bottom 10% and `windows100bp.gc_content.top10p.bed.gz` for the top 10%.
Make sure to sort each BED file using (`sort -k1,1 -k2,2n bedfile.bed`).

## Installation Questions

*I ran `git clone git@github.com:huishenlab/biscuit.git`, but when trying to
compile the executable, I get an error message.*

> First, make sure you included `--recursive` when running `git clone`. If that
does not solve your problem, make sure you have `zlib` and `ncurses` in your
PATH.
