---
title: BISCUIT Subcommands
nav_order: 10
has_children: true
permalink: docs/subcommands
---

# BISCUIT Subcommands

Usage for available BISCUIT subcommands (printed with either `biscuit` or `biscuit help`).

```bash

Program: BISCUIT (BISulfite-seq CUI Toolkit)
Version: 1.6.1
Contact: Jacob Morrison <jacob.morrison@vai.org>

Usage: biscuit <command> [options]

Command:
 -- Read mapping
    index        Index reference genome sequences in the FASTA format
    align        Align bisulfite treated short reads using adapted BWA-mem
                     algorithm

 -- BAM operation
    tview        Text alignment viewer with bisulfite coloring
    bsstrand     Validate/correct bisulfite conversion strand label (YD tag)
    bsconv       Summarize/filter reads by bisulfite conversion (ZN tag)
    cinread      Print cytosine-read pair in a long form

 -- Base summary
    pileup       Pileup cytosine and mutations
    vcf2bed      Convert VCF to BED file
    mergecg      Merge C and G in CpG context

 -- Epireads
    epiread      Convert BAM to epibed format
    rectangle    Convert epiread format to rectangle format
    asm          Test allele-specific methylation

 -- Other
    version      Print BISCUIT and library versions
    help         Print usage and exit
    qc           Generate QC files from BAM
    bc           Extract barcodes from FASTQ files

```
