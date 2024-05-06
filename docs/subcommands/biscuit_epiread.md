---
title: biscuit epiread
parent: BISCUIT Subcommands
nav_order: 10
permalink: /biscuit_epiread/
---

# biscuit epiread
```bash

Usage: biscuit epiread [options] <ref.fa> <in.bam>

Options:
    -B STR    Bed input for SNP display in epiread output
    -g STR    Region (optional, will process the whole bam if not specified)
    -s STR    Step of window dispatching [100000]
    -@ INT    Number of threads [3]

Output options:
    -o STR    Output file [stdout]
    -N        NOMe-seq mode [off]
    -L INT    maximum read length (will need to be increased for long reads) [302]
    -M        BAM file has modBAM tags (MM/ML) [off]
    -P        Pairwise mode [off]
    -O        Old BISCUIT epiread format, not compatible with -P [off]
    -A        Print all CpG and SNP locations in location column, ignored if -O not given [off]
    -v        Verbose (print additional info for diagnostics) [off]

Filter options:
    -b INT    Minimum base quality [20]
    -m INT    Minimum mapping quality [40]
    -a INT    Minimum alignment score (from AS-tag) [40]
    -t INT    Max cytosine retention in a read [999999]
    -l INT    Minimum read length [10]
    -5 INT    Minimum distance to 5' end of a read [3]
    -3 INT    Minimum distance to 3' end of a read [3]
    -E        NO filtering of empty epireads
    -d        Double count cytosines in overlapping mate reads (avoided
                  by default)
    -u        NO filtering of duplicate
    -p        NO filtering of improper pair
    -n INT    Maximum NM tag [999999]
    -y FLT    Minimum probability a modification is correct (0.0 - 1.0) [0.900000]
    -h        This help

Note, the -O (old epiread format) and -P (pairwise format for biscuit asm) are not guaranteed
    to match output from biscuit pileup. These file formats have been left in for legacy purposes.
    Default output with an unfiltered BISCUIT SNP BED file (biscuit pileup ...
    -> biscuit vcf2bed -t snp ...) should have the same results in the epiBED as in pileup.

```
