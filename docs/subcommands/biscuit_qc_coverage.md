---
title: biscuit qc_coverage
parent: BISCUIT Subcommands
nav_order: 14
permalink: /biscuit_qc_coverage/
---

# biscuit qc_coverage
```bash

Usage: biscuit qc_coverage <ref.fa> <cpgs.bed.gz> <in.bam>

Options:
    -P STR    Prefix for output file names
    -B STR    Bottom 10 percent GC content windows BED file
    -T STR    Top 10 percent GC content windows BED file
    -s INT    Step size of windows [100000]
    -@ INT    Number of threads [3]
Filter options:
    -b INT    Minimum base quality [20]
    -a INT    Minimum alignment score (from AS-tag) [40]
    -t INT    Max cytosine retention in a read [999999]
    -l INT    Minimum read length [10]
    -u        NO filtering of duplicate reads
    -p        NO filtering of improper pair
    -n INT    Maximum NM tag [999999]
    -h        Print usage

```
