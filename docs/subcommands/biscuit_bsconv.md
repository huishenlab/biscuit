---
title: biscuit bsconv
parent: BISCUIT Subcommands
nav_order: 5
permalink: /biscuit_bsconv/
---

# biscuit bsconv
```bash

Usage: biscuit bsconv [options] <ref.fa> <in.bam> [out.bam]

Options:
    -g STR      Region (optional, will process the whole bam if not specified)
    -u          Filter unclear bs-strand (YD:u) reads
    -m INT      Filter: maximum CpH retention [Inf]
    -f FLOAT    Filter: maximum CpH retention fraction [1.0]
    -y FLOAT    Filter: maximum CpY retention fraction [1.0]
    -a INT      Filter: maximum CpA retention [Inf]
    -c INT      Filter: maximum CpC retention [Inf]
    -t INT      Filter: maximum CpT retention [Inf]
    -x INT      Filter: maximum CpY retention [Inf]
    -p          Print in tab-separated format, print order:
                    CpA_R, CpA_C, CpC_R, CpC_C, CpG_R, CpG_C, CpT_R, CpT_C
    -v          Show filtered reads instead of remaining reads
    -h          This help

```
