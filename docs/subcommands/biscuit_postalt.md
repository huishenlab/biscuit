---
title: biscuit postalt
parent: BISCUIT Subcommands
nav_order: 6
permalink: /biscuit_postalt/
---

# biscuit postalt
```bash

Usage: biscuit postalt [options] <in.bam> [out.bam]

Options:
    --unplaced    Enable unplaced/alt-contig demotion (the only mode)
    -p STR        Primary-contig regex [^chr([0-9]{1,2}|X|Y|M|MT)$]
    -d INT        NM near-score delta (a primary hit within
                      bestNM+delta is a promotion candidate) [1]
    -q INT        MAPQ to set when unique on the primary assembly [60]
    -g STR        Process only this region (needs an index)
    -h            This help


NOTE: in.bam/out.bam may be '-' for stdin/stdout; input may be SAM or BAM.

```
