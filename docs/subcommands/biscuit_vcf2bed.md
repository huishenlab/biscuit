---
title: biscuit vcf2bed
parent: BISCUIT Subcommands
nav_order: 8
permalink: /biscuit_vcf2bed/
---

# biscuit vcf2bed
```bash

Usage: biscuit vcf2bed [options] <in.vcf>

Options:
    -t STR    Extract type {c, cg, ch, hcg, gch, snp} [cg]
    -k INT    Minimum coverage [3]
    -s STR    Sample, (takes "FIRST", "LAST", "ALL", or specific
                  sample names separated by ",") [FIRST]
    -e        Show context (reference base, context group {CG,CHG,CHH},
                  2-base {CA,CC,CG,CT} and 5-base context) before beta
                  value and coverage column
    -h        This help

```
