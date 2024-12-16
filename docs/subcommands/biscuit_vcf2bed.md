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
    -t STR    Extract type {c, cg, ch, hcg, gch, snp} [CG]
    -k INT    Minimum coverage (see Note 1) [1]
    -s STR    Sample, (takes "FIRST", "LAST", "ALL", or specific
                  sample names separated by ",") [FIRST]
    -e        Show context (reference base, context group {CG,CHG,CHH},
                  2-base {CA,CC,CG,CT} and 5-base context) before beta
                  value and coverage column
    -c        Output Beta-M-U instead of Beta-Cov.
    -h        This help

Note 1: Starting with version 1.6.0, the default minimum coverage was changed from
        three (3) to one (1) to better match other tools and serve as a better default
        for single-cell experiments.

```
