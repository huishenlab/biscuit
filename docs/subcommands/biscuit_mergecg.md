---
title: biscuit mergecg
parent: BISCUIT Subcommands
nav_order: 9
permalink: /biscuit_mergecg/
---

## biscuit mergecg
```bash

Usage: biscuit mergecg [options] <ref.fa> <in.bed>

Options:
    -N        NOMe-seq mode, only merge C,G both in HCGD context
    -k INT    Minimum depth after merging - applies to the maximum depth
                  across samples [0]
    -h        This help

Note, in.bed is a position sorted bed file with beta values and coverages found
    in columns 4 and 5, respectively. Additional beta value-coverage column
    pairs are added for each additional sample. This is the format that would be
    found in the output of biscuit vcf2bed without the '-e' flag included.

```
