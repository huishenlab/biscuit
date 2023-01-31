---
title: biscuit index
parent: BISCUIT Subcommands
nav_order: 1
permalink: /biscuit_index/
---

# biscuit index
```bash

Usage: biscuit index [options] <in.fasta>

Options:
    -a STR     BWT construction algorithm: bwtsw, div, or is [auto]
    -p STR     Prefix of the index [same as fasta name]
    -6         Index files named as <in.fasta>.64.* instead of <in.fasta>*
    -h         This help

Warning: '-a bwtsw' does not work for short genomes, while '-a is' and '-a div'
         do not work not for long genomes. Please choose '-a' according to the
         length of the genome.

```
