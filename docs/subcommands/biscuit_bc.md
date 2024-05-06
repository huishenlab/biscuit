---
title: biscuit bc
parent: BISCUIT Subcommands
nav_order: 14
permalink: /biscuit_bc/
---

# biscuit bc
```bash

Usage: biscuit bc [options] <FASTQ 1> [FASTQ 2]

Processing Options:
    -m, --mate INT             which mate the barcode is in (1 or 2) [1]
    -s, --bc-start INT         start position of barcode in read (1-based) [1]
    -l, --bc-length INT        length of barcode [8]
Output Options:
    -o, --output-prefix STR    prefix for output files (NULL writes to stdout) [NULL]
General Options:
    -h, --help             This help

Note 1: When writing to stdout, reads 1 and 2 will alternate (i.e., are interleaved)
Note 2: Also adds an artificial UMI (AAAAAAAA) for compatibility purposes and to serve
        as a placeholder for future UMI work

```
