---
title: biscuit cinread
parent: BISCUIT Subcommands
nav_order: 6
permalink: /biscuit_cinread/
---

# biscuit cinread
```bash

Usage: biscuit cinread [options] <ref.fa> <in.bam>

Options:
    -g STR    Region (optional, will process the whole bam if not specified)
    -t STR    Target (c, cg, ch, hcg, gch, hch) [cg]
    -p STR    Content to print, ","-delimited:
                  QNAME, QPAIR, STRAND, BSSTRAND, MAPQ
                  QBEG, QEND, CHRM, CRPOS, CGRPOS
                  CQPOS, CRBASE, CCTXT, CQBASE, CRETENTION
                      [QNAME,QPAIR,BSSTRAND,CRBASE,CQBASE]
    -s        Consider secondary mapping [off]
    -o STR    Output file [stdout]
    -h        This help

```
