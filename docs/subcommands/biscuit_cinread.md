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
    -p STR    Content to print, ","-delimited: [QNAME,QPAIR,BSSTRAND,CRBASE,CQBASE]
                  Options + Descriptions:
                      QNAME:      read name
                      QPAIR:      which read in pair
                      STRAND:     forward or reverse strand
                      BSSTRAND:   which original strand the read derives from
                      MAPQ:       MAPQ score
                      QBEG:       read start position
                      QEND:       read end position
                      CHRM:       chromosome
                      CRPOS:      cytosine position on reference
                      CGRPOS:     CpG position on reference (-1 if not applicable)
                      CQPOS:      cytosine position on read
                      CRBASE:     cytosine reference base
                      CCTXT:      cytosine context, strand flipped
                      CQBASE:     base called on read
                      CRETENTION: retention (R) or conversion (C))
    -s        Consider secondary mapping [off]
    -o STR    Output file [stdout]
    -h        This help

```
