---
title: biscuit pileup
parent: BISCUIT Subcommands
nav_order: 7
permalink: /biscuit_pileup/
---

# biscuit pileup
```bash

Usage: biscuit pileup [options] <ref.fa> <in1.bam> [in2.bam in3.bam ...]
Som. Mode Usage: biscuit pileup [options] <-S -T tum.bam -I norm.bam> <ref.fa>

Options:
    -g STR      Region (optional, will process the whole bam if not specified)
    -@ INT      Number of threads [3]
    -s INT      Step of window dispatching [100000]
    -N          NOMe-seq mode [off]
    -S          Somatic mode, must provide -T and -I arguments [off]
    -T STR      Somatic mode, tumor BAM
    -I STR      Somatic mode, normal BAM

Output options:
    -o STR      Output file [stdout]
    -w STR      Pileup statistics output prefix [same as output]
    -v INT      Verbosity level (0: no added info printed, 0<INT<=5: print
                    diagnostic info, INT>5: print diagnostic and debug info) [0]

Filter options:
    -b INT      Minimum base quality [20]
    -m INT      Minimum mapping quality [40]
    -a INT      Minimum alignment score (from AS-tag) [40]
    -t INT      Maximum cytosine retention in a read [999999]
    -l INT      Minimum read length [10]
    -5 INT      Minimum distance to 5' end of a read [3]
    -3 INT      Minimum distance to 3' end of a read [3]
    -r          NO redistribution of ambiguous (Y/R) calls in SNP genotyping
    -c          NO filtering secondary mapping
    -d          Double count cytosines in overlapping mate reads (avoided
                    by default)
    -u          NO filtering of duplicate flagged reads
    -p          NO filtering of improper pair flagged reads
    -n INT      Maximum NM tag [999999]

Genotyping options:
    -E FLOAT    Error rate [0.001]
    -M FLOAT    Mutation rate [0.001]
    -x FLOAT    Somatic mutation rate [0.001]
    -C FLOAT    Contamination rate [0.010]
    -P FLOAT    Prior probability for heterozygous variant [0.333]
    -Q FLOAT    Prior probability for homozygous variant [0.333]
    -h          This help
```
