---
title: biscuit align
parent: BISCUIT Subcommands
nav_order: 2
permalink: /biscuit_align/
---

# biscuit align
```bash

Usage: biscuit align [options] <fai-index base> <in1.fq> [in2.fq]

Algorithm options:
    -@ INT          Number of threads [1]
    -b INT          For PE, map read 1 to parent, read 2 to daughter
                        (1: directional library) or read1 and read 2 to both
                        (0: non-directional library [default]). For SE, parent
                        (1), daughter (3), or both (0 [default]). Note, parent
                        is the bisulfite treated strand and the daughter is
                        the complement strand.
    -f INT          1: BSW strand; 3: BSC strand; 0: both. Note, libraries
                        targeting either BSW or BSC are unseen so far! [0]
    -k INT          Minimum seed length [19]
    -w INT          Band width for banded alignment [100]
    -d INT          Off-diagonal X-dropoff [100]
    -r FLOAT        Look for internal seeds inside a seed longer than
                        {-k}*FLOAT [1.5]
    -y INT          Seed occurrence for the 3rd round of seeding [20]
    -J STR          Adaptor of read 1 (fastq direction)
    -K STR          Adaptor of read 2 (fastq direction)
    -z INT          Minimum base quality to keep from both ends of reads [0]
    -5 INT          Number of extra bases to clip from 5'-end [0]
    -3 INT          Number of extra bases to clip from 3'-end [0]
    -c INT          Skip seeds with more than INT occurrences [500]
    -D FLOAT        Drop chains shorter than FLOAT fraction of the longest
                        overlapping chain [0.50]
    -W INT          Discard a chain if seeded bases shorter than INT [0]
    -m INT          Perform at most INT rounds of mate rescues for a read [50]
    -S              Skip mate rescue
    -P              Skip pairing - mate rescue performed unless -S also given
    -e              Discard full-length exact matches

Scoring options:
    -A INT          Score for a sequence match, scales options -TdBOELU unless
                        overridden [1]
    -B INT          Penalty for a mismatch [2]
    -O INT[,INT]    Gap open penalties for deletions and insertions [6,6]
    -E INT[,INT]    Gap extension penalty; a gap of size, g, has penalty
                        {-O} + {-E}*g [1,1]
    -L INT[,INT]    Penalty for 5'- and 3'-end clipping [10,10]
    -U INT          Penalty for an unpaired read pair [17]

Input/output options:
    -1 STR          Align a single read STR
    -2 STR          Align a read STR paired with -1 read
    -i              Turn off autoinference of ALT chromosomes
    -p              Smart pairing (ignores in2.fq)
    -R STR          Read group header line (such as '@RG\\tID:foo\\tSM:bar')
    -F              Suppress SAM header output
    -H STR/FILE     Insert STR to header if it starts with @ or insert lines
                        in FILE
    -j              Treat ALT contigs as part of the primary assembly (i.e.
                        ignore <fai-index base>.alt file)
    -q              Do not modify mapQ of supplementary alignments
    -T INT          Minimum score to output [30]
    -g INT[,INT]    Maximum number of hits output in XA [5,5]
    -a              Output all alignments for SE or unpaired PE
    -C              Append FASTA/FASTQ comment to SAM output
    -V              Output the reference FASTA header in the XR tag
    -Y              Use soft clipping for supplementary alignments
    -M              Mark shorter split hits as secondary
    -I FLOAT[,FLOAT[,INT[,INT]]]
                    Specify the mean, standard deviation (10%% of the mean
                        if absent), maximum (4 sigma from the mean if absent)
                        and minimum of insert size distribution. FR orientation
                        only [inferred]
    -v INT          Verbosity level: 
                        1: error, 2: warning, 3: message, 4+: debugging [3]
    -h              This help

```
