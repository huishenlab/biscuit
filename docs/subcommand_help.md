---
title: BISCUIT Subcommands
nav_order: 9
---

# BISCUIT Subcommands

Usage for all BISCUIT subcommands as of BISCUIT version 1.0.0.

## biscuit
```bash

Program: BISCUIT (BISulfite-seq CUI Toolkit)
Version: 1.0.0
Contact: Jacob Morrison <jacob.morrison@vai.org>

Usage: biscuit <command> [options]

Command:
 -- Read mapping
    index        Index reference genome sequences in the FASTA format
    align        Align bisulfite treated short reads using adapted BWA-mem
                     algorithm

 -- BAM operation
    tview        Text alignment viewer with bisulfite coloring
    bsstrand     Validate/correct bisulfite conversion strand label (YD tag)
    bsconv       Summarize/filter reads by bisulfite conversion (ZN tag)
    cinread      Print cytosine-read pair in a long form

 -- Base summary
    pileup       Pileup cytosine and mutations
    vcf2bed      Convert VCF to BED file
    mergecg      Merge C and G in CpG context

 -- Epireads
    epiread      Convert BAM to epibed format
    rectangle    Convert epiread format to rectangle format
    asm          Test allele-specific methylation

 -- Other
    version      Print BISCUIT and library versions
    qc           Generate QC files from BAM

```

## biscuit index
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

## biscuit align
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

## biscuit tview
```bash

Usage: biscuit tview [options] <in.bam> <ref.fa>

Options:
    -g STR    Go directly to this position
    -m INT    Max number of reads to load per position [50]
    -n STR    Highlight the read(s) with STR as the read name
    -f INT    Flanking sequence length [100]
    -h        This help

```

## biscuit bsstrand
```bash

Usage: biscuit bsstrand [options] <ref.fa> <in.bam> [out.bam]

Options:
    -g STR    Region (optional, will process the whole bam if not specified)
    -y        Append count of C>T (YC tag) and G>A (YG tag) in out.bam
    -c        Correct bsstrand in out.bam, YD tag will be replaced if it exists
                  and created if not
    -h        This help

```

## biscuit bsconv
```bash

Usage: biscuit bsconv [options] <ref.fa> <in.bam> [out.bam]

Options:
    -g STR      Region (optional, will process the whole bam if not specified)
    -u          Filter unclear bs-strand (YD:u) reads
    -m INT      Filter: maximum CpH retention [Inf]
    -f FLOAT    Filter: maximum CpH retention fraction [1.0]
    -y FLOAT    Filter: maximum CpY retention fraction [1.0]
    -a INT      Filter: maximum CpA retention [Inf]
    -c INT      Filter: maximum CpC retention [Inf]
    -t INT      Filter: maximum CpT retention [Inf]
    -p          Print in tab-separated format, print order:
                    CpA_R, CpA_C, CpC_R, CpC_C, CpG_R, CpG_C, CpT_R, CpT_C
    -v          Show filtered reads instead of remaining reads
    -h          This help

```

## biscuit cinread
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

## biscuit pileup
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
    -P FLOAT    Prior probability for heterozygous variant [0.33333]
    -Q FLOAT    Prior probability for homozygous variant [0.33333]
    -h          This help

```

## biscuit vcf2bed
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

## biscuit epiread
```bash

Usage: biscuit epiread [options] <ref.fa> <in.bam>

Options:
    -B STR    Bed input for SNP display in epiread output
    -g STR    Region (optional, will process the whole bam if not specified)
    -s STR    Step of window dispatching [100000]
    -@ INT    Number of threads [3]

Output options:
    -o STR    Output file [stdout]
    -N        NOMe-seq mode [off]
    -P        Pairwise mode [off]
    -O        Old BISCUIT epiread format, not compatible with -P [off]
    -A        Print all CpG and SNP locations in location column, ignored if -O not given [off]
    -v        Verbose (print additional info for diagnostics) [off]

Filter options:
    -b INT    Minimum base quality [20]
    -m INT    Minimum mapping quality [40]
    -a INT    Minimum alignment score (from AS-tag) [40]
    -t INT    Max cytosine retention in a read [999999]
    -l INT    Minimum read length [10]
    -5 INT    Minimum distance to 5' end of a read [3]
    -3 INT    Minimum distance to 3' end of a read [3]
    -c        NO filtering secondary mapping
    -d        Double count cytosines in overlapping mate reads (avoided
                  by default)
    -u        NO filtering of duplicate
    -p        NO filtering of improper pair
    -n INT    Maximum NM tag [999999]
    -h        This help

```

## biscuit rectangle
```bash

Usage: biscuit rectangle [options] <ref.fa> <in.epiread>

Options:
    -o STR    Output file [stdout]
    -h        This help
Note, this is not currently compatible with epiread format when run with the -A flag
    i.e., biscuit epiread -A [-B snps.bed] <ref.fa> <in.bam>

```

## biscuit asm
```bash

Usage: biscuit asm [options] <in.epiread>

Options:
    -h    This help

```

## biscuit qc
```bash

Usage: biscuit qc [options] <ref.fa> <in.bam> <sample_name>

Options:
    -s    Run for single-end data
    -h    This help

Note, this currently only produces a subset of QC metrics. Use scripts/QC.sh for full QC

```

## biscuit version

Running `biscuit version` prints out the BISCUIT version and the commit hashes for each of the git submodules used by
BISCUIT.
