---
title: Quality Control
nav_order: 3
parent: Read Mapping
---
# Quality Control

A bash script is provided to simplify QC procedure. It generates QC information
that can be picked up by MultiQC.
```bash
$ ./scripts/QC.sh -v input.vcf setup_file sample.name input.bam
```
Setup files for [hg19](http://zwdzwd.io/BISCUITqc/hg19_QC_assets.zip),
[hg38](http://zwdzwd.io/BISCUITqc/hg38_QC_assets.zip) and
[mm10](http://zwdzwd.io/BISCUITqc/mm10_QC_assets.zip) are provided. Other genome
builds can be built following the same format.

## Validate Bisulfite Conversion Label

Sometimes, the bisulfite conversion labels in a given alignment are inaccurate,
conflicting, or ambiguous. The `bsstrand` subcommand summarizes these labels, given
the number of C&#8594;T and G&#8594;A substitutions. As an option, it can also
correct inaccurate labels.
```bash
$ biscuit bsstrand -g "chr1:1000000-1050000" my_reference.fa input.bam
```
returns something like
```
Mapped reads: 2865
Unmapped reads: 44
Corrected reads: 0 (0.00%)

Strand Distribution:
strand\BS      BSW (f)      BSC (r)
     R1 (f):   67           528
     R1 (r):   692          47
     R2 (f):   765          92
     R2 (r):   107          567


R1 mapped to converted:   114
R1 mapped to synthesized: 1220
R2 mapped to converted:   1332
R2 mapped to synthesized: 199

Confusion counts:
orig\infer      BSW (f)      BSC (r)      Conflict (c) Unknown (u)
     BSW (f):   1483         48           11           89
     BSC (r):   27           1084         6            117
Conflict (c):   0            0            0            0
 Unknown (u):   0            0            0            0
```

The inferred `YD` tag gives the following
  - f: foward/Waston strand
  - r: reverse/Crick strand
  - c: conflicting strand information
  - u: unintelligible strand source (unknown)

`YD` is inferred based on the count of `C>T` (`nCT`) and `G>A` (`nGA`) observations
in each read according to the following rules:

  - If both `nCT` and `nGA` are zero, `YD = u` and `s = min(nGA,nCT) / max(nGA,nCT)`.
  - If `nCT > nGA` and (`nGA == 0` or `s <= 0.5`), then `YD = f`.
  - If `nCT < nGA` and (`nCT == 0` or `s <= 0.5`), then `YD = r`.
  - All other scenarios give `YD = c`.

The flag `-y` appends `nCT` (YC tag) and `nGA` (YG tag) in the output BAM file.

### When `-b 1` was specified in mapping

`-b 1` forces the mapping to be conducted in a stranded manner. Using `bsstrand`,
you can validate whether this enforcement is successful. For example,
```bash
$ biscuit bsstrand -g "chr1:1000000-1050000" my_reference.fa stranded.bam
```
returns
```
Mapped reads: 2918
Unmapped reads: 93
Corrected reads: 0 (0.00%)

Strand Distribution:
strand\BS      BSW (f)      BSC (r)
     R1 (f):   840          0
     R1 (r):   0            551
     R2 (f):   0            565
     R2 (r):   962          0


R1 mapped to converted:   1391
R1 mapped to synthesized: 0
R2 mapped to converted:   0
R2 mapped to synthesized: 1527

Confusion counts:
orig\infer      BSW (f)      BSC (r)      Conflict (c) Unknown (u)
     BSW (f):   715          916          9            162
     BSC (r):   499          500          12           105
Conflict (c):   0            0            0            0
 Unknown (u):   0            0            0            0
```
In this case, read 1 is always mapped to the converted strand and read 2 is
always mapped to the synthesized strand.

## Summarize and Filter Reads by Bisulfite Conversion

For some library preparations, incomplete conversions are enriched in a subset
of reads that need to be filtered. The `bsconv` subcommand transforms the BAM
file into a BAM that contains the `ZN` tag (like `ZN:Z:CA_R0C11,CC_R1C14,CG_R0C2,CT_R1C5`).
This tag summarizes counts of retention and conversion for four different
cytosine contexts, `CpA`, `CpC`, `CpG` and `CpT`. In general, it contains a
minimum threshold of `CpA`, `CpC`, `CpT` or `CpH`. The `-b` option outputs the
summary in tables instead of as tags in the BAM file.
```bash
$ biscuit bsconv -g "chr1:1000000-1050000" my_reference.fa input.bam
```
