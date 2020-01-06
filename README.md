# BISCUIT [![Travis-CI Build Status](https://travis-ci.org/huishenlab/biscuit.svg?branch=master)](https://travis-ci.org/huishenlab/biscuit) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.48262.svg)](http://dx.doi.org/10.5281/zenodo.48262)
<!-- [![Downloads](https://img.shields.io/github/downloads/zwdzwd/biscuit/latest/total.svg)](https://github.com/zwdzwd/biscuit/releases) -->

BISulfite-seq CUI Toolkit (BISCUIT) is a utility for analyzing sodium bisulfite
conversion-based DNA methylation/modification data. It was written to perform
alignment, DNA methylation and mutation calling, and allele specific methylation
from bisulfite sequencing data.

User Guide is available at
[https://huishenlab.github.io/biscuit/](https://huishenlab.github.io/biscuit/).

# Download

All releases are available [here](https://github.com/huishenlab/biscuit/releases).

## Download Source Code and Compile

### Download Source Code

You can compile from source code using either `git` or `curl`.

Using `git`,
```bash
$ git clone --recursive git@github.com:huishenlab/biscuit.git
```
Note, after v0.2.0, if you choose to download via Git, make sure to use
`git clone --recursive` to retrieve the submodules.

Using `curl`,
```bash
$ curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest | 
    grep browser_download_url | grep release-source.zip | cut -d '"' -f 4)
```

### Compile

If you download the source code, BISCUIT can be compiled through,

From `git`,
```bash
$ cd biscuit
$ make
```

From `curl`
```bash
$ unzip release.zip
$ cd biscuit-release
$ make
```

The created `biscuit` binary is the main entry point.

## Download Precompiled Binaries
The precompiled binaries can be found on the
[latest release page](https://github.com/huishenlab/biscuit/releases/latest)
on Github (currently only supports latest versions of Linux and macOS). You can
also do this in terminal using the following one-liner:

For macOS,
```bash
$ curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep darwin | cut -d '"' -f 4)
```

For linux,
```bash
$ curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep linux_amd64 | cut -d '"' -f 4)
$ chmod +x biscuit_*
```

# Usage

## Index Reference for Alignment

```bash
biscuit index my_reference.fa
```
The index of BISCUIT is composed of the 2-bit packed reference (`*.bis.pac`,
`*.bis.amb`, `*.bis.ann`). The FM-index and suffix array of the parent strand
(`*.par.bwt` and `*.par.sa`) and the daughter strand (`*.dau.bwt` and `*.dau.sa`).
Here, `*` refers to the input reference file, `my_reference.fa`.

For more information regarding read mapping, run `biscuit index` in the terminal.

## Read Alignment

After creating an index of the reference genome, the sequenced reads can be
mapped. In conjunction with [samtools](https://github.com/samtools/samtools)
and [samblaster](https://github.com/GregoryFaust/samblaster), BISCUIT can be
used to map, mark duplicate, and sort the provided reads.  Additionally,
`samtools` can be used to create an indexed alignment BAM file.

```bash
$ biscuit align -M my_reference.fa read1.fq.gz read2.fq.gz |
    samblaster -M | samtools sort -o my_output.bam -O BAM -
$ samtools index my_output.bam
```
The first line is referred to as the *biscuitBlaster* pipeline and provides
the user with a sorted BAM file of aligned reads with duplicates marked, but not
removed. During the [Read Pileup](https://huishenlab.github.io/biscuit/docs/pileup.html)
state, BISCUIT will skip marked duplicates by default. If desired, `samblaster`
has a flag for removing duplicates during the marking process. Or, if desired,
there is a flag for BISCUIT to retain marked duplicates during the pileup stage.

For more information regarding read mapping, run `biscuit align` in the terminal.

## Visualize Alignment

The `tview` subroutine allows you to visualize mappings generated by BISCUIT.
`biscuit tview` is similar to `samtools tview`, but includes additional coloring
for reads that have undergone bisulfite conversion. One difference is that a 
reference FASTA file is required for BISCUIT so that bisulfite conversion can be
identified. Therefore, make sure you include the reference FASTA in the command.
```bash
$ biscuit tview -g chr19:7525080 my_input.bam my_reference.fa
```

An example of the output can be found
[here](https://huishenlab.github.io/biscuit/docs/alignment/visualization.html).

In the reference, CpGs are marked in red, while C's and G's in other contexts are
marked in blue. For reads, retained C's/G's are marked in red, converted C's/G's
are marked in blue. Mismatches not related to bisulfite conversion are marked in
yellow.

## Mark Duplicate Reads

This functionality is DEPRECATED and will be removed in the next release!!!!
We suggest using samblaster and the
[biscuitBlaster](https://huishenlab.github.io/biscuit/bamoperation/) pipeline
for duplicate marking. This step is optional. The mark duplicate of BISCUIT is
bisulfite strand aware.
```bash
$ biscuit markdup input.bam output.bam
```

## Pileup

Like samtools, BISCUIT can extract DNA methylation as well as genetic information.
The following shows how to produce a tabix-indexed VCF file:
```bash
$ biscuit pileup /path/to/my_reference.fa my_output.bam -o my_pileup.vcf
$ bgzip output.vcf
$ tabix -p vcf output.vcf.gz
```

## Make BED Files

After generating the VCF file with both genetic and methylation information,
beta values and coverage can be extracted to study the methylation levels at
any CpGs sequenced. The following shows how to extract this information from
the VCF file:

```bash
$ biscuit vcf2bed -t cg my_pileup.vcf.gz
```

This will create a BED file (`my_pileup.bed`), which can be used for further
analysis.

The `-t` flag can be used to retrieve mutation and other information, including

  * `snp` for SNP information
  * `c` for all cytosines
  * `hcg` for HCG (from NOMe-seq)
  * `gch` for GCH (from NOMe-seq)
          
Information about the sequence context for each row can be included by adding
the `-e` flag to the command, while filtering out low coverage rows can be done
using the `-k` flag.

## QC Your Run

A bash script is provided to simplify QC procedure. It generates QC information
that can be picked up by MultiQC.
```bash
$ ./scripts/QC.sh -v input.vcf assets_directory genome.fa sample.name input.bam
```
Setup files for [hg19](http://zwdzwd.io/BISCUITqc/hg19_QC_assets.zip),
[hg38](http://zwdzwd.io/BISCUITqc/hg38_QC_assets.zip) and
[mm10](http://zwdzwd.io/BISCUITqc/mm10_QC_assets.zip) are provided. Other genome
builds can be built following the same format.

## Epireads and Allele-specific Methylation

The following illustrates how to produce an epiread file, which carries
information about the epi-haplotype.
```bash
$ biscuit epiread -r my_reference.fa -i input.bam -B snp.bed
```

To test all SNP-CpG pairs,
```bash
$ biscuit epiread -r my_reference.fa -P -i input.bam -B snp.bed
```
More details can be found [here](https://huishenlab.github.io/biscuit/epiread_format/).

Allele-specific methylation can be extracted from the epiread files via:
```bash
sort -k1,1 -k2,2n -k3,3n in.epiread > out.epiread
biscuit asm out.epiread > out.asm
```

## Validate Bisulfite Conversion Label

Sometimes, the bisulfite conversion labels in a given alignment are inaccurate,
conflicting, or ambiguous. The `bsstrand` command summarizes these labels, given
the number of C&#8594;T and G&#8594;A substitutions. As an option, it can also
correct inaccurate labels.
```bash
$ biscuit bsstrand -g "chr1:1000000-1050000" my_reference.fa input.bam 
```
returns something like
```
Mapped reads: 16688
Unmapped reads: 29
Corrected reads: 0 (0.00%)

Confusion counts:
orig\infer     BSW (f)      BSC (r)      Conflict (c) Unknown (u)
     BSW (f):   8426         42           4            12
     BSC (r):   15           8167         3            19
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

# Acknowledgements

 * lib/aln was adapted from Heng Li's BWA-mem code.
 * lib/htslib was subtree-ed from the htslib library.
 * lib/klib was subtree-ed from Heng Li's klib.
