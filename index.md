---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults
title: Home
layout: default
description: "BISulfite-seq CUI Toolkit (BISCUIT) BISulfite-seq CUI Toolkit
  (BISCUIT) is a utility suite for analyzing sodium bisulfite
  conversion-based DNA methylation/modification data. It was written
  to perform alignment, DNA methylation and mutation calling, allele
  specific methylation from bisulfite sequencing data."
nav_order: 1
---

# BISCUIT - Understand Sequencing Data with Bisulfite Conversion
{: .fs-9 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/huishenlab/biscuit){: .btn .fs-5 .mb-4 .mb-md-0 }

---

BISulfite-seq CUI Toolkit (BISCUIT) is a utility suite for analyzing bulk and single-cell sodium bisulfite- or
enzyme-based DNA methylation/modification data, such as WGBS, capture bisulfite sequencing, RRBS, NOMe-seq, and EM-seq.
It was written to perform read alignment, DNA methylation and mutation calling, and allele specific methylation from
bisulfite or bisulfite-like sequencing data.

BISCUIT was developed by Wanding Zhou while he was a member of the [Shen Lab](https://shenlab.vai.org) at Van Andel
Institute. He now holds a faculty position at University of Pennsylvania and Children's Hospital of Philadelphia.
BISCUIT is currently maintained by Jacob Morrison (who also developed the User's Guide website) in the Shen Lab. Current
versions of BISCUIT are available at [https://github.com/huishenlab/biscuit](https://github.com/huishenlab/biscuit),
while legacy versions are located at [https://github.com/zhou-lab/biscuit](https://github.com/zhou-lab/biscuit).

## Quick Start

In order to get started with performing analyses with BISCUIT, precompiled binaries are available for download on the
[BISCUIT release page](https://github.com/huishenlab/biscuit/releases/latest). Note, binaries are only available for
Linux and macOS. (See [Download and Install](#download-and-install) for more information about downloading and
installing BISCUIT).

The basic workflow to align and extract methylation information using BISCUIT is:
1. Create an index of the reference genome (only needs to be done once for each reference).
2. Align sequencing reads to the reference.
3. Create a pileup VCF of DNA methylation and genetic information.
4. Extract DNA methylation into BED format.

Practically, the commands to run are:
```bash
# Create index of the reference genome (only needs to be run once for each reference)
# Gzipped FASTA references can also be used
biscuit index my_reference.fa

# Align sequencing reads to the reference
# Gzipped FASTQ files can also be used
biscuit align -@ NTHREADS -R "my_rg" /path/to/my_reference.fa read1.fastq read2.fastq |
    dupsifter /path/to/my_reference.fa | samtools sort -@ NTHREADS -o my_output.bam -O BAM -
samtools index my_output.bam

# Create a pileup VCF of DNA methylation and genetic information
# Also compresses and indexes the VCF
biscuit pileup -@ NTHREADS -o my_pileup.vcf /path/to/my_reference.fa my_output.bam
bgzip -@ NTHREADS my_pileup.vcf
tabix -p vcf my_pileup.vcf.gz

# Extract DNA methylation into BED format
# Also compresses and indexes the BED
biscuit vcf2bed my_pileup.vcf.gz > my_methylation_data.bed
bgzip my_methylation_data.bed
tabix -p bed my_methylation_data.bed.gz
```
This basic order of commands will produce all the necessary files needed to read data into R using the R/Bioconductor
companion package, [biscuiteer](https://www.bioconductor.org/packages/release/bioc/html/biscuiteer.html).

An overview of all available functionalities can be found below in the
[Overview of Functionalities](#overview-of-functionalities) section.

## Download and Install

BISCUIT is available as a [precompiled binary](#download-precompiled-binaries) (for macOS and Linux), as
[source code](#download-source-code-and-compile) for compilation on your own machine, as a
[conda recipe](#download-with-conda), or as a [Docker container](#download-the-docker-container).

### Download Precompiled Binaries

Precompiled binaries can be found on the [latest release page](https://github.com/huishenlab/biscuit/releases/latest) on
GitHub. Currently, there are only precompiled binaries for the latest versions of Linux and macOS. You can also download
the binaries directly from the terminal using the following one-liner:

On macOS,
```bash
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep darwin_amd64 | cut -d '"' -f 4)
mv biscuit_* biscuit
chmod +x biscuit
```

On Linux,
```bash
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep linux_amd64 | cut -d '"' -f 4)
mv biscuit_* biscuit
chmod +x biscuit
```

To download the scripts to generate the QC asset files, generate QC files, and flip PBAT strands post-alignment, run
```bash
# QC asset build
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep build_biscuit_QC_assets.pl | cut -d '"' -f 4

# QC bash script
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep QC.sh | cut -d '"' -f 4

# Flip PBAT strands script
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep flip_pbat_strands.sh | cut -d '"' -f 4
```
These commands work on both macOS and Linux.

### Download Source Code and Compile

The source code for BISCUIT can be downloaded using either `git` or `curl`. Compilation requires that `zlib` and
`ncurses` are in the PATH environment variable.

Using `git`,
```bash
git clone --recursive git@github.com:huishenlab/biscuit.git
cd biscuit
make
```
Note, after v0.2.0, if downloading via `git`, make sure to use the `--recursive` flag to get the submodules. If an SSH
key has not been set up, and you receive a "permission denied" error, replace the first line with
```bash
git clone --recursive https://github.com/huishenlab/biscuit.git
```

Using `curl`,
```bash
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep release-source.zip | cut -d '"' -f 4)
unzip release-source.zip
cd biscuit-release
make
```

The QC and strand-flipping scripts can be found in the `scripts/` directory.

### Download with Conda

Note, this requires that `conda` has been installed. To download with conda, run:
```bash
conda install -c bioconda biscuit
```

This will also install both `QC.sh` and `build_biscuit_QC_assets.pl`.

### Download the Docker Container

The Docker container can be downloaded from [GitHub](https://github.com/huishenlab/sv_calling_docker) via:
```bash
git clone git@github.com:huishenlab/sv_calling_docker.git
```

For more information about the docker container, see
[Structural Variant Calling]({{ site.baseurl }}{% link docs/structural_variants.md %}).

## Overview of Functionalities

The following list provides an overview of the different subcommands and the various functionalities provided by
`biscuit`. You can also find much of this by typing `biscuit` in the terminal. Help for each subcommand can be found on
the [BISCUIT Subcommands]({{ site.baseurl }}{% link docs/subcommands/subcommand_help.md %}) page or by typing
`biscuit (subcommand)` in the terminal.

### Read Mapping

  - `index` Index reference genome (see 
  [Read Mapping]({{ site.baseurl }}{% link docs/alignment/alignment.md %}))
  - `align` Map bisulfite converted short reads to reference (see
  [Read Mapping]({{ site.baseurl }}{% link docs/alignment/alignment.md %}))

### BAM Operation

  - `tview` View read mapping in terminal with bisulfite coloring (see
  [Visualization]({{ site.baseurl }}{% link docs/alignment/visualization.md %}) under the Read Mapping tab)
  - `bsstrand` Investigate bisulfite conversion strand label (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) under the Read Mapping tab)
  - `bsconv` Investigate bisulfite conversion rate (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) under the Read Mapping tab)
  - `cinread` Print cytosine-read pair in a long form (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) under the Read Mapping tab)

### Methylation and SNP Extraction

  - `pileup` Generate standard-compliant VCF (see 
  [Read Pileup]({{ site.baseurl }}{% link docs/pileup.md %}))
  - `vcf2bed` Extract mutation or methylation from VCF (see
  [Extracting Methylation and Mutation Information]({{ site.baseurl }}{% link docs/methylextraction.md %}))
  - `mergecg` Merge neighboring C and G in CpG context (see
  [Extracting Methylation and Mutation Information]({{ site.baseurl }}{% link docs/methylextraction.md %}))
  
### Epi-read & Epi-allele

  - `epiread` Convert BAM to epibed format (see
  [Epireads and the epiBED Format]({{ site.baseurl }}{% link docs/epiread/epiread_format.md %}))
  - `rectangle` Convert epiread format to rectangle format (see
  [Epireads and the epiBED Format]({{ site.baseurl }}{% link docs/epiread/old_epiread_format.md %}))
  - `asm` Test allele-specific methylation. (see
  [Allele-specific Methylation]({{ site.baseurl }}{% link docs/allele_meth.md %}))

### Other

  - `version` Print `biscuit` and library versions
  - `qc` Generate QC files from BAM (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}))

## About the project

This package is made by the folks from [Van Andel Institute](https://www.vai.org) with help from prior code base from
the internet.

## Acknowledgement

 - lib/aln was adapted from Heng Li's BWA-mem code.
 - lib/htslib was submoduled from the htslib library.
 - lib/klib was submoduled from Heng Li's klib.
 - This work is supported by NIH/NCI R37CA230748.

## Reference

In preparation
