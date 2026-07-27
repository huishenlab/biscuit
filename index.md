---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults
title: Home
layout: default
description: "BISulfite-seq CUI Toolkit (BISCUIT) is a utility suite for
  aligning and analyzing sodium bisulfite conversion-based DNA
  methylation/modification data. Analysis options include DNA methylation and
  mutation calling, allele specific methylation, and post-alignment filtering."
nav_order: 1
---

# BISCUIT - Understand Sequencing Data with Bisulfite Conversion
{: .fs-9 }

[Get started now](#quick-start){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/huishenlab/biscuit){: .btn .fs-5 .mb-4 .mb-md-0 }

---

BISulfite-seq CUI Toolkit (BISCUIT) is a utility suite for analyzing bulk and single-cell sodium bisulfite- or
enzyme-based DNA methylation/modification data, such as WGBS, capture bisulfite sequencing, RRBS, NOMe-seq, and EM-seq.
It was written to perform read alignment, DNA methylation and mutation calling, and allele specific methylation from
bisulfite or bisulfite-like sequencing data.

BISCUIT was developed by Wanding Zhou while he was a member of the [Shen Lab](https://shenlab.vai.org) at Van Andel
Institute. He now holds a faculty position at University of Pennsylvania and Children's Hospital of Philadelphia.
BISCUIT is currently maintained by Jacob Morrison (who also developed the User's Guide website). Current versions of
BISCUIT are available at [https://github.com/huishenlab/biscuit](https://github.com/huishenlab/biscuit), while legacy
versions are located at [https://github.com/zhou-lab/biscuit](https://github.com/zhou-lab/biscuit).

## Reference

Wanding Zhou, Benjamin K Johnson, Jacob Morrison, *et al.*, BISCUIT: an efficient, standards-compliant tool suite for
simultaneous genetic and epigenetic inference in bulk and single-cell studies, *Nucleic Acids Research*, Volume 52,
Issue 6, 12 April 2024, Page e32, [https://doi.org/10.1093/nar/gkae097](https://doi.org/10.1093/nar/gkae097).

## Quick Start

Instructions for downloading and installing BISCUIT can be found in the [Download and Install](#download-and-install)
section below.

The basic workflow to align and extract methylation information using BISCUIT is:
1. Create an index of the reference genome (only needs to be done once for each reference).
2. Align sequenced reads to the reference.
3. Create a pileup VCF of DNA methylation and genetic information.
4. Extract DNA methylation into BED format.

*Note, it is important to use the same reference throughout to ensure consistency in chromosome order and naming.*

A basic example workflow is:
```bash
# 1. Create index of the reference genome
# Gzipped FASTA references can also be used
biscuit index my_reference.fa

# 2. Align sequenced reads to the reference
# Gzipped FASTQ files can also be used
biscuit align \
    -@ NTHREADS \
    -R "my_rg" \
    /path/to/my_reference.fa \
    read1.fastq \
    read2.fastq | \
dupsifter /path/to/my_reference.fa | \
samtools sort \
    -@ NTHREADS \
    -o my_output.bam \
    -O BAM \
    -

samtools index my_output.bam

# 3. Create a pileup VCF of DNA methylation and genetic information
# Also compresses and indexes the VCF
biscuit pileup \
    -@ NTHREADS \
    -o my_pileup.vcf \
    /path/to/my_reference.fa \
    my_output.bam

bgzip -@ NTHREADS my_pileup.vcf
tabix -p vcf my_pileup.vcf.gz

# 4. Extract DNA methylation into BED format
# Also compresses and indexes the BED
biscuit vcf2bed my_pileup.vcf.gz > my_methylation_data.bed
bgzip my_methylation_data.bed
tabix -p bed my_methylation_data.bed.gz
```
This basic order of commands will produce all the necessary files needed to read data into R using the R/Bioconductor
companion package, [biscuiteer](https://www.bioconductor.org/packages/release/bioc/html/biscuiteer.html).

An overview of all available functionalities can be found in the
[Overview of Functionalities](#overview-of-functionalities) section.

## Download and Install

BISCUIT is available through several avenues:

- [Bioconda](#download-via-bioconda)
- [As a precompiled binary](#download-precompiled-binaries)
- [Source code on GitHub](#download-source-code-and-compile)
- [In a Docker container](#download-the-docker-container)

### Download via Bioconda

BISCUIT is available via [bioconda](https://anaconda.org/channels/bioconda/packages/biscuit/overview) for download via
`conda`, `mamba`, or other similar package managers. To install with `conda`, run:

```bash
conda install bioconda::biscuit
```

This will also install `QC.sh`, `build_biscuit_QC_assets.pl`, and `flip_pbat_strands.sh`.

### Download Precompiled Binaries

The [latest release page on GitHub](https://github.com/huishenlab/biscuit/releases/latest) holds precompiled binaries
and several auxiliary scripts for BISCUIT. Precompiled binaries are available for Linux and macOS.

*Note, versions 1.7.1 and earlier will not build on ARM-based CPUs (e.g., Apple M-series CPUs). As of version 1.8.0,
BISCUIT will build on both Intel (x86_64) and ARM CPUs. Further, the precompiled macOS binary for version 1.8.0 and
following is built on an Apple M-series (arm64) CPU. Users who need an Intel (x86_64) macOS binary may install BISCUIT
with bioconda or the source code. Precompiled Linux binaries continue to be built on an Intel (x86_64) CPU.*

Precompiled binaries can be downloaded from GitHub on the release page or directly from the terminal.

On macOS,
```bash
# Latest release
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep darwin_arm64 | cut -d '"' -f 4) --output biscuit
chmod +x biscuit

# Specify version (version 1.8.0 and newer)
curl \
    --output biscuit \
    -L https://github.com/huishenlab/biscuit/releases/download/vX.Y.Z.YYYYMMDD/biscuit_X_Y_Z_darwin_arm64
chmod +x biscuit

# Specify version (version 1.7.0 and older)
curl \
    --output biscuit \
    -L https://github.com/huishenlab/biscuit/releases/download/vX.Y.Z.YYYYMMDD/biscuit_X_Y_Z_darwin_amd64
chmod +x biscuit
```
where `X.Y.Z` is the version number and `YYYYMMDD` is the date of the release.

On Linux,
```bash
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep linux_amd64 | cut -d '"' -f 4) --output biscuit
chmod +x biscuit
```

Auxiliary scripts to generate QC asset files, run QC on BISCUIT outputs, and flip PBAT strands post-alignment can also
be downloaded from the terminal on both macOS and Linux:
```bash
# QC asset build
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep build_biscuit_QC_assets.pl | cut -d '"' -f 4)

# QC bash script
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep QC.sh | cut -d '"' -f 4)

# Flip PBAT strands script
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep flip_pbat_strands.sh | cut -d '"' -f 4)
```

### Download Source Code and Compile

#### Version 1.4.0 and Newer

BISCUIT uses a CMake-based build system as of version 1.4.0. In addition, to `cmake` (minimum version 3.21), BISCUIT has
several required dependencies:

- `zlib`
- `ncurses`
- `pthread`
- `curl`

Optionally, if `libdeflate` is installed, BISCUIT will build `htslib` against `libdeflate`, which provides speed
improvements over `zlib` in `htslib`.

The source code can be retrieved with either of these two commands:
```bash
# git
git clone git@github.com:huishenlab/biscuit.git
cd biscuit

# curl
curl -OL $(curl -s https://api.github.com/repos/huishenlab/biscuit/releases/latest |
    grep browser_download_url | grep release-source.zip | cut -d '"' -f 4)
unzip release-source.zip
cd biscuit-release
```

Regardless of how the source code is retrieved, building BISCUIT proceeds as follows:
```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../ ../
make install
```
This will create a `bin/` directory in the top level directory of BISCUIT with the following executables:

- `biscuit`
- `QC.sh`
- `build_biscuit_QC_assets.pl`
- `flip_pbat_strands.sh`

You may change the directory to install the BISCUIT executables with:
```bash
# During cmake invocation
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/your/install/location ../
make install

# After building
mkdir build && cd build
cmake ../
make
cmake --install --prefix /path/to/your/install/location
```

If BISCUIT is not installed during the build process, the BISCUIT binary can be found in `build/src/biscuit` (relative
to the top level directory of BISCUIT) and the other executables can be found in the `scripts/` directory.

#### Version 1.3.0 and Earlier

Users who want to build BISCUIT version 1.3.0 or earlier must download the `release-source.zip` file available on the
[GitHub releases page](https://github.com/huishenlab/biscuit/releases). `zlib` and `ncurses` are required dependencies
to build version 1.3.0 and earlier.

```bash
unzip release-source.zip
cd biscuit-release
make
```

The QC, asset creator, and strand-flipping scripts can be found in the `scripts/` directory.

### Download the Docker Container

The Docker container can be downloaded from [GitHub](https://github.com/huishenlab/sv_calling_docker) via:
```bash
git clone git@github.com:huishenlab/sv_calling_docker.git
```

For more information about the docker container, see
[Structural Variant Calling]({{ site.baseurl }}{% link docs/structural_variants.md %}).

## Overview of Functionalities

The various subcommands and functionalities available in BISCUIT can be viewed below or by running `biscuit` in your
terminal. Help for each subcommand can be found by running `biscuit <subcommand>` in your terminal or on the
[BISCUIT Subcommands page]({{ site.baseurl }}{% link docs/subcommands/subcommand_help.md %}).

### Read Mapping

  - `index` Index reference genome (see [Read Mapping]({{ site.baseurl }}{% link docs/alignment/alignment.md %}))
  - `align` Align reads to reference (see [Read Mapping]({{ site.baseurl }}{% link docs/alignment/alignment.md %}))

### BAM Operation

  - `tview` View aligned reads in terminal with bisulfite coloring (see
  [Visualization]({{ site.baseurl }}{% link docs/alignment/visualization.md %}) under the Read Mapping tab)
  - `bsstrand` Investigate bisulfite conversion strand label (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) under the Read Mapping tab)
  - `bsconv` Investigate bisulfite conversion rate (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) under the Read Mapping tab)
  - `postalt` Prefer primary chromosomes over unplaced/alt contigs (see
  [Understand the output BAM]({{ site.baseurl }}{% link docs/alignment/understand_bam.md %}) under the Read Mapping tab)
  - `cinread` Print cytosine-read pair in a long form (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) under the Read Mapping tab)

### Methylation and SNP Extraction

  - `pileup` Generate standard-compliant VCF (see [Read Pileup]({{ site.baseurl }}{% link docs/pileup.md %}))
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
  - `help` Print usage and exit
  - `qc` Generate QC files from BAM (see [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}))
  - `qc_coverage` Generate BISCUIT-specific coverage files for QC (see
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}))
  - `bc` Extract cell barcodes from reads (see [Extract Barcodes]({{ site.baseurl }}{% link docs/barcode_extract.md %}))

## About the project

This package is made by the folks from [Van Andel Institute](https://www.vai.org) with help from prior code base from
the internet.

## Acknowledgement

 - `lib/aln` was adapted from Heng Li's BWA-mem code.
 - This work is supported by NIH/NCI R37CA230748.
