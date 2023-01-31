---
title: Epireads and the epiBED Format
nav_order: 5
has_children: true
permalink: docs/epiread
---

# Epireads and the Epibed Format

Bisulfite sequencing data may be able to inform how methylation from neighboring positions are linked and how
methylation is linked to mutations (provided they can be unambiguously determined). We introduce the **epibed** format,
which combines the utility of the [epiread format](http://smithlabresearch.org/downloads/methpipe-manual.pdf) with the
standard-compliant BED format. The epibed format is a compact data format used to store CpG (and, in the case of
NOMe-seq, GpC) methylation, as well as SNP and indel information from the same read. It can be useful for estimating the
epiallele fraction and clonal structure/cell population, which is similar to
[MethylPurify](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0419-x).

## Generating Epireads in the Epibed Format

The original epiread format (proposed by the methpipe team), only includes information related to CpGs. BISCUIT builds
off this format to put it in a standard-compliant format, while also containing information related to filtered bases,
SNPs, indels, and, in NOMe-seq mode, GpC methylation. To generate an epibed file, run
```bash
# For WGBS data
biscuit epiread [-B snps.bed] /path/to/my_reference.fa my_output.bam | \
    sort -k1,1 -k2,2n > my_epireads.epibed

# For NOMe-seq data
biscuit epiread [-B snps.bed] -N /path/to/my_reference.fa my_output.bam | \
    sort -k1,1 -k2,2n > my_epireads.epibed
```
The `snps.bed` file can be obtained by running `biscuit vcf2bed -t snp my_pileup.vcf.gz`. If no SNP file is given, the
output will not include information related to SNPs. Note, it is suggested that you sort your epibed file before
running `bgzip` and `tabix` to compress and index your epibed file. This is due to how softclipping is handled in the
BAM file versus how it is handled in the epibed file. The difference can cause some reads to not be coordinate sorted,
which is required when bgzipping and tabixing.

For more help with `epiread`, run `biscuit epiread` in the terminal or check out the
[epiread help page]({{ site.baseurl }}{% link docs/subcommands/biscuit_epiread.md %}).

The columns in the BISCUIT epibed format indicate:

  1. Chromosome name
  2. Start position (0-based)
  3. End position
  4. Read name
  5. Read number in paired-end sequencing
  6. Bisulfite strand (bisulfite Watson (+) or bisulfite Crick (-))
  7. CpG [run length encoded](https://en.wikipedia.org/wiki/Run-length_encoding) string
  8. GpC run length encoded string (only included when run in NOMe-seq mode (`-N`))

An example of the epibed format is:
```
chr1    869996    870097    read_123    1    -    F3x2U1x17U1x1A1U1x7U1x16U1x5U1x1U1x7U1x4U1x17U1x7F3
```

The letters in the run length encoded string are:

  - **F:** A filtered base (filters are the same in `pileup` and `epiread` for consistency)
  - **D:** Base in reference was deleted in read at that location
  - **P:** Soft clipped base
  - **M:** Methylated CpG
  - **U:** Unmethylated CpG
  - **O:** Methylated GpC (open accessibility, only used in NOMe-seq mode)
  - **S:** Unmethylated GpC (shut/closed accessibility, only used in NOMe-seq mode)
  - **A/T/C/G:** SNP base seen in read relative to reference (NOTE, these are *upper* case)
  - **a/t/c/g:** Inserted base included in read, but no in reference (NOTE, these are *lower* case)
  - **x:** Ignored bases (i.e., bases that do not have any of the above occurring at them)

Some things that should be noted about the run length encoded strings:

  1. The M/U/O/S letters always occur at the cytosine location in the reference, regardless of if it is a CpG or GpC
  (the "G" will always get an "x") or what strand the read derives from.
  2. The start to end region (columns 2 and 3) is based on the read start position (properly accounting for soft
  clipping that occurs at the beginning of the read), the length of the read (including insertions), and the number of
  deletions. This means that the region length will always be equal to the length of the run length encoded string, but
  the actual genomic locations relative the reference may differ slightly from those in the read if insertions are
  present.
  3. The default output for `epiread` is the epibed format, but the old BISCUIT epiread format (`-O`) and the pairwise
  epiread format (`-P`) formats can also be accessed with `biscuit epiread`.
  4. As a general note, for reads deriving from the Watson strand, methylation is determined by looking for C&#8594;T
  substitutions (C = methylated, T = unmethylated), while reads deriving from the Crick strand determines methylation by
  looking for G&#8594;A substitutions (G = methylated, A = unmethylated). When running in NOMe-seq mode, BISCUIT only
  allows CpGs in the HCG (Watson) or CGH (Crick) context. This avoids ambiguity in the GCG (or CGC) context as the C (or
  G) can be methylated either due naturally or due to the GpC methyltransferase used in NOMe-seq. With respect to GpCs,
  accessibility is found in the GCH (Watson) or HGC (Crick) context. Accessibility is determined using the C (C = open,
  T = shut) for reads deriving from the Watson strand or the G (G = open, A = shut) for reads deriving from the Crick
  strand. A useful graphic for how methylation is determined can be found
  [below](#three-base-contexts-for-methylation-determination).
  5. In the case of a C&#8594;T SNP in a CpG context (G&#8594;A on the Crick strand), the SNP will take precedence over
  the methylation state. In other words, the T (or A) will appear in the run length encoded string, not a U.
  6. When performing NOMe-seq, methylation due to off-target activity of the M.CviPI enzyme in endogenous CCG contexts
  has been seen by [Kelly et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3514679/). The off-target methylation is a
  small effect (~5% according to Supplemental Figure S8 in Kelly et al), therefore BISCUIT does not filter out
  methylation occurring in a CCG context ([contrary to Bismark](http://felixkrueger.github.io/Bismark/Docs/)). If all
  CCGs were removed (in addition to the GCGs already removed due to ambiguity in NOMe-seq), the number of genome-wide
  cytosines available for analysis would be [cut in half](https://www.nature.com/articles/s41467-018-03149-4). It should
  be noted that, in cases where analyses are hyper-sensitive to methylation levels, CCGs may need to be filtered
  post-hoc (in hopefully a nuanced manner using) the epibed file and a reference genome.

The epibed format can be easily read into
[biscuiteer](https://www.bioconductor.org/packages/release/bioc/html/biscuiteer.html) using the function `readEpibed`,
which creates a GRanges object which seemlessly integrates with the R/Bioconductor suite of tools. `readEpibed` also
provides the option to collapse mate reads in a read pair into a single fragment, as DNA methylation from both reads can
be traced back to a single DNA molecule and represent a physically phased molecular event.

## Three Base Contexts for Methylation Determination

<object data="../assets/2021_09_09_methylation_contexts.pdf" width="525" height="482" type='application/pdf'></object>
