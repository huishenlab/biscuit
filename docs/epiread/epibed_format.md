---
title: The epiBED Format
parent: Epireads and the epiBED Format
nav_order: 1
permalink: /epibed_format/
---

# The epiBED Format

The epiBED format has gone through two iterations, versions 1.0 and 2.0. Version 1.0 is the output from BISCUIT versions
1.0-1.1, while version 2.0 is output from BISCUIT version 1.2 and following.

The epiBED format makes use of [run length encoded](https://en.wikipedia.org/wiki/Run-length_encoding) (RLE) strings
to store relevant information out of read sequences in an aligned BAM. An assortment of letters are used to specify
information:

  - **F:** A filtered base (filters are the same in `pileup` and `epiread` for consistency)
  - **D/d:** Base in reference was deleted in read at that location
  - **P:** Soft clipped base
  - **M:** Methylated CpG
  - **U:** Unmethylated CpG
  - **O:** Methylated GpC (open accessibility, used in NOMe-seq mode)
  - **S:** Unmethylated GpC (shut/closed accessibility, used in NOMe-seq mode)
  - **A/T/C/G/R/Y:** SNP base seen in read relative to reference (NOTE, these are *upper* case); R and Y represent A/G
  and C/T, respectively
  - **a/t/c/g/i:** Inserted base included in read, but no in reference (NOTE, these are *lower* case); **i** is used as
  a placeholder in version 2.0 CpG and GpC RLE strings
  - **x:** Ignored bases (i.e., bases that do not have any of the above occurring at them)

The epiBED format can be easily read into
[biscuiteer](https://www.bioconductor.org/packages/release/bioc/html/biscuiteer.html) using the function `readEpibed`,
which creates a GRanges object which seemlessly integrates with the R/Bioconductor suite of tools. `readEpibed` also
provides the option to collapse mate reads in a read pair into a single fragment, as DNA methylation from both reads can
be traced back to a single DNA molecule and represent a physically phased molecular event.

## epiBED version 2.0

The columns are:

  1. Chromosome name
  2. Start position (0-based)
  3. End position (1-based, non-inclusive)
  4. Read name
  5. Read number in paired-end sequencing
  6. Bisulfite strand (OT/CTOT (+) or OB/CTOB (-))
  7. CpG RLE string
  8. GpC RLE string (`.` when `-N` not included in command)
  9. Variant RLE string (SNPs and indels)

An example of the format is:
```
chr11    2132661    2132762    read_123    1    -    F3x22Fx72F3    .    F3x3RxRx16Fx11Rx60F3
```

Notes:

  1. The M/U/O/S now occur at the position they are found in the BAM. In other words, the methylation status deriving
  from the C in a CpG or GpC will be placed at the C's position, whereas the methylation status from a G would be placed
  at the G's position.
  2. The start position is based on the on the read start position listed in the BAM adjusted by any soft clipped bases
  at the start of the read (i.e., start - n_softclipped).
  3. The end position is set such that `end - start` would be the length of the reference genome the read spans. In
  other words, `end = start + read_length + n_deletions - n_insertions`, where `n_deletions` is the number of deletions
  in the read and `n_insertions` is the number of insertions. Relative to the run length encoded string,
  `end - start = rle_length + n_insertions`.
  4. As a general note, for reads deriving from the OT/CTOT strands, methylation is determined by looking for C&#8594;T
  substitutions (C = methylated, T = unmethylated), while reads deriving from the OB/CTOB strand determines methylation
  by looking for G&#8594;A substitutions (G = methylated, A = unmethylated). When running in NOMe-seq mode, BISCUIT only
  allows CpGs in the HCG (OT/CTOT) or CGH (OB/CTOB) context. This avoids ambiguity in the GCG (or CGC) context as the C
  (or G) can be methylated either naturally or due to the GpC methyltransferase used in NOMe-seq. With respect to GpCs,
  accessibility is found in the GCH (OT/CTOT) or HGC (OB/CTOB) context. Accessibility is determined using the C
  (C = open, T = shut) for reads deriving from the OT/CTOT strands or the G (G = open, A = shut) for reads deriving from
  the OB/CTOB strands. A useful graphic for how methylation is determined can be found below.
  5. When performing NOMe-seq, methylation due to off-target activity of the M.CviPI enzyme in endogenous CCG contexts
  has been seen by [Kelly et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3514679/). The off-target methylation is a
  small effect (~5% according to Supplemental Figure S8 in Kelly et al), therefore BISCUIT does not filter out
  methylation occurring in a CCG context ([contrary to Bismark](http://felixkrueger.github.io/Bismark/Docs/)). If all
  CCGs were removed (in addition to the GCGs already removed due to ambiguity in NOMe-seq), the number of genome-wide
  cytosines available for analysis would be [cut in half](https://www.nature.com/articles/s41467-018-03149-4). It should
  be noted that, in cases where analyses are hyper-sensitive to methylation levels, CCGs may need to be filtered
  post-hoc (in hopefully a nuanced manner) using the epiBED file and a reference genome.

<object data="../assets/2023_01_31_methylation_contexts.pdf" width="525" height="482" type='application/pdf'></object>

Differences from version 1.0:

  1. Methylation and SNPs match the output from `biscuit pileup` in version 2.0. Note, this requires using an unfiltered
  SNP BED file from BISCUIT.
  2. There are always nine columns in version 2.0 versus the 7 or 8 in version 1.0.
  3. Version 2.0 uses "i" and "d" to mark insertions and deletions as well as not including SNPs in the CpG and GpC RLE
  strings. All variants are collected in the variant RLE string instead.
  4. Methylation status is printed at the location it is found in `pileup` in version 2.0 (see above). Version 1.0 would
  always put the methylation status at the location of the C.
  5. It is possible in version 2.0 for methylation and a SNP to be called in the same location. This is possible when
  the alternate base is not a T (A) or when the alternate base is a T (A) and the allele frequency is < 0.05.

## epiBED version 1.0

The columns are:

  1. Chromosome name
  2. Start position (0-based)
  3. End position
  4. Read name
  5. Read number in paired-end sequencing
  6. Bisulfite strand (bisulfite Watson (+) or bisulfite Crick (-))
  7. CpG run length encoded string
  8. GpC run length encoded string (only included when run in NOMe-seq mode (`-N`))

An example of the epiBED format is:
```
chr1    869996    870097    read_123    1    -    F3x2U1x17U1x1A1U1x7U1x16U1x5U1x1U1x7U1x4U1x17U1x7F3
```

Notes:

  1. The M/U/O/S letters always occur at the cytosine location in the reference, regardless of if it is a CpG or GpC
  (the "G" will always get an "x") or what strand the read derives from.
  2. The start to end region (columns 2 and 3) is defined in the same way as Version 2.0.
  3. Methylation and accessibility is found in the same manner in both version 1.0 and 2.0.
  4. In the case of a C&#8594;T SNP in a CpG context (G&#8594;A on the Crick strand), the SNP will take precedence over
  the methylation state. In other words, the T (or A) will appear in the run length encoded string, not a U.

