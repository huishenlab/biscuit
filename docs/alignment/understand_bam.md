---
title: Understand the Output BAM
parent: Read Mapping
nav_order: 1
---

# Understanding BISCUIT BAMs

The BISCUIT BAM is slightly different from typical BAM files with regards to
tags and how insert sizes are computed.

### Useful Tags to Know

  - `NM` Number of mismatches, excluding cytosine conversions.
  - `XA` Location of suboptimal alignments.
  - `XB` Integer pair. The first integer indicates the number of suboptimal
  mappings in the primary/non-decoy chromosomes. The second integer in the pair
  indicates the number of suboptimal mappings in the ALT/decoy chromosomes. For
  example, `10,5` means ten suboptimal alignments exist on primary/non-decoy
  chromosomes and five exist on ALT/decoy chromosomes.
  - `ZC` Number of cytosine conversions.
  - `ZR` Number of cytosine retentions.
  - `AS` Best alignment score.
  - `XS` Suboptimal alignment score.  This is usually equal to or less than `AS`.
  In rare cases, pairing would cause a `XS` greater than `AS`.
  - `MD` Location of mismatches, following samtools conventions.
  - `PA` Ratio of best score to alternate score (`AS/XS`). The higher the ratio,
  the more accurate the position.
  - `SA` Other parts of a chimeric primary mapping.
  - `YD` Bisulfite conversion strand label.
    - `f` for forward/Watson strand (C<p>&#8594;</p>T from IGV)
    - `r` for reverse/Crick strand (G<p>&#8594;</p>A from IGV)

See [BAM Operation]({{ site.baseurl }}{% link docs/alignment/bam_operation.md %})
for how to add a `ZN` tag for cytosine conversion under `CpG` and other non-CpG
dinucleotide sequence contexts.

### Insert Size

In BISCUIT, the insert size/TLEN is printed in a different way than BWA. Here,
the TLEN is the actual insert size:

```
(right-most coordinate of reverse-mate read) - (left-most coordinate of forward-mate read)
```
