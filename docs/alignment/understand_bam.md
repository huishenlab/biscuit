---
title: Understand the Output BAM
parent: Read Mapping
nav_order: 2
---

# Understanding the Output BISCUIT BAM

BISCUIT follows the SAM/BAM specification for the SAM/BAM output from `biscuit align`. Definitions for optional tags can
be found below, along with how BISCUIT defines the TLEN column (which differs from how BWA defines it).

## Useful Tags to Know

  - `NM` Number of mismatches, excluding cytosine conversions.
  - `MD` Location of mismatches, following samtools conventions.
  - `ZC` Number of cytosine conversions.
  - `ZR` Number of cytosine retentions.
  - `AS` Best alignment score.
  - `XS` Suboptimal alignment score. This is usually equal to or less than `AS`.  In rare cases, pairing could cause
  `XS` to be greater than `AS`.
  - `RG` Read group.
  - `SA` Other parts of a chimeric primary mapping.
  - `PA` Ratio of best score to alternate score (`AS/XS`). The higher the ratio, the more accurate the position.
  - `XL` Read length excluding adapter.
  - `XA` Location of suboptimal alignments.
  - `XB` Integer pair. The first integer indicates the number of suboptimal mappings in the primary/non-decoy
  chromosomes. The second integer in the pair indicates the number of suboptimal mappings in the ALT/decoy chromosomes.
  For example, `10,5` means ten suboptimal alignments exist on primary/non-decoy chromosomes and five exist on ALT/decoy
  chromosomes.
  - `XR` Reference/chromosome annotation.
  - `YD` Bisulfite conversion strand label.
    - `f` for forward/Watson strand (C&#8594;T from IGV)
    - `r` for reverse/Crick strand (G&#8594;A from IGV)
  - `ZN` Cytosine retention and conversion. Not included by default, but can be added by running `biscuit bsconv`. See
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for more details.
  - `YC` Number of C&#8594;T observations. Not included by default, but can be added by running `biscuit bsstrand`. See
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for more details.
  - `YG` Number of G&#8594;A observations. Not included by default, but can be added by running `biscuit bsstrand`. See
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for more details.

## Insert Size

In BISCUIT, the insert size/TLEN column is given in a different way than BWA.

The insert size as defined by BISCUIT:
```
(right-most coordinate of reverse-mate read) - (left-most coordinate of forward-mate read)
```

The insert size as defined by BWA:
```
-(p0 - p1 + offset)
```
where
```
p0     = { left-most  coordinate of read if read on forward strand }
         { right-most coordinate of read if read on reverse strand }

p1     = { left-most  coordinate of mate if mate on forward strand }
         { right-most coordinate of mate if mate on reverse strand }

         { +1 if p0 > p1 }
offset = {  0 if p0 = p1 }
         { -1 if p0 < p1 }
```
