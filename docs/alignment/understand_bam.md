---
title: Understand the Output BAM
parent: Read Mapping
nav_order: 2
---

# Understanding the Output BISCUIT BAM

BISCUIT follows the [SAM/BAM](http://samtools.github.io/hts-specs/SAMv1.pdf) specification for the SAM/BAM output from
`biscuit align`. Definitions for SAM auxiliary tags can be found below, along with a discussion on how the insert size
is defined in BISCUIT, as this was corrected in version 1.10.0.

## Useful Tags to Know

  - `NM` Number of non-cytosine-conversion mismatches. Note, this does not match exactly with the SAM spec for the NM
  tag. To recreate the exact NM tag as defined by the spec, add the values in the NM and ZC tags.
  - `MD` Location of mismatches, following the SAM specification.
  - `ZC` Number of conversions (C&rarr;T for OT/CTOT reads or G&rarr;A for OB/CTOB).
  - `ZR` Number of retentions (retained C's for OT/CTOT reads or G's for OB/CTOB).
  - `AS` Best alignment score.
  - `XS` Suboptimal alignment score. This is usually equal to or less than `AS`.  In rare cases, pairing could cause
  `XS` to be greater than `AS`.
  - `RG` Read group.
  - `SA` Other parts of a chimeric primary mapping.
  - `PA` Ratio of best score to alternate score (`AS/XS`). The higher the ratio, the more accurate the position.
  - `XL` Read length excluding adapter.
  - `XA` Location of suboptimal alignments.
  - `XB` Integer pair. The first integer indicates the number of suboptimal mappings in the primary/non-alternate
  chromosomes. The second integer in the pair indicates the number of suboptimal mappings in the ALT/decoy chromosomes.
  For example, `10,5` means ten suboptimal alignments exist on primary/non-decoy chromosomes and five exist on ALT/decoy
  chromosomes.
  - `XR` Reference/chromosome annotation.
  - `YD` Bisulfite conversion strand label.
    - `f` for OT/CTOT strands (C&rarr;T from IGV)
    - `r` for OB/CTOB strands (G&rarr;A from IGV)
  - `MC` CIGAR string for mate/next segment
  - `MQ` Mapping quality for mate/next segment
  - `ZN` Cytosine retention and conversion. Not included by default, but can be added by running `biscuit bsconv`. See
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for more details.
  - `YC` Number of C&rarr;T observations. Not included by default, but can be added by running `biscuit bsstrand`. See
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for more details.
  - `YG` Number of G&rarr;A observations. Not included by default, but can be added by running `biscuit bsstrand`. See
  [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for more details.
  - `CB` Extracted cell barcode. See [Barcode Extraction]({{ site.baseurl }}{% link docs/barcode_extract.md %}) for
  more details.
  - `RX` Extracted unique molecular index (UMI). See
  [Barcode Extraction]({{ site.baseurl }}{% link docs/barcode_extract.md %}) for more details.

## Insert Size

With the release of version 1.10.0, BISCUIT now follows the BWA convention for defining the insert size (the TLEN column
in the SAM specification).

The insert size as defined by BWA (and now BISCUIT) is:
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

Prior to version 1.10.0, the insert size/TLEN column for BISCUIT was given in a different way than BWA. In these
versions, the insert size was defined as:
```
(right-most coordinate of reverse-mate read) - (left-most coordinate of forward-mate read)
```

From a practical standpoint, the magnitude of the insert size is the same between the two methods. However, BWA follows
the SAM spec which defines the insert size as `TLEN = (right-most coordinate) - (left-most coordinate)`. The right-most
read is given an insert size of `-TLEN`, while the left-most read is given an insert size of `+TLEN`. The signed TLEN
value will cut the insert size distribution in half for `QC.sh` and `biscuit qc` as the quality control will no longer
be double counting the insert sizes for paired-end reads.
