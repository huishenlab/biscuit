---
title: Glossary of Terms
nav_order: 8
---

# Glossary of Terms
{: .no_toc }

## Table of Content
{: .no_toc .text-delta }

1. TOC
{:toc}

This glossary serves as a reference for terms that can be easily misunderstood,
are used multiple ways in bioinformatics, or just need more explaining than is
given in earlier sections. Care has been taken to try and use the definitions
provided here in the same way throughout the BISCUIT documentation.

## General Terms

  - **Coverage:** Fraction of the genome that has been sequenced to a given
  *depth*.
  - **Depth:** Average number of reads covering each each base in the genome.
  - **Watson/Forward/Top Strand:** DNA strand where reading from the 5' to 3'
  end proceeds from left to right. The choice of the Watson and Crick strands
  is arbitrary, but a consistent definition should be used across an entire
  chromosome.
  - **Crick/Reverse/Bottom Strand:** DNA strand where reading from the 5' to 3'
  end proceeds from right to left. The choice of the Watson and Crick strands
  is arbitrary, but a consistent definition should be used across an entire
  chromosome.
  - **CpG (CG) Dinucleotide:** A cytosine-guanine dinucleotide (often referred
  to as just CpG). Most DNA methylation occurs in at CpG sites.
  - **CpG Terminology:**
    - *CpG Island:* Region of about 200 bp (or longer) in length with a GC
    percentage greater than 50% and has an observed to expected CpG ratio
    greater than 60%.
    - *CpG Shore:* Region up to 2 kb upstream/downstream from a CpG island.
    - *CpG Shelf:* Region between 2 and 4 kb upstream/downstream from a CpG
    island.
    - *CpG Open Sea:* Isolated CpGs in the genome.

## Bisulfite Sequencing Terms

  - **Whole Genome Bisulfite Sequencing (WGBS):** Method of whole genome DNA
  sequencing that provides additional information regarding the DNA methylation
  status of single cytosine bases. In this method, DNA is treated with sodium
  bisulfite prior to sequencing, which converts unmethylated cytosines into
  uracils and leaves the methylated cytosines unconverted. Following
  sequencing, the unmethylated cytosines (which were converted to uracils),
  appear as thymines.
  - **Bisulfite Watson (BSW) / Original Top (OT):** Reads originating from the
  bisulfite converted Watson strand.
  - **Bisulfite Crick (BSC) / Original Bottom (OB):** Reads originating from
  the bisulfite converted Crick strand.
  - **Bisulfite Watson Reverse (BSWR) / Complement to Original Top (CTOT):**
  Complement to the bisulfite converted Watson strand.
  - **Bisulfite Crick Reverse (BSCR) / Complement to Original Bottom (CTOB):**
  Complement to the bisulfite converted Crick strand.
  - **Under (or Incomplete) Conversion:** Occurs when unmethylated cytosines are
  not converted during the sodium bisulfite treatment, which can lead to an
  overestimate of the methylation levels in the sample.
  - **Over Conversion:** Occurs when the sodium bisulfite treatment begins
  converting methylated cytosines, which can lead to an underestimate of the
  methylation levels in the sample.
  - **Directional Library:** Library preparation format where only only reads
  from the original top and original bottom strands can be sequenced.
  - **Non-directional Library:** Library preparation format where reads can be
  sequenced from any of the following: the original top strand, the original
  bottom strand, the complement to the original top strand, or the complement
  to the original bottom strand. Reads from each strand will be mapped with
  approximately the same frequency.
  - **Post-Bisulfite Adapter Tagging (PBAT):** In general, bisulfite treatment
  can be harsh on DNA. PBAT is a DNA library preparation method for WGBS where
  bisulfite treatment precedes adapter tagging and two rounds of random primer
  extension. This method is less harsh on DNA and allows for an increased
  number of unamplified reads from lower quantities of DNA than other methods.


## Methylation Analysis Terms

  - **Conversion Rate:** Fraction of thymines that BISCUIT recognizes as
  converted cytosines relative to the total number of cytosines (unconverted and
  converted). Refers to the number of unmethylated cytosines in the genome.
  - **Retention Rate:** Fraction of cytosines that BISCUIT recognizes as
  unconverted cytosines relative to the total number of cytosines (unconverted
  and converted). Refers to the number of methylated cytosines in the genome.
  - **M-Bias:** In theory, the methylation of cytosines within a read should be
  independent of position. Therefore, a plot of the methylation fraction as a
  function of the read position should be a flat, horizontal line. In practice,
  systematic sequencing and base-calling errors can lead to biases in the
  methylation level at either the 5' or 3' ends of the read. This bias is
  referred to as the *M-bias*.
  - **Epi-Allele:** Similar to genetic alleles, epigenetic marks (like DNA
  methylation) can be heritable at specific regions in the genome of some
  species. To distinguish the heritable epigenetic marks from their genetic
  counterparts, the term *Epi-allele (of epiallele)* is used.
  - **Allele Specific Methylation:** Methylation that occurs on one allele in a
  diploid cell, but not the other.
  - **Read Level Retention:** In BISCUIT, is calculated by directly averaging
  the retention rate across reads, without grouping the reads by their mapped
  locations.
  - **Base Pair Retention:** In BISCUIT, is calculated by first averaging the
  retention rate over each base/CpG and then averaging over all bases.
