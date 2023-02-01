---
title: Glossary of Terms
nav_order: 10
---

# Glossary of Terms

This glossary serves as a reference for terms that can be easily misunderstood, are used multiple ways in
bioinformatics, or just need more explaining than is given in earlier sections. Care has been taken to try and use the
definitions provided here in the same way throughout the BISCUIT documentation.

## General Terms

  - **Coverage:** Average number of reads covering each base in the genome.
  - **Depth:** Number of reads covering a given base in the genome.
  - **Top/Forward/Watson Strand:** DNA strand where reading from the 5' to 3' end proceeds from left to right. The
  choice of the Top and Bottom strands is arbitrary, but a consistent definition should be used across an entire
  chromosome.
  - **Bottom/Reverse/Crick Strand:** DNA strand where reading from the 5' to 3' end proceeds from right to left. The
  choice of the Top and Bottom strands is arbitrary, but a consistent definition should be used across an entire
  chromosome.
  - **CpG (CG) Dinucleotide:** A cytosine-guanine dinucleotide (often referred to as just CpG). Most DNA methylation
  occurs in at CpG sites.
  - **CpH (CH) Dinucleotide:** A dinucleotide of a cytosine followed by a adenine, thymine, or another cytosine (i.e.,
  not a guanine). Most of these are unmethylated, although there is some biological instances of CpA methylation.
  - **CpG Terminology:**
    - *CpG Island:* Region of about 200 bp (or longer) in length with a GC percentage greater than 50% and has an
    observed to expected CpG ratio greater than 60%.
    - *CpG Shore:* Region up to 2 kb upstream/downstream from a CpG island.
    - *CpG Shelf:* Region between 2 and 4 kb upstream/downstream from a CpG island.
    - *CpG Open Sea:* Isolated CpGs in the genome between CpG shelves

## Bisulfite Sequencing Terms

  - **Whole Genome Bisulfite Sequencing (WGBS):** Method of whole genome DNA sequencing that provides additional
  information regarding the DNA methylation status of single cytosine bases. In this method, DNA is treated with sodium
  bisulfite prior to sequencing, which converts unmethylated cytosines into uracils and leaves the methylated cytosines
  unconverted. Following sequencing, the unmethylated cytosines (which were converted to uracils), appear as thymines.
  - **Whole Genome (DNA) Methylation Sequencing (WGMS):** With the invention of EM-seq, sodium bisulfite is no longer
  the only method for targeting DNA methylation using whole genome sequencing. WGMS serves as an umbrella term for
  describing both WGBS and EM-seq.
  - **Original Top (OT) / Bisulfite Watson (BSW):** Reads originating from the original top strand.
  - **Original Bottom (OB) / Bisulfite Crick (BSC):** Reads originating from the original bottom strand.
  - **Complement to Original Top (CTOT) / Bisulfite Watson Reverse (BSWR):** Complement to the original top strand.
  - **Complement to Original Bottom (CTOB) / Bisulfite Crick Reverse (BSCR):** Complement to the original bottom strand.
  - **Under (or Incomplete) Conversion:** Occurs when unmethylated cytosines are not converted during the sodium
  bisulfite or enzyme treatment, which can lead to an overestimate of the methylation levels in the sample.
  - **Over Conversion:** Occurs when the sodium bisulfite or enzyme treatment begins converting methylated cytosines,
  which can lead to an underestimate of the methylation levels in the sample.
  - **Directional Library:** Library preparation format where only read 1 derives from the OT (or OB) strands and read 2
  derives from the CTOT (or CTOB) strands.
  - **Non-directional Library:** Library preparation format where read 1 derives from any of the four possible strands
  (OT/OB/CTOT/CTOB) and read 2 derives from the corresponding complement strand (CTOT/CTOB/OT/OB, respectively). Reads
  from each strand will be mapped with approximately equal frequency.
  - **Post-Bisulfite Adapter Tagging (PBAT):** In general, bisulfite treatment can be harsh on DNA. PBAT is a DNA
  library preparation method for WGBS where bisulfite treatment precedes adapter tagging and two rounds of random primer
  extension. This method is less harsh on DNA and allows for an increased number of unamplified reads from lower
  quantities of DNA than other methods.

## Methylation Analysis Terms

  - **Conversion Rate:** Fraction of thymines that BISCUIT recognizes as converted cytosines relative to the total
  number of cytosines (unconverted and converted). Refers to the number of unmethylated cytosines in the genome.
  - **Retention Rate:** Fraction of cytosines that BISCUIT recognizes as unconverted cytosines relative to the total
  number of cytosines (unconverted and converted). Refers to the number of methylated cytosines in the genome.
  - **M-Bias:** In theory, the methylation of cytosines within a read should be independent of position. Therefore, a
  plot of the methylation fraction as a function of the read position should be a flat, horizontal line. In practice,
  systematic sequencing and base-calling errors can lead to biases in the methylation level at either the 5' or 3' ends
  of the read. This bias is referred to as the *M-bias*.
  - **Epi-Allele:** Similar to genetic alleles, epigenetic marks (like DNA methylation) can be heritable at specific
  regions in the genome of some species. To distinguish the heritable epigenetic marks from their genetic counterparts,
  the term *Epi-allele (or epiallele)* is used.
  - **Allele Specific Methylation:** Methylation that occurs on one allele in a diploid cell, but not the other. A
  common place this occurs is in imprinted regions.
  - **Read Level Retention:** In BISCUIT, read level retention is calculated by averaging the retention rate across
  reads, without grouping the reads by their mapped locations.
  - **Base Pair Retention:** In BISCUIT, base level retention is calculated by first averaging the retention rate over
  each base/CpG and then averaging over all bases.
