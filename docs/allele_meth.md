---
title: Allele-specific Methylation
nav_order: 6
---

# Allele-specific Methylation

Allele-specific methylation can be extracted from the pairwise epiread files via:
```bash
biscuit epiread -@ NTHREADS -P -B snps.bed \
    /path/to/my_reference.fa my_output.bam | \
sort -k1,1 -k2,2n -k3,3n my_output.pairwise > my_output.sorted.pairwise

biscuit asm my_output.sorted.pairwise > my_output.asm
```

An example of the `.asm` file format is:
```
chr11    2132802    2132844    T/C    T/C    0    1    2    0    3.333333e-01    8.326452e-02
chr11    2132851    2132804    T/C    T/C    1    0    0    1    1.000000e+00    1.572992e-01
```

The columns represent:

  1. Chromosome name
  2. SNP location
  3. CpG location
  4. SNP<sub>1</sub>/SNP<sub>2</sub> (SNPs for row<sub>1</sub>/row<sub>2</sub> in contingency table)
  5. CpG<sub>1</sub>/CpG<sub>2</sub> (CpGs for col<sub>1</sub>/col<sub>2</sub> in contingency table)
  6. SNP<sub>1</sub>-CpG<sub>1</sub> count (contingency table value in row<sub>1</sub>-col<sub>1</sub>)
  7. SNP<sub>1</sub>-CpG<sub>2</sub> count (contingency table value in row<sub>1</sub>-col<sub>2</sub>)
  8. SNP<sub>2</sub>-CpG<sub>1</sub> count (contingency table value in row<sub>2</sub>-col<sub>1</sub>)
  9. SNP<sub>2</sub>-CpG<sub>2</sub> count (contingency table value in row<sub>2</sub>-col<sub>2</sub>)
  10. Two-tail Fisher's exact p-value
  11. &chi;<sup>2</sup> p-value
