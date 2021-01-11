---
title: Allele-specific Methylation
nav_order: 6
---

# Allele-specific Methylation

Allele-specific methylation can be extracted from the epiread files via:
```bash
$ biscuit epiread -P -o my_output.epiread -B snps.bed \
      /path/to/my_reference.fa my_output.bam |
$ sort -k1,1 -k2,2n -k3,3n my_output.epiread > my_output.sorted.epiread
$ biscuit asm my_output.sorted.epiread > my_output.asm
```

An example of the `.asm` file format is:
```
chr1    17360   17406   T/C T/C 4   2   0   1   0.428571    0.459426
chr1    17380   17378   A/G T/C 1   9   0   1   1.000000    0.946485
```

The columns represent:

  1. Chromosome name
  2. SNP location
  3. CpG location
  4. SNP_1/SNP_2 (SNPs for row_1/row_2 in contingency table)
  5. CpG_1/CpG_2 (CpGs for col_1/col_2 in contingency table)
  6. SNP_1-CpG_1 count (contingency table value in row_1-col_1)
  7. SNP_1-CpG_2 count (contingency table value in row_1-col_2)
  8. SNP_2-CpG_1 count (contingency table value in row_2-col_1)
  9. SNP_2-CpG_2 count (contingency table value in row_2-col_2)
  10. Two-tail Fisher's exact p-value
  11. $\chi^2$ p-value
