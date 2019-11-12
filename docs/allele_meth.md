---
title: Allele-specific Methylation
nav_order: 7
---

# Allele-specific Methylation

Allele-specific methylation can be extracted from the epiread files via:
```bash
$ sort -k1,1 -k2,2n -k3,3n my_epiread_file.epiread > my_output.epiread
$ biscuit asm my_output.epiread > my_output.asm
```
