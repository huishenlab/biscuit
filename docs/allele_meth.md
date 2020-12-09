---
title: Allele-specific Methylation
nav_order: 6
---

# Allele-specific Methylation

Allele-specific methylation can be extracted from the epiread files via:
```bash
$ sort -k1,1 -k2,2 -k3,3n my_epiread_file.epiread > my_output.epiread
$ biscuit asm my_output.epiread > my_output.asm
```
