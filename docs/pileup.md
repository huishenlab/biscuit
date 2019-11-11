---
title: Read Pileup
nav_order: 3
---

# Generate Standard VCF Output

## The pileup subcommand

Like samtools, BISCUIT can extract DNA methylation as well as genetic information.
The following shows how to produce a tabix-indexed VCF file:
```bash
$ biscuit pileup /path/to/my_reference.fa my_output.bam -o my_pileup.vcf
$ bgzip my_pileup.vcf
$ tabix -p vcf my_pileup.vcf.gz
```
