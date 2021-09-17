---
title: Extracting Methylation and Mutation Information
nav_order: 4
---

# Extract Methylation and Mutation Information

## Make BED files

After generating the VCF file with both genetic and methylation information, beta values and coverage can be extracted
to study the methylation levels at any sequenced CpGs. The following shows how to extract this information from the VCF
file:
```bash
$ biscuit vcf2bed -t cg my_pileup.vcf.gz > my_pileup.bed
```

The `-t` flag can be used to retrieve mutation and other information, including

  - `snp` for SNP information
  - `c` for all cytosines
  - `hcg` for HCG (from NOMe-seq)
  - `gch` for GCH (from NOMe-seq)
  
Information about the sequence context for each row can be included by adding the `-e` flag to the command, while
filtering out low coverage rows can be done using the `-k` flag.

For more help on available flags, run `biscuit vcf2bed` in the terminal or visit the
[vcf2bed help page]({{ site.baseurl }}{% link docs/subcommand_help.md %}).

## Merge Neighboring C and G in CpG Context

At times, it might be desired that the methylation status be merged across the C and G in a CpG dinucleotide context.
BISCUIT provides an easy-to-use subcommand, `mergecg`, to do this:
```bash
$ biscuit mergecg /path/to/my_reference.fa my_pileup.bed
```

For more help on available flags, run `biscuit mergecg` in the terminal or visit the
[mergecg help page]({{ site.baseurl }}{% link docs/subcommand_help.md %}).
