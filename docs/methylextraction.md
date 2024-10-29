---
title: Extracting Methylation and Mutation Information
nav_order: 5
---

# Extract Methylation and Mutation Information

## Make BED files

After generating the VCF file with both genetic and methylation information, beta values and coverage can be extracted
to study the methylation levels at any sequenced CpGs. The following shows how to extract this information from the VCF
file:
```bash
$ biscuit vcf2bed -t cg my_pileup.vcf.gz > my_pileup.bed
```
For BISCUIT version 1.5.0 and earlier, the minimum coverage by default is 3. This setting was good for traditional WGBS;
however, for low input and single cell protocols, this option should be changed to 1 (`biscuit vcf2bed -k 1 ...`).
BISCUIT version 1.6.0 and later have the default value switched to 1 to avoid user confusion.

The `-t` flag can be used to retrieve mutation and other information, including

  - `snp` for SNP information
  - `c` for all cytosines
  - `hcg` for HCG (from NOMe-seq)
  - `gch` for GCH (from NOMe-seq)

The columns for `-t snp` are:

  1. Chromosome
  2. Start position (0-based)
  3. End position
  4. Reference base
  5. Alternate base
  6. Genotype (GT tag)
  7. Base support (SP tag)
  8. Coverage (AC tag)
  9. Allele frequency (AF1 tag)

Columns 6-9 are repeated for each sample found in the VCF file.

The columns for the other `-t` options (i.e., all methylation related options) are:

  1. Chromosome
  2. Start position (0-based)
  3. End position
  4. Methylation fraction
  5. Coverage

Columns 4 and 5 are repeated for each sample in the VCF file. If the `-e` option is included, then four columns are
inserted between columns 3 and 4 listed above:

  - Reference base
  - Cytosine group (CG, CHG, CHH)
  - 2-base context (CA, CC, CG, CT)
  - 5-base context
  
Information about the sequence context for each row can be included by adding the `-e` flag to the command, while
filtering out low coverage rows can be done using the `-k` flag.

For more help on available flags, run `biscuit vcf2bed` in the terminal or visit the
[vcf2bed help page]({{ site.baseurl }}{% link docs/subcommands/biscuit_vcf2bed.md %}).

## Merge Neighboring C and G in CpG Context

At times, it might be desirable to merge the methylation status across the C and G in a CpG dinucleotide context.
BISCUIT provides an easy-to-use subcommand, `mergecg`, to do this:
```bash
$ biscuit mergecg /path/to/my_reference.fa my_pileup.bed
```

For more help on available flags, run `biscuit mergecg` in the terminal or visit the
[mergecg help page]({{ site.baseurl }}{% link docs/subcommands/biscuit_mergecg.md %}).

## Towards the Bismark COV format

Starting in version 1.3.0, BISCUIT can output a format that can be easily worked into the Bismark COV file format for
easy use in downstream tools that expect such a file. To produce this file, run
```
biscuit vcf2bed -c my_pileup.vcf.gz > my_beta_m_u.bed
```
which will produce a BED-compliant file with these columns:

  1. Chromosome
  2. Start position (0-based)
  3. End position
  4. Methylation percentage
  5. M (number of methylated reads covering locus)
  6. U (number of unmethylated reads covering locus)

The Bismark COV file format is not BED-compliant, so to create a "true" COV file, you can use AWK in conjunction with
`vcf2bed`:
```
# Pipe straight from vcf2bed
biscuit vcf2bed -c my_pileup.vcf.gz | \
awk -v OFS='\t' '{ print $1, $2+1, $3, $4, $5, $6 }' > my_beta_m_u.cov

# Create from already processed file
biscuit vcf2bed -c my_pileup.vcf.gz > my_beta_m_u.bed
awk -v OFS='\t' '{ print $1, $2+1, $3, $4, $5, $6 }' my_beta_m_u.bed > my_beta_m_u.cov
```

This file will keep methylation across the strands separate, as is normally done in `biscuit vcf2bed`. To merge
methylation across strands, run:
```
biscuit vcf2bed my_pileup.vcf.gz | \
biscuit mergecg -c /path/to/my_reference.fa - | \
awk -v OFS='\t' '{ print $1, $2+1, $3-1, $4, $5, $6 }' > my_merged_beta_m_u.cov
```

Note, the `-c` flag is *only* included in the `biscuit mergecg` call. If you would prefer a BED-compliant file, remove
the AWK command.
