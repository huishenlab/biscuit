---
title: Read Pileup
nav_order: 3
---

# Generate Standard VCF Output

The `biscuit pileup` subcommand allows the user to compute cytosine retention and callable SNP mutations.

BISCUIT tries to minimize false positive methylation calls and can be stringent in terms of mutations. If a read has a
high mapping confidence score and shows a mutation with a high base quality score, `biscuit pileup` starts its SNP
processing and *may* revert the methylation call if the SNP interferes with the determination of cytosine retention or
conversion.

## Generating Pileup for a Single Sample

BISCUIT can create a tabix-indexed VCF file for a single sample by running:
```bash
biscuit pileup -o my_pileup.vcf /path/to/my_reference.fa my_output.bam
bgzip my_pileup.vcf
tabix -p vcf my_pileup.vcf.gz
```

To get additional diagnostic information, include `-v 1` in the `biscuit pileup` call. Additional debug information can
be including by using `-v 6` instead.

For more help on available flags, run `biscuit pileup` in the terminal or visit the
[pileup help page]({{ site.baseurl }}{% link docs/subcommand_help.md %}).

### Helpful Flag Definitions

A few flags that can be found in the VCF and are useful to understand are:

  - **DP:** raw number of reads covering that position
  - **CV:** strand-specific coverage on cytosine
  - **BT:** beta value [No. methylated / (No. unmethylated + No. methylated)]
  - **GT:** genotype relative to normal
  - **GL1:** comma-separated list of likelihoods for each of the three genotypes (0/0, 0/1, 1/1)
  - **GQ:** genotype quality (for the called genotype, GT)
  - **SP:** shows allelic support for SNPs
  - **NS:** number of samples with data
  - **CX:** cytosine context (CG, CHH, etc.)
  - **N5:** 5 base sequence context

The VCF header contains descriptions of all INFO and FORMAT flags contained within the VCF file.

### Example Output with Diagnostic Information

The following example shows output `biscuit pileup` in a low coverage, low quality region, with additional diagnostic
information included.

```
# Example command run in BISCUIT
biscuit pileup -o diagnostic.vcf -v 1 /path/to/my_reference.fa output.bam

# Output line from VCF file
chr1    2361154    .    C    .    5    LowQual    NS=1;CX=CG;N5=ACCGG \
    GT:GL1:GQ:DP:SP:CV:BT    0/0:-1,-2,-6:5:1:Y1:1:0.00 \
    DIAGNOSE;RN=0;CN=1;Bs0=T;Sta0=1;Bq0=F;Str0=+;Pos0=44;Rret0=23
```

When diagnostic information is included, the additional tags, with short descriptions, will be included in the VCF
header. In brief, the DIAGNOSE data in the above example shows:

  - The number of high quality reads showing retention (`RN=0`)
  - The number of high quality reads showing conversion (`CN=1`)
  - Retention/conversion pattern for cytosines on the Watson strand (`Bs0=T`)
  - Retention-mutation status internally determined by BISCUIT (`Sta0=1`) See
  [Codes for Retention Mutation Status](#codes-for-retention-mutation-status) for a description of possible values.
  - Phred-scaled base qualities for bases in question (`Bq0=F`)
  - Whether the read is on the forward (`+`) or reverse (`-`) (`Str0=+`)
  - Position along the read where the bases are (`Pos0=44`)
  - Number of retained (i.e., methylated) cytosines on each read (`Rret0=23`)
  - While not shown, information for reads on the Crick strand can be found in the `Bs1`, `Sta0`, etc. tags

When running with the `-v 1` option, BISCUIT will also print positions with no SNP or cytosine methylation, which allows
for differentiating between "no mutation" and "no coverage."

Note, when `biscuit pileup` is run with diagnostic information included the resulting VCF is not standard-compliant.
Therefore, it is not recommended that `-v 1` is included when performing analyses.

### Codes for Retention Mutation Status

In a strictly methylation context, the possible values are:

  - **0:** retention
  - **1:** conversion
  - **2:** unknown

Otherwise, the possible values are:

  - **0:** mutation to A
  - **1:** mutation to C
  - **2:** mutation to G
  - **3:** mutation to T
  - **4:** mutation to Y
  - **5:** mutation to R
  - **6:** retention
  - **7:** conversion
  - **8:** reference base

## Generating Pileup for Multiple Samples

BISCUIT has the ability to put mutation calls and DNA methylation measurements from multiple samples next to each other
in the output VCF by providing `biscuit pileup` with more than one input BAM.
```bash
biscuit pileup -o my_combined_pileup.vcf /path/to/my_reference.fa \
    my_output_1.bam my_output_2.bam [...]
bgzip my_combined_pileup.vcf
tabix -p vcf my_combined_pileup.vcf.gz
```

### Somatic Mode

BISCUIT can call somatic mutations by providing a tumor and matched normal BAM to `pileup`. To run in *somatic mode*,
run `biscuit pileup` with the `-S` flag:
```bash
biscuit pileup -S -o somatic_mode.vcf /path/to/my_reference.fa \
    -T tumor.bam -I normal.bam
bgzip somatic_mode.vcf
tabix -p vcf somatic_mode.vcf.gz
```
Note, the `-T` and `-I` flags must be used to specify the tumor and matched normal BAMs, respectively, when running in
somatic mode.

When running in somatic mode, two additional flags are added to the INFO string:

  - **SS:** somatic state
  - **SC:** somatic score

The following table shows a brief summary of the somatic states.

| State | Description                       | Exists in Tumor | Exists in Normal |
|:-----:|-----------------------------------|:---------------:|:----------------:|
|   0   | Wild type                         |        No       |        No        |
|   1   | Germline                          |        Yes      |        Yes       |
|   2   | Somatic                           |        No       |        Yes       |
|   3   | LOH                               |        Yes      |        No        |
|   4   | Post-transcriptional modification |        NA       |        NA        |
|   5   | Unknown                           |        NA       |        NA        |

The somatic score is calculated based on the posterior probability that, given the read counts in the tumor and normal
samples, the mutation only exists in the tumor sample.

To control the contamination rate, include the `-x` option when running `biscuit pileup`. A higher contamination rate
means a more conservative approach will be used when making somatic calls (i.e., there will be fewer SS=2 calls).

## Ambiguous Alternative Alleles

At times, it is possible that BISCUIT is unable to determine what the alternative allele should be, either because the
alternative is completely unknown or because there is a thymine in the converted strand. When the latter occurs, BISCUIT
is unable to determine whether the thymine is a true thymine or a converted cytosine. If BISCUIT is unable to
distinguish what the alternative allele should be, it will provide a `N` in the ALT column of the VCF.

It should also be noted that ambiguous alleles can appear in the allelic support entry of the FORMAT string. In this
case, either a `Y` (representing a C or T) or a `R` (representing a G or A) will appear. The `Y` and `R` representations
are based on the [IUPAC codes](https://www.bioinformatics.org/sms/iupac.html). Depending on the alternate allele
provided in the VCF, these ambiguous alleles may or may not contribute to the variant allele frequency that is
calculated for this particular variant.

## Features

  - Fast, multi-way pileup of multiple samples with VCF output
  - Computes beta values
  - Computes genotype, genotype likelihood, and genotype quality
  - Computes somatic score and somatic state when tumor and matched normal samples are provided (`-S` option)
  - Calls ambiguous alternative allele
  - Distinguishes between ambiguous alternative allele and multiple alternative alleles
  - Generates a coordinate-sorted output VCF when provided a sorted input BAM
  - Flexible read filtering based on retention number, mapping quality, duplicate marking, and mate pairing
  - Flexible base filtering using base quality and distance to the read ends
