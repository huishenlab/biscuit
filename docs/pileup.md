---
title: Read Pileup
nav_order: 3
---

# Generate Standard VCF Output

The `biscuit pileup` subcommand allows the user to computes cytosine retention
and callable SNP mutations.

BISCUIT tries to minimize false positive methylation calls and can be stringent
in terms of mutations. If a read has a high mapping confidence score and shows a
mutation with a high base quality score, `biscuit pileup` starts its SNP
processing and *may* revert the methylation call if the SNP interferes with the
determination of cytosine retention or conversion.

## Generating Pileup for a Single Sample

Like samtools, BISCUIT can extract DNA methylation as well as genetic information.
The following shows how to produce a tabix-indexed VCF file:
```bash
$ biscuit pileup /path/to/my_reference.fa my_output.bam -o my_pileup.vcf
$ bgzip my_pileup.vcf
$ tabix -p vcf my_pileup.vcf.gz
```

A few flags that can be found in the VCF and are useful to understand are:

  - **DP:** the total read coverage
  - **CV:** cytosine strand coverage
  - **BT:** Beta value [No. methylated / (No. unmethylated + No. methylated)]
  - **GT:** Genotype
  - **GP:** Comma-separated list of likelihoods for each of the three genotypes
  - **GQ:** Genotype quality (for the called genotype)
  - **SP:** Shows allelic support for SNPs
  - **NS:** Number of samples with data
  - **CX:** Cytosine context (CG, CHH, etc.)
  - **N5:** 5 base sequence environment

To get additional diagnostic information, add `-v 1` to the `biscuit pileup`
call. If you would prefer to receive additional debut information, use `-v 6`
instead.

### Example Output with Diagnostic Information

The following example shows output `biscuit pileup` in a low coverage, low
quality region, including additional diagnostic information.

```
chr20   47419734    .   C   .   19  PASS    NS=1;CX=CG;N5=CACGG DP:GT:GP:GQ:SP:CV:BT    12:0/0:1,6,19:19:C7:2:0.00
```
```
DIAGNOSE;
RN=0;CN=2;
Bs0=CCCTT; Sta0=66677; Bq0=;7);<; Str0=++---; Pos0=74,68,51,44,23; Rret0=21,19,21,22,18;
Bs1=CCCCCCC; Sta1=8888888; Bq1=;;):<<2; Str1=+++----; Pos1=69,63,63,62,59,20,1; Rret1=19,19,15,17,19,26,25
```

When `biscuit pileup` is run with diagnostic information printed out, the
additional tags will be included in the VCF header, with short descriptions
included. In brief, the diagnostic information shows the the number of high
quality reads that show conversion (`CN=2`) and the number of high quality
reads that show retention (`RN=0`). Information about the bisulfite Watson
strand (BSW) can be found on the line starting with *Bs0*, while information
about the bisulfite Crick (BSC) strand can be found on the line starting with
*Bs1*. The BSW strand has five reads, with base identities of CCCTT (see `Bs0`
flag). On the other hand, the BSC strand has seven reads, with base identities
of CCCCCCC (`Bs1` flag). The `Sta0` and `Sta1` flags show the retention-mutation
status code that is internally determined by BISCUIT (see *Codes for
retention-mutation status* below for a description on each possibly value). The
`Bq0` and `Bq1` flags show the Phred-scaled base qualities for the specific
bases in question. The `Str0` and `Str1` flags show whether the reads are on the
forward (`+`) or reverse (`-`) strand. The `Pos0` and `Pos1` flags show where
the bases are positioned in the read, while the `Rret0` and `Rret1` flags show
the number of retention (i.e. methylated cytosines) on each read.

With regard to the five bases on the BSW strand, two reads showing retention
have a low base quality, while the base in the other read is positioned in the
second to last base in read, so it is not counted as a retained read (this
filtering is a tunable option using the `-e` flag in `biscuit pileup`). The
two bases which suggest conversion both have high qualities, so these are
counted towards the number of converted reads.

Note, when running with the `-v 1` option, BISCUIT will also print positions
with no SNP or cytosine methylation, which allows for differentiation between
"no mutation" and "no coverage."

#### Codes for Retention-Mutation Status

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

BISCUIT has the ability to put mutation calling and DNA methylation measurements
from multiple samples next to each other by providing `biscuit pileup` with more
than one input BAM.
```bash
$ biscuit pileup /path/to/my_reference.fa my_output_1.bam my_output_2.bam [...] -o my_combined_pileup.vcf
$ bgzip my_combined_pileup.vcf
$ tabix -p vcf my_combined_pileup.vcf.gz
```

### Somatic Mode

If provided BAMs from a tumor and matched normal, it is possible to call somatic
mutations with BISCUIT. To run in "somatic mode," run `biscuit pileup` with the
`-T` flag:
```bash
$ biscuit pileup /path/to/my_reference.fa tumor.bam normal.bam -T -o somatic_mode.vcf
$ bgzip somatic_mode.vcf
$ tabix -p vcf somatic_mode.vcf.gz
```
Note, the tumor BAM must be given before the normal BAM, as BISCUIT treats the
first BAM as the tumor sample.

When running in somatic mode, two additional flags are added to the INFO string:

  - **SS:** somatic state
  - **SC:** somatic score

The somatic score is calculated based on the posterior probability that, given
the read counts in the tumor and normal samples, the mutation only exists in the
tumor sample.

The following table shows a brief summary of the somatic states.

| State | Description                       | Exists in Tumor | Exists in Normal |
|:-----:|-----------------------------------|:---------------:|:----------------:|
|   0   | Wild type                         |        No       |        No        |
|   1   | Germline                          |        Yes      |        Yes       |
|   2   | Somatic                           |        No       |        Yes       |
|   3   | LOH                               |        Yes      |        No        |
|   4   | Post-transcriptional modification |        NA       |        NA        |
|   5   | Unknown                           |        NA       |        NA        |

To control the contamination rate, include the `-x` option when running `biscuit
pileup`. A higher contamination rate means a more conservative approach will be
used when making somatic calls (i.e. there will be fewer SS=2 calls).

## Ambiguous Alternative Alleles

At times, it is possible that BISCUIT is unable to determine what the
alternative allele should be, either because the alternative is completely
unknown or because there is a thymine in the converted strand. When the latter
occurs, BISCUIT is unable to determine whether the thymine is a true thymine or
a converted cytosine. If BISCUIT is unable to distinguish what the alternative
allele should be, it will provide a `N` in the ALT column of the VCF.

It should also be noted that ambiguous alleles can appear in the allelic support
entry of the FORMAT string. In this case, either a `Y` (representing a C or T)
or a `R` (representing a G or A) will appear. The `Y` and `R` representations
are based on the [IUPAC codes](https://www.bioinformatics.org/sms/iupac.html).
Depending on the alternate allele provided in the VCF, these ambiguous alleles
may or may not contribute to the variant allele frequency that is calculated for
this particular variant.

## Features

  - Fast, multi-way pileup of multiple samples with VCF output
  - Computes beta values
  - Computes genotype, genotype likelihood, and genotype quality
  - Computes somatic score and somatic status when tumor and matched normal
  samples are provided (-T option)
  - Calls ambiguous alternative allele
  - Distinguishs between ambiguous alternative allele and multiple alternative
  alleles
  - Generates a coordinate-sorted output VCF when provided a sorted input BAM
  - Flexible read filtering based on retention number, mapping quality,
  duplicate marking, and mate pairing
  - Flexible base filtering using base quality and distance to the read ends
