---
title: Epi-read & Epi-allele
permalink: /epiread_format/
nav_order: 5
---

# Epi-read and Epi-allele
{: .no_toc }

Bisulfite sequencing data may be able to inform how methylation from neighboring
positions are linked and how methylation is linked to mutations (provided they
can be unambiguously determined). The epiread format is a commonly used, compact
data format used to store the CpG retention pattern, as well as SNP information
on the same read. It is useful for estimating the epiallele fraction and clonal
structure/cell population
([Li et al. Genome Biology 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4),
[Zheng et al. Genome Biology 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0419-x)).

## Generating Epireads

The original epiread format (proposed by the methpipe team), only includes
information related to CpGs. BISCUIT extends this format to also include SNP
information. The columns in the BISCUIT epiread format indicate:

  1. Chromosome name
  2. Read name
  3. Read position in paired-end sequencing
  4. Bisulfite strand (bisulfite Watson (+) or bisulfite Crick (-))
  5. Position of the cytosine in the first CpG (0-based)
  6. Retention pattern ("C" for retention or "T" for conversion) for all
  CpGs covered
  7. Position of the first SNP, if a SNP location file is provided
  8. Base call of all SNPs covered

An example of the epiread format is:
```
chr19  NS500653:8:HF5FGBGXX:3:12402:11299:9856   1  +  3040315  CCCCTCCC  .         .
chr19  NS500653:8:HF5FGBGXX:3:12609:17196:5738   1  +  3055491  T         .         .
chr19  NS500653:8:HF5FGBGXX:1:23304:20253:14257  1  +  3078472  CC        3078510   T
chr19  NS500653:8:HF5FGBGXX:4:22509:19067:5776   1  +  3078472  CC        3078510   C
chr19  NS500653:8:HF5FGBGXX:4:22611:9688:19912   2  +  3078946  C         .         .
chr19  NS500653:8:HF5FGBGXX:4:22602:25913:14920  2  +  3078946  CC        3078982   A
```

To produce an epiread formatted file, you need to run the `epiread` subcommand
with a BAM file, a FASTA file of the reference genome, and, optionally, a BED
file for SNPs.
```bash
biscuit epiread [other options] [-B snps.bed] /path/to/my_reference.fa my_output.bam
```

The SNP BED file can be obtained by running `biscuit vcf2bed -t snp
my_pilefup.vcf.gz`.  If no SNP file is supplied, the output does not include the
extra columns related to SNPs. To get back the original epiread format, run `cut
-f 1,5,6` on the output epiread file.

To test all SNP-CpG pairs, include the `-P` flag in your command prompt. Note,
if looking for allele-specific methylation, you will need to run with the `-P`
flag and include a SNP BED file. For more help on available flags, run
`biscuit epiread` in the terminal.

## Paired-end Epireads

DNA methylation information from both mate reads are physically "phased"
molecular events that can be traced back to the same DNA molecule. The read name
and read position in the single-end epiread output can be used to collate mate
reads in a read pair. The default behavior of the `epiread` subcommand focuses
only on the primary mapping. The following `awk` command gives a nice, compact
file for a paired-end epiread format:
```bash
$ sort -k2,2 -k3,3n single_end.epiread |
$     awk 'BEGIN{ qname="" ; rec="" }
$          qname == $2 { print rec"\t"$5"\t"$6"\t"$7"\t"$8 ; qname="" }
$          qname != $2 { qname=$2 ; rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 ; pair=$3}'
```

In the following output, the columns represent:

  1. Chromosome name
  2. Bisulfite strand (bisulfite Watson (+) or bisulfite Crick (-))
  3. Location of the cytosine in the first CpG covered (0-based) in read 1
  4. Retention pattern of all CpGs covered in read 1
  5. Location of the first SNP covered in read 1
  6. Base call of all SNPs covered in read 1
  7. Location of the cytosine in the first CpG covered (0-based) in read 2
  8. Retention pattern of all CpGs covered in read 2
  9. Location of the first SNP covered in read 2
  10. Base call of all SNPs covered in read 2

```
chr19  -  3083513  CCCCCCC        3083495  ATT  3083513  CCCCCCC      3083495  ATT
chr19  -  3083545  CCTCCCCCCT     .        .    3083527  CCCCTCCCCCC  3083523  AG
chr19  +  3083616  TTTTTTTTTTT    .        .    3083722  TTTTTTTTTTT  .        .
chr19  -  3083630  CCCCCTCCCCCCT  .        .    3083616  TCCCCCTCCCC  .        .
chr19  +  3083638  TTTTTTTTTTTT   .        .    3083705  TTTTTTTTTTT  .        .
```

## NOMe-seq Epireads

Including the `-N` option in the `epiread` subcommand allows the GCH and HCG
retention states in NOMe-seq data to be listed side by side. For example,
```bash
$ biscuit vcf2bed -t snp my_pileup.vcf.gz > snp.bed
$ biscuit epiread -B snp.bed -N -q 20 -n 3 -N /path/to/my_reference.fa my_output.bam |
    gzip -c > single_end.epiread.gz
$ # Collating paired epireads
$ zcat single_end.epiread.gz |
$ sort -k2,2 -k3,3n |
$ awk 'BEGIN{ qname="" ; rec="" }
$      qname == $2 { print rec"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 ; qname="" }
$      qname != $2 { qname=$2 ; rec=$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 ; pair=$3}' |
$ sort -k1,1 -k3,3n | gzip -c > paired.epiread.gz
```

The paired output looks like,
```
chr18  +  10689  CCC     10663  TTCT         10696  T  10689  CCC    10663  TTCT      10696  T
chr18  +  10689  CCC     10694  TT           10696  T  10689  CCC    10694  TT        10696  T
chr18  -  10689  CCC     10699  TTT          10696  T  10689  CCC    10699  TTTT      10696  T
chr18  +  10703  CCCCCC  10694  CCCCCTCTTTT  10696  T  10753  CCCCC  10743  CTCTTTTT  .      .
chr18  -  11134  CCCC    11131  CTTT         .      .  11134  CCCC   11131  CTTT      .      .
chr18  -  11344  CC      11339  CCTTT        .      .  11344  CC     11339  CCTT      .      .
```

The output columns are:
  1. Chromosome name
  2. Bisulfite strand (bisulfite Watson (+) or bisulfite Crick (-))
  3. Location of the cytosine in the first CpG covered (0-based) in read 1
  4. Retention pattern of all CpGs covered in read 1
  5. Location of the cytosine in the first GpCpH covered (0-based) in read 1
  6. Retention pattern of all GpCpHs covered in read 1
  7. Location of the first SNP covered in read 1
  8. Base call of all SNPs covered in read 1
  9. Location of the cytosine in the first CpG covered (0-based) in read 2
  10. Retention pattern of all CpGs covered in read 2
  11. Location of the cytosine in the first GpCpH covered (0-based) in read 2
  12. Retention pattern of all GpCpHs covered in read 2
  13. Location of the first SNP covered in read 2
  14. Base call of all SNPs covered in read 2

To be precise, BISCUIT only counts HpCpG for CpG. However, BISCUIT records the
location of C, regardless of whether it is the G or the C that is measured. In
other words, if there is a TGCGA and a read was bisulfite-treated on the Crick
strand, BISCUIT will only record the fourth G and not the third C, but will still
use the position of the third C as the position for the CpG. This is done to be
consistent with standard BS-seq.

## Generate Rectangular Forms

By combining the `epiread` and `rectangle` subcommands, you can generate a matrix
with columns representing CpGs and rows representing reads. This is important for
understanding [Methylation Haplotype Load](https://dx.doi.org/10.1038%2Fng.3805),
epi-polymorphism, methylation entropy, and allele-specific methylation.

```bash
$ biscuit rectangle /path/to/my_reference.fa my_epiread_file.epiread
```

Note, `my_epiread_file.epiread` can only have a single chromosome in it for
`rectangle` to work. For more details on the available `rectangle` flags, run
`biscuit rectangle` in the terminal.
