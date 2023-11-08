---
title: Barcode Extraction
nav_order: 2
---

# Extracting Cell Barcodes from scWGBS Protocols

Cell barcodes are commonly used in scRNA-seq protocols to overcome DNA input issues and increase throughput by merging
cDNA from separate wells into a single well.  In a similar manner, the single-cell WGBS (scWGBS) protocol, snmc-seq2,
uses the same methodology for increasing throughput. While there aren't many scWGBS protocols that use cell barcoding,
it seems likely barcodes (and unique molecular indexes (UMIs) down the road) are used more frequently in protocols.

BISCUIT provides the `biscuit bc` subcommand as a way to remove uncorrected cell barcodes from FASTQ files, which can
either be written to new FASTQ files or fed directly into `biscuit align`. Further, the output from `biscuit bc` matches
the output from the commonly used tool, [umi-tools](https://github.com/CGATOxford/UMI-tools), which allows for the user
to perform barcode correction with `umi-tools` and then use the resulting FASTQs as input to `biscuit align`.

## Barcode Extraction with BISCUIT

As an overview, `biscuit bc` takes FASTQs as input and outputs the processed reads to standard out:
```
# Single-end data
biscuit bc read1.fq.gz

# Paired-end data
biscuit bc read1.fq.gz read2.fq.gz
```
For paired-end data, the output is printed in an interleaved form (i.e., read 1 followed by read 2 for each read pair).
**It is important to include both FASTQs as input for paired-end data, as the barcode and an artificial UMI (needed for
compatibility with `umi-tools`) are written to the read name of both read 1 and read 2.**

You can set the geometry of where the cell barcode is with the `-m/--mate`, `-s/--bc-start`, and `-l/--bc-length`
options.
```
biscuit bc -m 1 -s 1 -l 8 read1.fq.gz read2.fq.gz
```
`-m/--mate` specifies the mate the barcode is in (1 or 2), `-s/--bc-start` is the start position of the barcode (given
as a 1-based number), and `l/--bc-length` is the length of the barcode.

Rather than writing the processed reads to standard out, you can instead write them to a FASTQ file with
`-o/--output-prefix`, which is the basename of your output FASTQ files. For single-end data, `.fq.gz` will be appended
to your output prefix. Paired-end data will have `_R1.fq.gz` or `_R2.fq.gz` appended to the prefix for reads 1 and 2,
respectively.
```
biscuit bc -o my_processed_reads read1.fq.gz read2.fq.gz
```

For more help on available flags, run `biscuit bc` in the terminal or visit the
[bc help page]({{ site.baseurl }}{% link docs/subcommands/biscuit_bc.md %}).

### Passing Ouput to `biscuit align`

You can pass the output of `biscuit bc` straight to `biscuit align`:
```
# Single-end data
biscuit bc read1.fq.gz | \
biscuit align -9 -@ NTHREADS /path/to/my_reference.fa -

# Paired-end data
biscuit bc read1.fq.gz read2.fq.gz | \
biscuit align -p -9 -@ NTHREADS /path/to/my_reference.fa -
```
The `-9` flag tells `biscuit align` to extract barcodes from the read name, while the `-p` flag for paired-end data is
needed to handle the interleaved input from `biscuit bc`.

You can also align FASTQs that were created by `biscuit bc`:
```
# Single-end data
biscuit bc -o reads my_se_data.fq.gz
biscuit align -9 -@ NTHREADS /path/to/my_reference.fa reads.fq.gz

# Paired-end data
biscuit bc -o reads my_pe_data_R1.fq.gz my_pe_data_R2.fq.gz
biscuit align -9 -@ NTHREADS /path/to/my_reference.fa read_R1.fq.gz read_R2.fq.gz
```

In either case, the cell barcode is placed in the `CB` SAM tag ("the optionally corrected cellular barcode sequence"),
while the UMI is placed in the `RX` SAM tag (the "[s]equence bases from the unique molecular identifier ... [where] the
value may be non-unique in the file.")

## Barcode Correction with `umi-tools`

Barcode correction and extraction can be performed using `umi-tools` with the help of
[synthbar](https://github.com/jamorrison/synthbar) (here shown for paired-end sequencing).

### Step 1: Add synthetic UMI with synthbar

`umi-tools` expects both a cell barcode and a UMI. At this time, there are no scWGBS protocols that include UMIs;
therefore, a synthetic barcode must be added to the read with the cell barcode in it. `synthbar` is a tool that can add
synthetic barcodes to cells, but for our purposes can also be used to add a synthetic UMI:
```
synthbar -b AAAAAAAA barcoded_reads_R1.fastq.gz | \
gzip > barcoded_with_umi_R1.fastq.gz
```
This will prepend `AAAAAAAA` to the start of each read in the input FASTQ file, which can be treated as the UMI for each
read.

### Step 2: Create Whitelist

Using `umi-tools` requires to steps. The first step is to process your FASTQ with the cell barcodes in it and try to
find the most likely true cell barcode:
```
umi_tools whitelist \
    --log2stderr \
    --bc-pattern=NNNNNNNNCCCCCCC \
    --stdin barcoded_with_umi_R1.fastq.gz \
> whitelist.txt
```
Here, `NNNNNNNNCCCCCCC` specifies a an 8 basepair UMI followed by a 7 basepair cell barcode.

### Step 3: Correct and Extract Barcodes

The results from `whitelist.txt` are then fed into the code that extracts the UMI and cell barcode into the readname:
```
umi_tools extract \
    --bc-pattern=NNNNNNNNCCCCCCC \
    --whitelist=whitelist.txt \
    --error-correct-cell \
    --stdin barcoded_with_umi_R1.fastq.gz \
    --stdout barcoded_with_umi_R1_corrected.fastq.gz \
    --read2-in barcoded_reads_R2.fastq.gz \
    --read2-out barcoded_reads_R2_corrected.fastq.gz
```

### Passing to `biscuit align`

The resulting FASTQs can then be passed into `biscuit align` in a similar manner to extracting barcodes with
`biscuit bc`:
```
biscuit align -9 -@ NTHREADS /path/to/my_reference.fa \
    barcoded_with_umi_R1_corrected.fastq.gz \
    barcoded_reads_R2_corrected.fastq.gz
```
