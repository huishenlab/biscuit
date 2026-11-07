# BISCUIT Testing Framework

## Initial Framework

The initial framework for testing BISCUIT compares the output from the version of BISCUIT built in the top level
directory of your BISCUIT repository (`../` relative to the `test/` directory you're reading this in) against a fixed
set of outputs from a previous version of BISCUIT. It assumes these previous outputs are correct; therefore, if a bug is
being fixed that affects outputs, it will be expected that one or more tests will fail. In this instance, you will have
to update the file being compared against to make the test pass before continuing. In time, unit tests will hopefully be
added that will test BISCUIT functions against known, expected results.

### Test Data

Tests for BISCUIT are run with chromosome 22. The reference FASTA (`chr22.fa.gz`) was downloaded from UCSC, decompressed
for use in Sherman (see below), and the recompressed using BGZIP compression to enable a FASTA index to be created.

Simulated reads used to run the tests were created with [Sherman](https://github.com/FelixKrueger/Sherman). 1000 paired
end reads with 150 bp length were created and used throughout.

Illumina 5-base reads (`5base_pe.bam`, used to test `pileup --chemistry 5base`) are a 10 kb slice of the public NA12878
100 ng demo sample (`Demo-5base-gDNA-Sample1-NA12878-100ng-B-B01`) from Illumina's "Illumina 5-Base DNA" BaseSpace
project, aligned by DRAGEN 4.4.6 to GRCh38. Its chr22 matches the UCSC chr22 used here outside the masked
pericentromere. Only the chr22 `@SQ` line is kept, `@PG` lines are dropped, as are reads whose mate maps to another
contig. DRAGEN's `XG`/`XR`/`XM` tags are removed so the test covers taking the strand from read orientation:

```
samtools view --no-PG -h -T GRCh38_full_analysis_set_plus_decoy_hla.fa.gz NA12878-100ng-5base.cram \
        chr22:20020000-20030000 \
    | awk -F'\t' -v OFS='\t' '
        /^@SQ/ { if ($2 == "SN:chr22") print $1, $2, $3; next }
        /^@PG/ { next }
        /^@/   { print; next }
        $7 == "=" || $7 == "*" || $7 == "chr22" { gsub(/\t(XG|XR|XM):Z:[^\t]*/, ""); print }' \
    | samtools sort --no-PG -o 5base_pe.bam --write-index -
```

To download and extract files for performing tests, run `python setup_tests.py`. For those interested in knowing how
these files were created, the code is recreated below.

```
# Directory for storing data files
mkdir -p data/ref
cd data/ref

# Download and decompression
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz
cd ../

# Create simulated data
Sherman \
    -l 150 \
    -n 1000 \
    --genome_folder ref/ \
    -pe

# Compression and FASTA index creation
bgzip chr22.fa
samtools faidx chr22.fa.gz

# Setup data directories
mkdir dynamic && cd dynamic
mkdir 00_index 01_align 02_pileup 03_vcf2bed 04_mergecg 05_bsconv 06_bsstrand 07_cinread

# Output files for each subcommand were created following the command line calls given in ../../dynamic/run_*.py
# To avoid having to update these if the commands change, please see the commands in those files for specifics
```

### Dynamic Tests

To run the tests found in the `dynamic` directory, run

```
cd dynamic

python dynamic_tests.py
```

Dynamic tests are designed to test the output of BISCUIT against an expected result. In this case, we expect continuity
of the outputs across BISCUIT versions. Therefore, we will compare updates to BISCUIT against previous versions. The
dynamic tests are run in the `dynamic` directory and are compared against results found in `data/dynamic`. Files for
testing to be aware of are described below.

- `config.toml`: Controls which subcommands to test and whether to force recreating the files. Options set to `true`
will be run (or outputs will be recreated regardless of if they exist or not), while options set to `false` will not be
run. In cases where `run.* = true` and `force.* = false`, the files will only be created the first time, but comparisons
against expected results will always be run.
- `dynamic_tests.py`: Runs the tests specified in `config.toml`.
- `run_*.py`: Tests for specific subcommands of BISCUIT.
