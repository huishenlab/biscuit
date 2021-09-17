---
title: Snakemake Pipeline
nav_order: 8
---

# Snakemake Pipeline

## Dependencies

  - Required:
    - BISCUIT version 1.0.0 or greater
    - Samtools
    - Snakemake
    - Samblaster
    - htslib
    - bedtools
    - pigz
    - GNU parallel
    - FastQC
    - MultiQC
  - Optional:
    - R with tidyverse, patchwork, and viridis (required for plotting methylation controls)
    - Bismark (required when running fastq_screen)
    - fastq_screen (required when running fastq_screen)

## Components of Workflow

In (approximate) running order:

  - Rename FASTQ files based on sample name in sample sheet
  - [default off] Modify and index genome reference to include methylation controls
  - Trim adapters
  - Run fastq_screen on FASTQ files
  - Align, duplicate mark, sort, index, and run `samtools flagstat` on BAM
  - Create pileup VCF
  - Extract methylation and merge CpGs
  - Create BISCUIT QC files
  - Extract SNPs from pileup VCF
  - Create epibed files
  - Run MultiQC over QC files
  - [default off] Perform methylation control QC and create figure

## Running the Workflow

  1. Clone the repo `git clone git@github.com:huishenlab/WGBS_Biscuit_Snakemake.git`
  2. Place *gzipped* FASTQ files into `raw_data/` directory
  3. Replace the example `bin/samples.tsv` with your own `bin/samples.tsv` sample sheet containing:
    - A row for each sample
    - Three columns per row (separated by a tab - "\t"):
      - `sample_name` name of sample to be used throughout processing
      - `fq1` name of R1 file for `sample_name` in `raw_data/`
      - `fq2` name of R2 file for `sample_name` in `raw_data/`
    - Any additional columns are ignored
  4. Modify `bin/config.yaml` to specify your
    - Reference genome
    - BISCUIT index
    - BISCUIT QC assets (see [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for details)
    - Environmental modules
      - If not using the modules, replace the module path with "null/path" or a similar string
      - It not using the modules, the required executable *must* be in PATH
      - If modules are not available (whether "null/path" or an incorrect path), snakemake gives a warning but will run
      successfully if the required executable is in PATH
    - Toggle optional workflow components
    - Specify other run parameters
  5. Run the first rule of Snakemake on the command line: `snakemake --cores 2 --use-envmodules --until get_R1_R2_files`
    - **This rule needs to be run separately first for the correct R1 and R2 files to be passed to `biscuit_align`**
    - Collects the list of comma separated R1 and R2 files in `bin/samples.tsv` and renames them
      - Only a few seconds per sample and allows quick debugging of missing input files
  6. Submit the workflow to an HPC using command similar to `bin/run_snakemake_workflow.sh`
    - `bin/run_snakemake_workflow.sh` works for submitting to a PBS/Torque queue system
      - Submit via: `qsub -q [queue_name] bin/run_snakemake_workflow.sh`
      - Make sure the queue you submit to is able to submit jobs from the nodes available to that queue
      - If the nodes are not able to submit jobs, the snakemake pipeline will not be able to run properly
    - `bin/run_snakemake_workflow.sh` can be easily modified for submission on other queue systems
  7. Snakemake can also be run on the command line: `snakemake --use-envmodules --cores 1`
    - When running on the command line the `--use-envmodules` is required

## After Workflow Completion

  - Analysis-related output can be found in the `analysis/` directory
    - `BISCUITqc/` output from `QC.sh`
    - `align/` output from biscuitBlaster pipeline
    - `epiread/` epibed files from `biscuit epiread` that can be used as input to `biscuiteer::readEpibed()` (included
    if `epiread: 1` in `config.yaml`)
    - `fastq_screen` reports from fastq_screen
    - `multiqc/` MultiQC output with BISCUIT, fastq_screen (if run), and trim_galore (if run) reports
    - `pileup/` VCF and merged CpG BED files that can be used as inputs to `biscuiteer::readBiscuit()`
    - `qc_vectors/` methylation control BED files and beta value/coverage figure
    - `snps/` SNP BED files (included if generate_snps: 1 in `config.yaml`)
    - `trim_reads` trimmed FASTQ files and FastQC reports
  - Log files can be found in the `logs/` directory
  - Benchmarking files can be found in the `benchmarks/` directory

## Example Dataset

The cloned Snakemake repository comes with a five sample test dataset to see how this workflow works on your system. To
run the test dataset, copy the ten `.fq.gz` files in `bin/working_example_dataset` into `raw_data/` and use the default
`bin/samples.tsv` file. This set of files should be mapped to the human genome.

## Useful Commands

  - To perform a dry run of the commands that will be run by snakemake
    - `snakemake -npr`
  - To unlock the pipeline after a manually aborted run
    - `snakemake --unlock --cores 1`
  - To create a workflow diagram of your run
    - snakemake --dag | dot -Tpng > my_dag.png`
