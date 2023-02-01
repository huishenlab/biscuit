---
title: Snakemake Pipeline
nav_order: 8
---

# Snakemake Pipeline

A Snakemake pipeline is [available on GitHub](https://github.com/huishenlab/Biscuit_Snakemake_Workflow) for streamlining
analyses.

## Dependencies

  - Required:
    - `snakemake` (version 6.0+)
    - `biscuit` (version 1.2.0+)
    - `samtools` (version 1.12+)
    - `htslib` (version 1.12+)
    - `samblaster` (version 0.1.26+)
    - bedtools
    - pigz
    - GNU parallel
    - FastQC
    - MultiQC
    - Python 3.7+ with pandas, numpy, matplotlib, and seaborn
  - Optional:
    - TrimGalore! (required for trimming adaptors)
    - R with tidyverse, ggplot, patchwork, and viridis (required for plotting methylation controls)
    - Bismark (required when running fastq_screen)
    - fastq_screen (required when running fastq_screen)
    - preseq (required for finding library complexity, version 3.1.2+ must be compiled with htslib enabled)

Two things of note, 1) it is easiest when working with `snakemake` to install `mamba` using `conda` when running with
`--use-conda`, and 2) it is preferable to install `snakemake` using `conda`, rather than using a module. This is due to
potential conflicts between packages (such as `matplotlib`) that can be found in the python distribution associated with
the snakemake module.

## Components of Workflow

The following components are generally in order, but may run in a different order, depending on exact dependencies
needed.
  - [default off] Generate asset files used during QC related rules
  - [default off] Modify and index genome reference to including methylation controls
  - [default off] Trim adapters and/or hard clip R2
  - [default off] Run Fastq Screen in bisulfite mode
  - Run FastQC on raw FASTQ files
  - Alignment, duplicate tagging, indexing, flagstat of input data (biscuitBlaster v1 and v2)
  - Methylation information extraction (BED Format)
  - Merge C and G beta values in CpG dinucleotide context
  - [default off] SNP and epiBED extraction
  - [default off] Run Preseq on aligned BAM
  - MultiQC with BICUIT QC modules specifically for methyaltion data
  - [default off] Generate plots of the observed / expected coverage ratio for different genomic features
  - [default off] Generate percentage of covered CpGs and CpG island coverage figures
  - [default off] QC methylated and unmethylated controls
  - [default off] Find binned average methylation
  - [default off] Find binned methylation centered on provided regions

Many options can be easily specified in the `config.yaml`! Otherwise, the commands in the Snakefile can also be modified
to meet different needs.

## Running the Workflow

  - Clone the repo `git clone git@github.com:huishenlab/Biscuit_Snakemake_Workflow.git`
  - Place *gzipped* FASTQ files into `raw_data/` directory
  - Replace the example `bin/samples.tsv` with your own `config/samples.tsv` sample sheet containing:
    - One row for each sample
    - Three columns per row (separated by a tab - "\t"):
      - `sample_name` name of sample to be used throughout processing
      - `fq1` name of R1 file for `sample_name` in `raw_data/`, multiple FASTQs can be specified (comma-separated)
      - `fq2` name of R2 file for `sample_name` in `raw_data/`, multiple FASTQs can be specified (comma-separated)
    - Any additional columns are ignored
  - Modify `config/config.yaml` to specify your:
    - Reference genome
    - BISCUIT index
    - BISCUIT QC assets (see [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for details)
    - If you are using environmental modules on your system, you can set the locations in the corresponding location. By
      default, the pipeline will use `conda`/`mamba` to download the required packages. Note, if using the modules and a
      module is not available, snakemake gives a warning but will run successfully *as long as the required executables
      are in the path*.
    - Toggle optional workflow components (change from False to True)
    - Specify other run parameters
  - Submit the workflow to an HPC using command similar to `bin/run_snakemake_workflow.sh`
    - `bin/run_snakemake_workflow.sh` works for submitting to a PBS/Torque queue system
      - Submit via: `qsub -q [queue_name] bin/run_snakemake_workflow.sh`
      - Make sure the queue you submit to is able to submit jobs from jobs running on that queue
      - If the nodes are not able to submit jobs, the snakemake pipeline will not be able to run properly
    - `bin/run_snakemake_workflow.sh` can be easily modified for submission on other queue systems
  - Snakemake can also be run on the command line:
    - `snakemake --use-envmodules --cores 1`

## After Workflow Completion

  - Analysis-related output can be found in the directory specified by `config["output_directory"]` (`analysis/` by
  default)
    - `BISCUITqc/` output from `QC.sh`
    - `align/` output from biscuitBlaster pipeline
    - `epiread/` epiBED files from `biscuit epiread` that can be used as input to `biscuiteer::readEpibed()` (included
    if `epiread: True` in `config.yaml`)
    - `fastq_screen/` reports from fastq_screen
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
    - `snakemake --dag | dot -Tpng > my_dag.png`
