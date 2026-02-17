---
title: Snakemake Pipeline
nav_order: 9
---

# Snakemake Pipeline

A Snakemake pipeline is [available on GitHub](https://github.com/huishenlab/Biscuit_Snakemake_Workflow) for streamlining
analyses.

## Dependencies

The following dependencies are downloaded when running with `--use-conda`, otherwise you must have these in your PATH.
  - `snakemake` (version 8.0+, needed before running pipeline)
  - `snakemake-executor-plugin-cluster-generic` (version 1.0.9, needed before running pipeline)
  - `mamba` (version 2.1.0, needed before running pipeline)
  - `biscuit` (version 1.8.0)
  - `htslib` (version 1.22.1)
  - `samtools` (version 1.22.1)
  - `dupsifter` (version 1.3.0)
  - `parallel` (version 20230322)
  - `bedtools` (version 2.31.1)
  - `perl` (version 5.32.1, only required if running `beta_bigwigs`)
  - `ucsc-bedgraphtobigwig` (version 455, only required if running `beta_bigwigs`)
  - `preseq` (version 3.2.0, must be compiled with htslib enabled)
  - `fastqc` (version 0.12.1)
  - `trim_galore` (version 0.6.10)
  - `fastq_screen` (version 0.16.0, only required if running `fastq_screen`)
  - `bismark` (version 0.25.1, only required if running `fastq_screen`)
  - `pigz` (version 2.8)
  - `python` (version 3.13.3)
  - `pandas` (version 2.3.3)
  - `numpy` (version 2.3.5)
  - `matplotlib` (version 3.10.7)
  - `seaborn` (version 0.13.2)
  - `multiqc` (version 1.33)
  - `R` (version 4.4.3)
  - `tidyverse` (version 2.0.0, only required for plotting methylation controls)
  - `ggplot2` (version 3.5.2, only required for plotting methylation controls)
  - `patchwork` (version 1.3.0, only required for plotting methylation controls)
  - `viridislite` (version 0.4.2, only required for plotting methylation controls)

Two things of note, 1) it is easiest when working with `snakemake` to install `mamba` using `conda` when running with
`--use-conda`, and 2) it is preferable to install `snakemake` using `conda`, rather than using a module. This is due to
potential conflicts between packages (such as `matplotlib`) that can be found in the python distribution associated with
the snakemake module.

## Components of Workflow

The following components are generally in order, but may run in a different order, depending on exact dependencies
needed.
  - [default off] Generate asset files used during QC related rules
  - [default off] Modify and index reference genome to include methylation controls (lambda phage and pUC19)
  - [default off] Trim FASTQ files
  - [default off] Run Fastq Screen in bisulfite mode
  - Run FastQC on raw FASTQ files
  - Alignment, duplicate marking, and indexing of input data (biscuitSifter pipeline)
  - Samtools flagstat of input data
  - Methylation information extraction (BED Format)
  - Merge C and G beta values in CpG dinucleotide context
  - [default off] SNP and epiBED extraction
  - [default off] Run Preseq on aligned BAM
  - MultiQC with BICUIT QC modules specifically for methylation data
  - [default off] Generate plots of the observed / expected coverage ratio for different genomic features
  - [default off] Generate percentage of covered CpGs and CpG island coverage figures
  - [default off] Find coverage uniformity across genome
  - [default off] Plot percentage of genome covered
  - [default off] Find average methylation values in bins across genome
  - [default off] Find average methylation values in bins centered on specified regions
  - [default off] QC methylated and unmethylated controls
  - [default off] Generate bigWigs from BISCUIT BED files

Many options can be easily specified in the `config.yaml`! Otherwise, the commands in the Snakefile can also be modified
to meet different needs.

## Running the Workflow
For ease of reference, the configuration file `config/config.yaml` will be referred to throughout as the file to define
any configuration needed for your pipeline run. That said, you can copy this config file to another file and use that
config file in your pipeline with `snakemake --configfile /my/new/config.yaml` or by changing the `CONFIG_FILE` variable
in the SLURM submit script.

  - [Clone the repo](https://github.com/huishenlab/Biscuit_Snakemake_Workflow)
    - `git clone git@github.com:huishenlab/Biscuit_Snakemake_Workflow.git`
  - Place *gzipped* FASTQ files into `raw_data/`. Alternatively, you can specify the location of your *gzipped* FASTQ
  files in `config/config.yaml`.
  - Replace the example `config/samples.tsv` with your own sample sheet containing:
    - One row for each sample
    - The following three columns for each row (separated by a tab):
      - A. `sample` (name of the sample used throughout processing)
      - B. `fq1` (name of R1 file for `sample` in your raw data directory, multiple FASTQs can be specified with a
      comma-separated list)
      - C. `fq2` (name of R2 file for `sample` in your raw data directory, multiple FASTQs can be specified with a
      comma-separated list)
      - D. Any other columns included are ignored
    - Note, you can either edit `config/samples.tsv` in place or specify the path to your sample sheet in
    `config/config.yaml`. If you create your own sample sheet, make sure to include the header line as is seen in the
    example file.
  - Modify `config/config.yaml` to specify the appropriate
    - Reference genome
    - BISCUIT index
    - BISCUIT QC assets (see [Quality Control]({{ site.baseurl }}{% link docs/alignment/QC.md %}) for details)
    - Toggle optional workflow components
    - Set other run parameters in `config/config.yaml`
    - Turn on optional rules in `config/config.yaml` (change from False to True)
    - If you are using environmental modules on your system, you can set the locations in the corresponding location. By
    default, the pipeline will use `conda`/`mamba` to download the required packages. Note, if using the modules and a
    module is not available, snakemake gives a warning but will run successfully *as long as the required executables are
    in PATH*.
  - Modify SLURM submit script as needed (new config file in `CONFIG_FILE`, etc.).
  - Then submit the workflow to an HPC using something similar to `bin/run_snakemake_workflow.slurm` (e.g.,
  `sbatch bin/run_snakemake_workflow.slurm`). `bin/run_snakemake_workflow.slurm` works for a SLURM queue
  system. A PBS/Torque version is available in a previous release on GitHub for those who need it.

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
    - `snps/` SNP BED files (included if `generate_snps: True` in `config.yaml`)
    - `trim_reads` trimmed FASTQ files and FastQC reports
  - Log files can be found in the `logs/` directory
  - Benchmarking files can be found in the `benchmarks/` directory

## Example Dataset

The cloned Snakemake repository comes with a five sample test dataset to see how this workflow works on your system. To
run the test dataset, copy the ten `.fq.gz` files in `bin/working_example_dataset` into `raw_data/` and use the default
`bin/samples.tsv` file. This set of files should be mapped to the human genome.

## Useful Commands
For more information on Snakemake: <https://snakemake.readthedocs.io/en/stable/>

  - Perform a dry run of the commands that will be run by snakemake: `snakemake --dry-run`
  - Unlock the pipeline after a manually aborted run: `snakemake --unlock --cores 1`
  - Create a workflow diagram of your run: `snakemake --dag | dot -Tpng > my_dag.png`
  - Snakemake can also be run on the command line: `snakemake --use-conda --cores 1`
