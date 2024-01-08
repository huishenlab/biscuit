---
title: Snakemake Pipeline
nav_order: 9
---

# Snakemake Pipeline

A Snakemake pipeline is [available on GitHub](https://github.com/huishenlab/Biscuit_Snakemake_Workflow) for streamlining
analyses.

## Dependencies

The following dependencies are downloaded when running with `--use-conda`, otherwise you must have these in your PATH.
| Package        | Conda Version Downloaded | Notes                                           |
|:---------------|:------------------------:|:------------------------------------------------|
| `snakemake`    | 7.0+                     | Needed before running pipeline                  |
| `biscuit`      | 1.2.0                    |                                                 |
| `htslib`       | 1.17                     |                                                 |
| `samtools`     | 1.17                     |                                                 |
| `dupsifter`    | 1.2.0                    |                                                 |
| `parallel`     | 20230322                 |                                                 |
| `bedtools`     | 2.30.0                   |                                                 |
| `preseq`       | 3.2.0                    | Must be compiled with htslib enabled            |
| `fastqc`       | 0.12.1                   |                                                 |
| `trim_galore`  | 0.6.10                   |                                                 |
| `fastq_screen` | 0.15.3                   | Only required if running `fastq_screen`)        |
| `bismark`      | 0.24.0                   | Only required if running `fastq_screen`)        |
| `pigz`         | 2.6                      |                                                 |
| `python`       | 3.11.3                   |                                                 |
| `pandas`       | 2.0.0                    |                                                 |
| `numpy`        | 1.24.2                   |                                                 |
| `matplotlib`   | 3.7.1                    |                                                 |
| `seaborn`      | 0.12.2                   |                                                 |
| `multiqc`      | 1.14                     |                                                 |
| `R`            | 4.2.3                    |                                                 |
| `tidyverse`    | 2.0.0                    | Only required for plotting methylation controls |
| `ggplot2`      | 3.4.2                    | Only required for plotting methylation controls |
| `patchwork`    | 1.1.2                    | Only required for plotting methylation controls |
| `viridislite`  | 0.4.1                    | Only required for plotting methylation controls |

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

Many options can be easily specified in the `config.yaml`! Otherwise, the commands in the Snakefile can also be modified
to meet different needs.

## Running the Workflow
For ease of reference, the configuration file `config/config.yaml` will be referred to throughout as the file to define
any configuration needed for your pipeline run. That said, you can copy this config file to another file and use that
config file in your pipeline with `snakemake --configfile /my/new/config.yaml` or by changing the `CONFIG_FILE` variable
in the SLURM submit script.

  - [Clone the repo](https://github.com/huishenlab/Biscuit_Snakemake_Workflow/tree/master)
    - SSH: `git clone git@github.com:huishenlab/Biscuit_Snakemake_Workflow.git`
    - HTTPS: `git clone https://github.com/huishenlab/Biscuit_Snakemake_Workflow.git`
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
    - `snps/` SNP BED files (included if generate_snps: 1 in `config.yaml`)
    - `trim_reads` trimmed FASTQ files and FastQC reports
  - Log files can be found in the `logs/` directory
  - Benchmarking files can be found in the `benchmarks/` directory

## Example Dataset

The cloned Snakemake repository comes with a five sample test dataset to see how this workflow works on your system. To
run the test dataset, copy the ten `.fq.gz` files in `bin/working_example_dataset` into `raw_data/` and use the default
`bin/samples.tsv` file. This set of files should be mapped to the human genome.

## Useful Commands
For more information on Snakemake: https://snakemake.readthedocs.io/en/stable/

  - Perform a dry run of the commands that will be run by snakemake: `snakemake -npr`
  - Unlock the pipeline after a manually aborted run: `snakemake --unlock --cores 1`
  - Create a workflow diagram of your run: `snakemake --dag | dot -Tpng > my_dag.png`
  - Snakemake can also be run on the command line: `snakemake --use-conda --cores 1`
