# Nanopore-Workflow
Snakemake workflow to process and filter long read data from Oxford Nanopore Technologies.  It is designed to compare whole human genome tumor/normal pairs of data sets, but can also run individual samples.  Reports and plots are genereated for differentially methylated regions, copy number variants, and structural variants.  Filtering heuristics typically reduce the reported translocations to the break points. It is suggested to have at least 15x - 20x of coverage, and a median read length of at least 5kbp - 6kbp.

## Installation instructions
Download the latest code from GitHub, and change to the workflow directory:

```
git clone https://github.com/mike-molnar/nanopore-workflow.git
cd nanopore-workflow
```

Copy the `Snakefile` and `config.yaml` files to the directory that you want to run the workflow:

```
cp Snakefile config.yaml /path/to/samples
```

Modify the `config.yaml` to point the necessary files and directories. The workflow is currently designed to have a single FASTQ, and a single sequencing summary file in a folder named `fastq` that is within a folder named after the sample.  The `config.yaml` file provides two examples of how to format the initial files and directories before running the workflow.
