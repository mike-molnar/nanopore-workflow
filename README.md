# Nanopore-Workflow
Snakemake workflow to process and filter long read data from Oxford Nanopore Technologies.  It is designed to compare whole human genome tumor/normal pairs, but can also run individual samples.  Reports and plots are generated for de novo genome assembly, differentially methylated regions, copy number variants, and structural variants.  Filtering heuristics typically reduce the reported translocations to the break points. It is suggested to have at least 15x - 20x of coverage, and a median read length of at least 5kbp - 6kbp.

![nanopore_workflow](https://user-images.githubusercontent.com/39533525/162601899-af7a5476-ced0-49a0-8108-71e8df757839.png)

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

Modify the `config.yaml` file to represent the information for the necessary files and directories of your sample(s). The workflow is currently designed to have a single FASTQ, and a single sequencing summary file in a folder named `fastq` that is in a folder named after the sample.  The `config.yaml` file provides an example of how to format the initial files and directories before running the workflow.

## Dependencies

There are many dependencies, so it is best to create a new Conda environment using the `environment.yml` file in the main directory of the workflow:

```
conda env create -n workflow -f environment.yml
conda activate workflow
```
