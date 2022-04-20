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

## To run on a grid engine

There are a few different grid engines, so the exact format may be different for your particular grid engine.  To run everything except the de novo assembly on a Univa grid engine:

```
snakemake --jobs 500 --rerun-incomplete --keep-going --latency-wait 30 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log -q queue_name -P project_name -pe smp {threads} -l h_vmem={params.memory_per_thread} -l h_rt={params.run_time} -b y" all_but_assembly
```

You will have to replace `queue_name` and `project_name` with the necessary values to run on your grid.

## Dependencies

There are many dependencies, so it is best to create a new Conda environment using the YAML files in the `env` directory.  There is a YAML file for the workflow, and another for Medaka.  You will need to install a separate environment for QUAST if you are going to run the de novo assembly portion of the workflow. Change to the `env` directory and create the environments with Conda:

```
cd /path/to/nanopore-worflow/env
conda env create -n nanopore-workflow -f nanopore-workflow_env.yml
conda env create -n medaka -f medaka_env.yml
conda env create -n quast -f quast_env.yml
conda env create -n R_env -f R_env.yml
conda activate nanopore-workflow
```

Before running the workflow you will need to `export` the paths of the four environments to your `PATH` variable:

```
export PATH="/path/to/conda/envs/nanopore-workflow/bin:$PATH"
export PATH="/path/to/conda/envs/medaka/bin:$PATH"
export PATH="/path/to/conda/envs/quast/bin:$PATH"
export PATH="/path/to/conda/envs/R_env/bin:$PATH"
```

### nanopore-workflow dependencies:
- bcftools
- bedtools
- cutesv
- flye
- longshot
- nanofilt v2.8.0
- nanoplot v1.20.0
- nanopolish
- seaborn v0.10.0
- snakemake
- sniffles
- survivor
- svim
- whatshap
- winnowmap

### R_env dependencies:
- bioconductor-karyoploter
- bioconductor-txdb.hsapiens.ucsc.hg38.knowngene
- bioconductor-org.hs.eg.db
- bioconductor-dss
- r-tidyverse
