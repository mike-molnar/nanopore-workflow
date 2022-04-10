# Nanopore-Workflow
Snakemake workflow to process and filter long read data from Oxford Nanopore Technologies.  It is designed to compare whole human genome tumor/normal pairs of data sets, but can also run individual samples.  Reports and plots are genereated for differentially methylated regions, copy number variants, and structural variants.  Filtering heuristics typically reduce the reported translocations to the break points. It is suggested to have at least 15x - 20x of coverage, and a median read length of at least 5kbp - 6kbp.

## Installation instructions
Download the latest code from GitHub:
`git clone https://github.com/mike-molnar/nanopore-workflow.git`
