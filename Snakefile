import os
import re

# To run all of the analysis on Univa:
#   snakemake --jobs 500 --keep-going --latency-wait 30 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log 
#   -q queue_name -P project_name -pe smp {threads} -l h_vmem={params.memory_per_thread} -l h_rt={params.run_time} -b y" run_all_analysis

configfile: "config.yaml"

# Get the path information for the programs and scripts
root_dir = os.getcwd()
conda_dir = config["conda_dir"]
scripts_dir = config["workflow_dir"] + "/scripts"
calls_header = config["workflow_dir"] + "/headers/nanopolish_call_header.tsv"
frequency_header = config["workflow_dir"] + "/headers/nanopolish_frequency_header.tsv"
dmr_header = config["workflow_dir"] + "/headers/dmr_header.tsv"

# Declare the variables and files needed for the reference genome
reference: config["workflow_dir"] + "/reference/GRCh38_no_alt_analysis_set.fna"
chromosome_sizes: config["workflow_dir"] + "/reference/GRCh38_chromosome_sizes.tsv"
high_frequency_kmers: config["workflow_dir"] + "/reference/GRCh38_high_frequency_kmers.txt"
genome_gaps: config["workflow_dir"] + "/reference/GRCh38_gaps.bed"
gene_list: config["workflow_dir"] + "/reference/GRCh38_genes.bed"
ctcf_sites: config["workflow_dir"] + "/reference/GRCh38_CTCF_binding_sites.bed"
genome_length: 3100000000
chromosomes: ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
               'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

#===============================================================================================
# Functions
#===============================================================================================

# Get the mean coverage for the dataset
def get_coverage(sample_name):
    file_name = str(sample_name) + "/analysis/nanoplot/" + str(sample_name) + "NanoStats.txt"
    if os.path.exists(file_name):
        f = open(file_name)
#        pattern = "number_of_bases"
        pattern = "Total bases"
        for line in f:
            if re.search(pattern, line):
#                total_bases = line.split()[1].strip()
                total_bases = line.split(":")[1].strip().replace(',', '')
                return float(total_bases)/3100000000

# Get the list of regions for the workflow
def get_regions():
    f = open(config["workflow_dir"] + "/reference/GRCh38_regions_list.txt")
    regions = list()
    for line in f:
        regions.append(line.rstrip())
    return regions
    
regions = get_regions()

# Include the rules for the workflow
include: config["workflow_dir"] + "/rules/analysis.smk"
include: config["workflow_dir"] + "/rules/assembly.smk"
include: config["workflow_dir"] + "/rules/dmrs.smk"
include: config["workflow_dir"] + "/rules/mapping.smk"
include: config["workflow_dir"] + "/rules/methylation.smk"
include: config["workflow_dir"] + "/rules/phasing.smk"
include: config["workflow_dir"] + "/rules/structural_variants.smk"

# Default rule to make everything
rule all_but_assembly:
    input:
        expand("{sample}/fastq/{sample}.fastq.gz", sample=config["samples"]),
        expand("{sample}/mapped/{sample}.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/methylation/{sample}.nanopolish_frequency.tsv", sample=config["samples"]),
        expand("{sample}/methylation/{sample}.nanopolish_calls.tsv.gz", sample=config["samples"]),
        expand("{sample}/mapped/{sample}.phased.methylated.bam.bai", sample=config["samples"]),
        expand("{tumor}/analysis/methylation/dmrs/{normal}/{tumor}.annotated_dmrs.bed", normal=config["normals"], tumor=config["tumors"]),
        expand("{sample}/analysis/structural_variants/{sample}.insertions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.deletions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.duplications.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.inversions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.translocations.bedpe", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.CNVs.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.1000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.500000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.small_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.large_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.{chromosome}.pdf", normal=config["normals"], sample=config["tumors"], chromosome=chromosomes)

rule filter_fastq:
    input:
        expand("{sample}/fastq/{sample}.fastq.gz", sample=config["samples"])

rule assembly:
    input:
        expand("{sample}/analysis/assembly/quast_flye/report.tsv", sample=config["samples"]),
        expand("{sample}/analysis/assembly/quast_racon/report.tsv", sample=config["samples"]),
        expand("{sample}/analysis/assembly/quast_medaka/report.tsv", sample=config["samples"])

rule mapping:
    input:
        expand("{sample}/analysis/coverage/{sample}.b_allele_frequency.longshot.bed", sample=config["samples"])
        
rule call_methylation:
    input:
        expand("{sample}/methylation/{sample}.nanopolish_frequency.tsv", sample=config["samples"]),
        expand("{sample}/methylation/{sample}.nanopolish_calls.tsv.gz", sample=config["samples"])

rule methylation_mapping:
    input:
        expand("{sample}/mapped/{sample}.phased.methylated.bam.bai", sample=config["samples"])

rule structural_variants:
    input:
        expand("{sample}/structural_variants/{sample}.insertions.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.deletions.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.duplications.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.inversions.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.translocations.bedpe", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.CNVs.bed", sample=config["samples"])
        
rule DMRs:
    input:
        expand("{tumor}/analysis/dmrs/{normal}/plots", normal=config["normals"], tumor=config["tumors"]),  
        expand("{tumor}/analysis/methylation/dmrs/{normal}/{tumor}.annotated_dmrs.bed", normal=config["normals"], tumor=config["tumors"])

rule run_SVs:
    input:
        expand("{sample}/structural_variants/{sample}.insertions.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.deletions.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.duplications.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.inversions.bed", sample=config["samples"]),
        expand("{sample}/structural_variants/{sample}.translocations.bedpe", sample=config["samples"])
        
rule run_SV_analysis:
    input:
        expand("{sample}/analysis/structural_variants/{sample}.insertions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.deletions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.duplications.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.inversions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.translocations.bedpe", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.CNVs.bed", sample=config["samples"])

rule run_coverage_analysis:
    input:
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.{chromosome}.pdf", normal=config["normals"], sample=config["samples"], chromosome=chromosomes),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.1000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.500000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.small_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.large_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.{chromosome}.pdf", normal=config["normals"], sample=config["tumors"], chromosome=chromosomes)

rule run_all_analysis:
    input:
        expand("{sample}/analysis/structural_variants/{sample}.insertions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.deletions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.duplications.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.inversions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.translocations.bedpe", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.CNVs.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.1000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.500000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.small_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.large_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.{chromosome}.pdf", normal=config["normals"], sample=config["tumors"], chromosome=chromosomes)

      
