# Index the fastq file with Nanopolish
rule nanopolish_index:
    input:
        "{sample}/fastq/{sample}.fastq.gz"
    output:
        protected("{sample}/fastq/{sample}.fastq.gz.index"),
        protected("{sample}/fastq/{sample}.fastq.gz.index.readdb"),
        protected("{sample}/fastq/{sample}.fastq.gz.index.fai"),
        protected("{sample}/fastq/{sample}.fastq.gz.index.gzi")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0",
        fast5 = lambda wildcards: config[wildcards.sample]["fast5"],
        summary = lambda wildcards: config[wildcards.sample]["summary"]        
    log:
        "{sample}/analysis/logs/fastq/{sample}.nanopolish_index.log"
    benchmark:
        "{sample}/analysis/benchmarks/fastq/{sample}.nanopolish_index.txt"
    threads: 1
    shell:
        "{conda_dir}/nanopolish index {params.fast5} {params.summary} {input} &> {log}"

# Methylation calling with Nanopolish
rule nanopolish_call_methylation:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai",
        fastq = "{sample}/fastq/{sample}.fastq.gz",
        index = "{sample}/fastq/{sample}.fastq.gz.index",
        fai = "{sample}/fastq/{sample}.fastq.gz.index.fai",
        readdb = "{sample}/fastq/{sample}.fastq.gz.index.readdb",
        gzi = "{sample}/fastq/{sample}.fastq.gz.index.gzi"
    output:
        temp("{sample}/temp_files/methylation_calls_split/{sample}.{regions}.tsv")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0"
    log:
        "{sample}/analysis/logs/temp_files/methylation_calls_split/{sample}.{regions}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/methylation_calls_split/{sample}.{regions}.txt"
    threads: 4
    shell:
        "{conda_dir}/nanopolish call-methylation -t {threads} -b {input.bam} \
        -r {input.fastq} -g {reference} -w {wildcards.regions} > {output} 2> {log}"

# Calculate frequency from Nanopolish call files
rule nanopolish_methylation_frequency:
    input:
        "{sample}/temp_files/methylation_calls_split/{sample}.{regions}.tsv"
    output:
        temp("{sample}/temp_files/methylation_frequency_split/{sample}.{regions}.tsv")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        split_groups="-s"
    log:
        "{sample}/analysis/logs/temp_files/methylation_frequency_split/{sample}.{regions}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/methylation_frequency_split/{sample}.{regions}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/calculate_methylation_frequency.py \
        {params.split_groups} {input} > {output} 2> {log}"

# Concatenate the methylation call files
rule concatenate_call_files:
    input:
        expand("{{sample}}/temp_files/methylation_calls_split/{{sample}}.{regions}.tsv", regions=regions)
    output:
        temp_file = temp("{sample}/methylation/{sample}.temp_calls.tsv"),
        final_out = protected("{sample}/methylation/{sample}.nanopolish_calls.tsv.gz")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        awk=lambda wildcards: "'FNR==1 {next;} {print}'",
        sort="-k 1,1 -k2,2n -k3,3n"
    log:
        "{sample}/analysis/logs/methylation/{sample}.concatenate_calls.log"
    benchmark:
        "{sample}/analysis/benchmarks/methylation/{sample}.concatenate_calls.txt"
    threads: 1
    shell:
        """
        awk {params.awk} {input} 2> {log} | sort {params.sort} > {output.temp_file} 2>> {log}
        cat {calls_header} {output.temp_file} 2>> {log} | gzip > {output.final_out} 2>> {log}
        """

# Concatenate the frequency files
rule concatenate_frequency_files:
    input:
        expand("{{sample}}/temp_files/methylation_frequency_split/{{sample}}.{regions}.tsv", regions=regions)
    output:
        temp_file = temp("{sample}/methylation/{sample}.temp_frequency.tsv"),
        final_out = protected("{sample}/methylation/{sample}.nanopolish_frequency.tsv")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        awk=lambda wildcards: "'FNR==1 {next;} {print}'",
        sort="-k 1,1 -k2,2n -k3,3n"
    log:
        "{sample}/analysis/logs/methylation/{sample}.concatenate_frequency.log"
    benchmark:
        "{sample}/analysis/benchmarks/methylation/{sample}.concatenate_frequency.txt"
    threads: 1
    shell:
        """
        awk {params.awk} {input} 2> {log} | sort {params.sort} > {output.temp_file} 2>> {log}
        cat {frequency_header} {output.temp_file} > {output.final_out} 2>> {log}
        """

# Convert the call files to a bedGraph
rule make_bedGraph:
    input:
        "{sample}/temp_files/methylation_calls_split/{sample}.{regions}.tsv"
    output:
        temp("{sample}/temp_files/bedGraph_split/{sample}.{regions}.bed")
    params:
        memory_per_thread="8G",
        run_time="0:2:0:0",
        sort="-k 1,1 -k2,2n -k3,3n"
    log:
        "{sample}/analysis/logs/temp_files/bedGraph_split/{sample}.{regions}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/bedGraph_split/{sample}.{regions}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/mtsv2bedGraph.py -i {input} 2> {log} | \
        sort {params.sort} > {output} 2>> {log}"

rule make_bedGraph_index:
    input:
        "{sample}/temp_files/bedGraph_split/{sample}.{regions}.bed"
    output:
        gz = temp("{sample}/temp_files/bedGraph_split/{sample}.{regions}.bed.gz"),
        tbi = temp("{sample}/temp_files/bedGraph_split/{sample}.{regions}.bed.gz.tbi")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        std_out="-c",
        out_type="bed"
    log:
        "{sample}/analysis/logs/temp_files/bedGraph_split/{sample}.{regions}.index.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/bedGraph_split/{sample}.{regions}.index.txt"
    threads: 4
    shell:
        """
        {conda_dir}/bgzip {params.std_out} --threads {threads} {input} > {output.gz} 2> {log}
        {conda_dir}/tabix -p {params.out_type} {output.gz} &>> {log}
        """

# Convert the bedGraph files to a bam file
rule bedGraph_to_bam:
    input:
        bed = "{sample}/temp_files/bedGraph_split/{sample}.{regions}.bed.gz",
        tbi = "{sample}/temp_files/bedGraph_split/{sample}.{regions}.bed.gz.tbi",
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        temp("{sample}/temp_files/converted_bam_split/{sample}.{regions}.bam")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0"
    log:
        "{sample}/analysis/logs/temp_files/converted_bam_split/{sample}.{regions}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/converted_bam_split/{sample}.{regions}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/convert_bam.py -b {input.bam} \
        -c {input.bed} -o {output} -w {wildcards.regions} &> {log}"

# Merge the converted bam files
rule merge_converted_bam:
    input:
        expand("{{sample}}/temp_files/converted_bam_split/{{sample}}.{regions}.bam", regions=regions)
    output:
        protected("{sample}/mapped/{sample}.phased.methylated.bam")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0"
    log:
        "{sample}/analysis/logs/mapped/{sample}.merge_converted_bam.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.merge_converted_bam.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools merge {output} {input} &> {log}"

# Index the converted bam files
rule index_converted_bam:
    input:
        "{sample}/mapped/{sample}.phased.methylated.bam"
    output:
        protected("{sample}/mapped/{sample}.phased.methylated.bam.bai")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/mapped/{sample}.index_converted_bam.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.index_converted_bam.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools index {input} &> {log}"
