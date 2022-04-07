# Assembly the genome with Flye
rule flye:
    input:
        "{sample}/fastq/{sample}.fastq.gz"
    output:
        assembly = protected("{sample}/assembly/assembly.fasta"),
        gfa = temp("{sample}/assembly/assembly_graph.gfa"),
        gv = temp("{sample}/assembly/assembly_graph.gv"),
        info = temp("{sample}/assembly/assembly_info.txt"),
        log = temp("{sample}/assembly/flye.log"),
        json = temp("{sample}/assembly/params.json")
    params:
        memory_per_thread="37G",
        run_time="14:0:0:0",
        genome_size="3.1g"
    log:
        "{sample}/analysis/logs/assembly/{sample}.flye.log"
    benchmark:
        "{sample}/analysis/benchmarks/assembly/{sample}.flye.txt"
    threads: 20
    shell:
        "{conda_dir}/flye --nano-raw {input} -o {wildcards.sample}/assembly \
        -g {params.genome_size} -t {threads} &> {log}"


# Map nanopore reads to the assembly
rule minimap2_racon:
    input:
        assembly = "{sample}/assembly/assembly.fasta",
        reads = "{sample}/fastq/{sample}.fastq.gz"
    output:
        temp("{sample}/assembly/{sample}.minimap2.paf")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        read_type="map-ont"
    log:
        "{sample}/analysis/logs/assembly/{sample}.minimap2_racon.log"
    benchmark:
        "{sample}/analysis/benchmarks/assembly/{sample}.minimap2_racon.txt"
    threads: 8
    shell:
        "{conda_dir}/minimap2 -x {params.read_type} -t {threads} \
        {input.assembly} {input.reads} > {output} 2> {log}"

# Polish assembly with racon
rule racon:
    input:
        assembly = "{sample}/assembly/assembly.fasta",
        reads = "{sample}/fastq/{sample}.fastq.gz",
        paf = "{sample}/assembly/{sample}.minimap2.paf"
    output:
        protected("{sample}/assembly/{sample}.racon.fasta"),
    params:
        memory_per_thread="25G",
        run_time="7:0:0:0",
        output_unpolished="-u"
    log:
        "{sample}/analysis/logs/assembly/{sample}.racon.log"
    benchmark:
        "{sample}/analysis/benchmarks/assembly/{sample}.racon.txt"
    threads: 20
    shell:
        "{conda_dir}/racon {params.output_unpolished} -t {threads} \
        {input.reads} {input.paf} {input.assembly} > {output} 2> {log}"

# Polish racon polished assembly with Medaka
rule medaka_consensus:
    input:
        assembly = "{sample}/assembly/{sample}.racon.fasta",
        reads = "{sample}/fastq/{sample}.fastq.gz"
    output:
        consensus = protected("{sample}/assembly/consensus.fasta"),
        bam = temp("{sample}/assembly/calls_to_draft.bam"),
        bai = temp("{sample}/assembly/calls_to_draft.bam.bai"),
        hdf = temp("{sample}/assembly/consensus_probs.hdf")
    params:
        memory_per_thread="31G",
        run_time="7:0:0:0",
        out_folder = "{sample}/assembly",
    log:
        "{sample}/analysis/logs/assembly/{sample}.medaka_consensus.log"
    benchmark:
        "{sample}/analysis/benchmarks/assembly/{sample}.medaka_consensus.txt"
    threads: 8
    shell:
        "{medaka_dir}/medaka_consensus -i {input.reads} -d {input.assembly} \
        -o {params.out_folder} -t {threads} &> {log}"

# Calculate Flye assembly statistics with QUAST
rule quast_flye:
    input:
        "{sample}/assembly/assembly.fasta",
    output:
        protected("{sample}/analysis/assembly/quast_flye/report.tsv")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        genome_size="--large",
        del_temp="--space-efficient",
        output_folder="{sample}/analysis/assembly/quast_flye"
    log:
        "{sample}/analysis/logs/analysis/assembly/quast_flye.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/assembly/quast_flye.txt"
    threads: 8
    shell:
        "{quast_dir}/python {quast_dir}/quast {params.genome_size} -r {reference} \
		-o {params.output_folder} -t {threads} {params.del_temp} {input} &> {log}"

# Calculate Racon polished assembly statistics with QUAST
rule quast_racon:
    input:
        "{sample}/assembly/{sample}.racon.fasta"
    output:
        protected("{sample}/analysis/assembly/quast_racon/report.tsv")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        genome_size="--large",
        del_temp="--space-efficient",
        output_folder="{sample}/analysis/assembly/quast_racon"
    log:
        "{sample}/analysis/logs/analysis/assembly/quast_racon.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/assembly/quast_racon.txt"
    threads: 8
    shell:
        "{quast_dir}/python {quast_dir}/quast {params.genome_size} -r {reference} \
		-o {params.output_folder} -t {threads} {params.del_temp} {input} &> {log}"

# Calculate Medaka polished assembly statistics with QUAST
rule quast_medaka:
    input:
        "{sample}/assembly/consensus.fasta"
    output:
        protected("{sample}/analysis/assembly/quast_medaka/report.tsv")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        genome_size="--large",
        del_temp="--space-efficient",
        output_folder="{sample}/analysis/assembly/quast_medaka"
    log:
        "{sample}/analysis/logs/analysis/assembly/quast_medaka.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/assembly/quast_medaka.txt"
    threads: 8
    shell:
        "{quast_dir}/python {quast_dir}/quast {params.genome_size} -r {reference} \
		-o {params.output_folder} -t {threads} {params.del_temp} {input} &> {log}"
