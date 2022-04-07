# Filter out low quality reads, short reads, and really long reads
rule nanofilt:
    input:
        "{sample}/fastq/{sample}.fastq"
    output:
        protected("{sample}/fastq/{sample}.fastq.gz")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        project_name="simpsonlab",
        quality="7",
        min_length="500",
        max_length="500000"
    log:
        "{sample}/analysis/logs/fastq/{sample}.nanofilt.log"
    benchmark:
        "{sample}/analysis/benchmarks/fastq/{sample}.nanofilt.txt"
    threads: 1
    shell:
        "{conda_dir}/NanoFilt --logfile {log} -q {params.quality} -l {params.min_length} \
        --maxlength {params.max_length} {input} | {conda_dir}/bgzip -c > {output} 2>> {log}"

# Index the zipped fastq file
rule bgzip_index:
    input:
        "{sample}/fastq/{sample}.fastq.gz"
    output:
        "{sample}/fastq/{sample}.fastq.gz.gzi"
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        project_name="simpsonlab"
    log:
        "{sample}/analysis/logs/fastq/{sample}.bgzip_index.log"
    benchmark:
        "{sample}/analysis/benchmarks/fastq/{sample}.bgzip_index.txt"
    threads: 1
    shell:
        "{conda_dir}/bgzip -r {input} 2> {log}"

# Map the sample to the reference genome
rule winnowmap:
    input:
        fastq = "{sample}/fastq/{sample}.fastq.gz",
        index = "{sample}/fastq/{sample}.fastq.gz.gzi"
    output:
        temp("{sample}/mapped/{sample}.bam")
    params:
        memory_per_thread="8G",
        run_time="2:0:0:0",
        project_name="simpsonlab",
        preset_options="-ax map-ont",
        include_MD_tag="--MD"
    log:
        "{sample}/analysis/logs/mapped/{sample}.winnowmap.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.winnowmap.txt"
    threads: 12
    shell:
        "{conda_dir}/winnowmap {params.preset_options} {params.include_MD_tag} \
        -t {threads} -W {high_frequency_kmers} {reference} {input.fastq} 2> {log} | \
        {conda_dir}/samtools sort -o {output} &>> {log}"

# Index the mapped reads
rule winnowmap_index:
    input:
        "{sample}/mapped/{sample}.bam"
    output:
        temp("{sample}/mapped/{sample}.bam.bai")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        project_name="simpsonlab"
    log:
        "{sample}/analysis/logs/mapped/{sample}.winnowmap_index.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.winnowmap_index.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools index {input} &> {log}"
        
       
# Find reads with low mapping quality
rule find_low_mapped_reads:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        temp("{sample}/mapped/{sample}.low_mapped_reads.sam")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        quality="30",
        max="255"
    log:
        "{sample}/analysis/logs/mapped/{sample}.find_low_mapped_reads.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.find_low_mapped_reads.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools view -h {input.bam} 2> {log} | \
        awk '$5<={params.quality} || $5>={params.max} {{print $0}}' > {output} 2>> {log}"

# Find the coverage of regions with low mapping quality reads
rule find_low_mapping_coverage:
    input:
        "{sample}/mapped/{sample}.low_mapped_reads.sam"
    output:
        temp("{sample}/mapped/{sample}.low_mapped_reads.cov"),
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        grep="-v -E 'random|chrUn|alt|fix|chrM|chrEBV'"
    log:
        "{sample}/analysis/logs/mapped/{sample}.find_low_map_regions.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.find_low_map_regions.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools view -S -b -h {input} 2> {log}| \
        {conda_dir}/samtools depth - 2>> {log} | grep {params.grep} > {output} 2>> {log}"
      
# Make a bed file of the regions with low mapping quality      
rule make_low_mapped_bed:
    input:
        "{sample}/mapped/{sample}.low_mapped_reads.cov"
    output:
        protected("{sample}/mapped/{sample}.low_mapped_reads.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        slop="500",
        merge_length="10000",
        distance="100",
        min_coverage="2"
    log:
        "{sample}/analysis/logs/mapped/{sample}.make_low_mapped_bed.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.make_low_mapped_bed.txt"
    threads: 1
    shell:
        "{conda_dir}/SURVIVOR bincov {input} {params.distance} {params.min_coverage} 2> {log} | \
        {conda_dir}/bedtools slop -i stdin -g {chromosome_sizes} -b {params.slop} 2>> {log} | \
        {conda_dir}/bedtools merge -i stdin -d {params.merge_length} > {output} 2>> {log}"
