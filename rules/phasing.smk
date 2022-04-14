# Find the variants from the mapped reads
rule medaka_variant:
    input:
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        temp1 = temp("{sample}/temp_files/medaka_variant/{regions}/round_0_hap_mixed_phased.vcf"),
        temp2 = temp("{sample}/temp_files/medaka_variant/{regions}/round_0_hap_mixed_unphased.vcf"),
        temp3 = temp("{sample}/temp_files/medaka_variant/{regions}/round_1_hap_1.vcf"),
        temp4 = temp("{sample}/temp_files/medaka_variant/{regions}/round_1_hap_2.vcf"),
        temp5 = temp("{sample}/temp_files/medaka_variant/{regions}/round_1_phased.vcf"),
        temp6 = temp("{sample}/temp_files/medaka_variant/{regions}/round_1_unfiltered.vcf"),
        temp7 = temp("{sample}/temp_files/medaka_variant/{regions}/round_1.vcf"),
        temp8 = temp("{sample}/temp_files/medaka_variant/{regions}/{sample}.bam"),
        temp9 = temp("{sample}/temp_files/medaka_variant/{regions}/{sample}.bam.bai"),
        temp10 = temp("{sample}/temp_files/medaka_variant/{regions}/round_0_hap_mixed_phased.bam"),
        temp11 = temp("{sample}/temp_files/medaka_variant/{regions}/round_0_hap_mixed_phased.vcf.gz"),
        temp12 = temp("{sample}/temp_files/medaka_variant/{regions}/round_0_hap_mixed_probs.hdf"),
        temp13 = temp("{sample}/temp_files/medaka_variant/{regions}/round_1_hap_2_probs.hdf"),
        temp14 = temp("{sample}/temp_files/medaka_variant/{regions}/round_0_hap_mixed_phased.bam.bai"),
        temp15 = temp("{sample}/temp_files/medaka_variant/{regions}/round_0_hap_mixed_phased.vcf.gz.tbi"),
        temp16 = temp("{sample}/temp_files/medaka_variant/{regions}/round_1_hap_1_probs.hdf")
    params:
        memory_per_thread="8G",
        run_time="2:0:0:0",
        project_name="simpsonlab",
        phased_output="-p"
    log:
        "{sample}/analysis/logs/temp_files/medaka_variant/{sample}.{regions}.medaka_variant.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/medaka_variant/{sample}.{regions}.medaka_variant.txt"
    threads: 8
    shell:
        "{medaka_dir}/medaka_variant {params.phased_output} -f {reference} -i {input.bam} -t {threads} \
        -o {wildcards.sample}/temp_files/medaka_variant/{wildcards.regions} -r {wildcards.regions} &> {log}"

# Zip the vcf for merging
rule zip_medaka_vcf:
    input:
        "{sample}/temp_files/medaka_variant/{regions}/round_1_phased.vcf"
    output:
        temp("{sample}/temp_files/medaka_variant/{regions}/round_1_phased.vcf.gz")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        project_name="simpsonlab"
    log:
        "{sample}/analysis/logs/temp_files/medaka_variant/{sample}.{regions}.zip_medaka_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/medaka_variant/{sample}.{regions}.zip_medaka_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/bgzip --threads {threads} -c {input} > {output} 2> {log}"
        
# Index the zipped vcf for merging
rule index_medaka_vcf:
    input:
        "{sample}/temp_files/medaka_variant/{regions}/round_1_phased.vcf.gz"
    output:
        temp("{sample}/temp_files/medaka_variant/{regions}/round_1_phased.vcf.gz.tbi")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        project_name="simpsonlab",
        file_type="vcf"
    log:
        "{sample}/analysis/logs/temp_files/medaka_variant/{sample}.{regions}.index_medaka_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/medaka_variant/{sample}.{regions}.index_medaka_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/tabix -p {params.file_type} {input} &> {log}"

# Merge the vcf files
rule merge_medaka_vcf:
    input:
        zip = expand("{{sample}}/temp_files/medaka_variant/{regions}/round_1_phased.vcf.gz", regions=regions),
        index = expand("{{sample}}/temp_files/medaka_variant/{regions}/round_1_phased.vcf.gz.tbi", regions=regions)
    output:
        protected("{sample}/mapped/{sample}.medaka.vcf.gz")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        project_name="simpsonlab",
        concat="-a -O v",
        sort="-O z"
    log:
        "{sample}/analysis/logs/mapped/{sample}.merge_medaka_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.merge_medaka_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/bcftools concat {params.concat} {input.zip} 2> {log} | \
        {conda_dir}/bcftools sort {params.sort} -o {output} &>> {log}"
        
# Index the merged vcf file
rule index_merged_medaka_vcf:
    input:
        "{sample}/mapped/{sample}.medaka.vcf.gz"
    output:
        protected("{sample}/mapped/{sample}.medaka.vcf.gz.tbi")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        project_name="simpsonlab",
        file_type="vcf"
    log:
        "{sample}/analysis/logs/mapped/{sample}.index_merged_medaka_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.index_merged_medaka_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/tabix -p {params.file_type} {input} &> {log}"
        
# Tag the bam file with haplotype information
rule haplotag_bam:
    input:
        vcf = "{sample}/mapped/{sample}.medaka.vcf.gz",
        tbi = "{sample}/mapped/{sample}.medaka.vcf.gz.tbi",
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        temp("{sample}/temp_files/whatshap/{regions}/{sample}.{regions}.phased.bam")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0",
        project_name="simpsonlab",
        RG="--ignore-read-groups"
    log:
        "{sample}/analysis/logs/temp_files/whatshap/{regions}/{sample}.{regions}.phased.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/whatshap/{regions}/{sample}.{regions}.phased.txt"
    threads: 1
    shell:
        "{conda_dir}/whatshap haplotag {params.RG} --regions {wildcards.regions} --reference {reference} \
        {input.vcf} {input.bam} 2> {log} | {conda_dir}/samtools sort -o {output} &>> {log}"

# Merge the phased bam files
rule merge_bam:
    input:
        expand("{{sample}}/temp_files/whatshap/{regions}/{{sample}}.{regions}.phased.bam", regions=regions)
    output:
        protected("{sample}/mapped/{sample}.phased.bam")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0",
        project_name="simpsonlab"
    log:
        "{sample}/analysis/logs/mapped/{sample}.merge_bam.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.merge_bam.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools merge {output} {input} &> {log}"
        
# Index the phased bam file
rule haplotag_bam_index:
    input:
        "{sample}/mapped/{sample}.phased.bam"
    output:
        protected("{sample}/mapped/{sample}.phased.bam.bai")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        project_name="simpsonlab"
    log:
        "{sample}/analysis/logs/mapped/{sample}.haplotag_index.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.haplotag_index.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools index {input} &> {log}"
       
#################################
# Variant calling with longshot #
#################################

# Find the variants from the mapped reads
rule longshot:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        temp("{sample}/temp_files/longshot/{sample}.{regions}.vcf")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0",
        project_name="simpsonlab",
        thresholds="-c 2 -C 200 -I 10000 -e 1"
    log:
        "{sample}/analysis/logs/temp_files/longshot/{sample}.{regions}.longshot.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/longshot/{sample}.{regions}.longshot.txt"
    threads: 1
    shell:
        "{conda_dir}/longshot {params.thresholds} -r {wildcards.regions} --bam {input.bam} \
        --ref {reference} --out {output} &> {log}"

# Zip the phased bams
rule zip_phased_bam:
    input:
        "{sample}/temp_files/longshot/{sample}.{regions}.vcf"
    output:
        temp("{sample}/temp_files/longshot/{sample}.{regions}.vcf.gz")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        project_name="simpsonlab"
    log:
        "{sample}/analysis/logs/temp_files/longshot/{sample}.{regions}.zip.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/longshot/{sample}.{regions}.zip.txt"
    threads: 1
    shell:
        "{conda_dir}/bgzip --threads {threads} -c {input} > {output} 2> {log}"
        
# Index the phased bams
rule index_phased_bam:
    input:
        "{sample}/temp_files/longshot/{sample}.{regions}.vcf.gz"
    output:
        temp("{sample}/temp_files/longshot/{sample}.{regions}.vcf.gz.tbi")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        project_name="simpsonlab",
        file_type="vcf"
    log:
        "{sample}/analysis/logs/temp_files/longshot/{sample}.{regions}.index.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/longshot/{sample}.{regions}.index.txt"
    threads: 1
    shell:
        "{conda_dir}/tabix -p {params.file_type} {input} &> {log}"
                
# Merge the vcf files
rule merge_longshot_vcf:
    input:
        zip = expand("{{sample}}/temp_files/longshot/{{sample}}.{regions}.vcf.gz", regions=regions),
        index = expand("{{sample}}/temp_files/longshot/{{sample}}.{regions}.vcf.gz.tbi", regions=regions)
    output:
        protected("{sample}/mapped/{sample}.longshot.vcf.gz")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        project_name="simpsonlab",
        concat="-a -O v",
        sort="-O z"
    log:
        "{sample}/analysis/logs/mapped/{sample}.merge_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.merge_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/bcftools concat {params.concat} {input.zip} 2> {log} | \
        {conda_dir}/bcftools sort {params.sort} -o {output} &>> {log}"
        
# Merge the vcf files
rule index_longshot_vcf:
    input:
        "{sample}/mapped/{sample}.longshot.vcf.gz"
    output:
        protected("{sample}/mapped/{sample}.longshot.vcf.gz.tbi")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        project_name="simpsonlab",
        file_type="vcf"
    log:
        "{sample}/analysis/logs/mapped/{sample}.index_longshot_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.index_longshot_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/tabix -p {params.file_type} {input} &> {log}"
               
# Extract the allele counts from the vcf
rule extract_longshot_vcf:
    input:
        zip = "{sample}/mapped/{sample}.longshot.vcf.gz",
        index = "{sample}/mapped/{sample}.longshot.vcf.gz.tbi"
    output:
        temp("{sample}/mapped/{sample}.longshot.tsv")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        project_name="simpsonlab",
        view_param="'PASS'",
        query_param="'%CHROM\t%POS\t%AC\n'"
    log:
        "{sample}/analysis/logs/mapped/{sample}.extract_longshot_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.extract_longshot_vcf.txt"
    threads: 1
    shell:
        "zcat {input.zip} 2> {log} | {conda_dir}/bcftools view -f {params.view_param} \
        2>> {log} | {conda_dir}/bcftools query -f {params.query_param} > {output} 2> {log}"
        
# Calcualte the allele counts from the vcf
rule calculate_longshot_b_allele_frequency:
    input:
        "{sample}/mapped/{sample}.longshot.tsv"
    output:
        protected("{sample}/mapped/{sample}.b_allele_frequency.longshot.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        project_name="simpsonlab"
    log:
        "{sample}/analysis/logs/mapped/{sample}.calculate_longshot_b_allele_frequency.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.calculate_longshot_b_allele_frequency.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/calculate_allele_frequency.py -i {input} -o {output} 2> {log}"
