# Find DMR's with DSS
rule find_DMRs:
    input:
        normal = "{normal}/methylation/{normal}.nanopolish_frequency.tsv",
        tumor = "{tumor}/methylation/{tumor}.nanopolish_frequency.tsv"
    output:
        protected("{tumor}/methylation/dmrs/{normal}/{tumor}.dmrs.csv")
    params:
        memory_per_thread="15G",
        run_time="1:0:0:0",
        vanilla="--vanilla",
        p_threshold="1e-3",
        dis_merge="1000",
        minCG="5",
        minlen="50",
        smoothing_span="100"
    log:
        "{tumor}/analysis/logs/methylation/dmrs/{normal}/{tumor}.DSS.log"
    benchmark:
        "{tumor}/analysis/benchmarks/methylation/dmrs/{normal}/{tumor}.DSS.txt"
    threads: 16
    shell:
        "{conda_dir}/Rscript {params.vanilla} {scripts_dir}/dss.R {input.normal} \
        {input.tumor} {output} {threads} {params.p_threshold} {params.dis_merge} \
        {params.minCG} {params.minlen} {params.smoothing_span} &> {log}"

rule convert_dmrs:
    input:
        "{tumor}/methylation/dmrs/{normal}/{tumor}.dmrs.csv"
    output:
        protected("{tumor}/methylation/dmrs/{normal}/{tumor}.dmrs.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        awk="'NR>1'",
        cut="-d, -f1 --complement",
        tr="',' '\\t'",
        sed="'s/\"//g'"
    log:
        "{tumor}/analysis/logs/methylation/dmrs/{normal}/{tumor}.convert_dmrs.log"
    benchmark:
        "{tumor}/analysis/benchmarks/methylation/dmrs/{normal}/{tumor}.convert_dmrs.txt"
    threads: 1
    shell:
        "awk {params.awk} {input} 2> {log} | cut {params.cut} 2>> {log} | \
        tr {params.tr} 2>> {log} | sed {params.sed} 2>> {log} | \
        {conda_dir}/bedtools sort -i > {output} 2>> {log}"        

rule annotate_dmrs:
    input:
        "{tumor}/methylation/dmrs/{normal}/{tumor}.dmrs.bed"
    output:
        protected("{tumor}/analysis/methylation/dmrs/{normal}/{tumor}.annotated_dmrs.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        write="-loj"
    log:
        "{tumor}/analysis/logs/analysis/methylation/dmrs/{normal}/{tumor}.annotate_dmrs.log"
    benchmark:
        "{tumor}/analysis/benchmarks/analysis/methylation/dmrs/{normal}/{tumor}.annotate_dmrs.txt"
    threads: 1
    shell:
        "{conda_dir}/bedtools intersect {params.write} -a {input} -b {genes} > {output} 2>> {log}"

rule dmr_regions:
    input:
        "{tumor}/methylation/dmrs/{normal}/{tumor}.filtered_dmrs.bed"
    output:
        protected("{tumor}/analysis/methylation/dmrs/{normal}/{tumor}.dmr_regions.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        options="-counts",
        names="genes promoters exons introns 5_utr 3_utr coding_exon ctcf LINEs SINEs repeats"
    log:
        "{tumor}/analysis/logs/analysis/methylation/dmrs/{normal}/{tumor}.dmr_regions.log"
    benchmark:
        "{tumor}/analysis/benchmarks/analysis/methylation/dmrs/{normal}/{tumor}.annotate_dmrs.txt"
    threads: 1
    shell:
        "{conda_dir}/bedtools annotate {params.options} -i {input} -names {params.names} -files {genes} {promoters} \
        {exons} {introns} {utr_5} {utr_3} {coding_exon} {ctcf} {LINEs} {SINEs} {repeats} > {output} 2>> {log}"

rule merge_dmrs_to_plot:
    input:
        "{tumor}/methylation/dmrs/{normal}/{tumor}.filtered_dmrs.bed"
    output:
        protected("{tumor}/methylation/dmrs/{normal}/{tumor}.merged_dmrs.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        awk="'$5>20'",
        merge_length="30000"
    log:
        "{tumor}/analysis/logs/analysis/methylation/dmrs/{normal}/{tumor}.merge_dmrs_to_plot.log"
    benchmark:
        "{tumor}/analysis/benchmarks/analysis/methylation/dmrs/{normal}/{tumor}.merge_dmrs_to_plot.txt"
    threads: 1
    shell:
        "awk {params.awk} {input} 2> {log} | {conda_dir}/bedtools merge -i stdin \
        -d {params.merge_length} > {output} 2>> {log}"
                
rule plot_dmrs:
    input:
        normal = "{normal}/methylation/{normal}.nanopolish_frequency.tsv",
        tumor = "{tumor}/methylation/{tumor}.nanopolish_frequency.tsv",
        dmrs = "{tumor}/methylation/dmrs/{normal}/{tumor}.filtered_dmrs.bed",
        merged_dmrs = "{tumor}/methylation/dmrs/{normal}/{tumor}.merged_dmrs.bed"
    output:
        directory("{tumor}/analysis/methylation/dmrs/{normal}/plots")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0",
        slop="15000"
    log:
        "{tumor}/analysis/logs/analysis/methylation/dmrs/{normal}/plots/{tumor}.plot_dmrs.log"
    benchmark:
        "{tumor}/analysis/benchmarks/analysis/methylation/dmrs/{normal}/plots/{tumor}.plot_dmrs.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/karyoploter_methylation.R {input.normal} {input.tumor} \
        {input.dmrs} {ctcf} {input.merged_dmrs} {wildcards.normal} {wildcards.tumor} \
        {params.slop} {root_dir}/{output} &> {log}"

