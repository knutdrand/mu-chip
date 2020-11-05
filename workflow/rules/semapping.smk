rule cutadapt:
    input:
        "results/reads/{sample}_R1_001.fastq.gz"
    output:
        fastq="results/trimmed/{sample}.fastq.gz",
        qc="results/trimmed/{sample}.qc.txt"
    params:
        adapters='-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AATGATACGGCGACCACCGAGATCTACAC',
        others='--nextseq-trim=20 -m 10'
    log:
        "logs/cutadapt/{sample}.log"
    threads: 16
    script:
        "0.66.0/bio/cutadapt/se"

rule bwa_mem:
    input:
        reads="results/trimmed/{sample}.fastq.gz"
    output:
        "results/{species}/mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{species}/{sample}.log"
    params:
        index=config["data_dir"]+"{species}/{species}.fa.gz",
        sort="samtools",
    threads: 16
    wrapper:
        "0.49.0/bio/bwa/mem"

rule filter:
    input:
        "results/{species}/mapped/{sample}.bam"
    output:
        "results/{species}/mapped_filtered/{sample}.bam"
    params:
        "-Bb -q 30" 
    wrapper:
        "0.50.4/bio/samtools/view"

rule remove_duplicates:
    input:
        "results/{species}/mapped_filtered/{sample}.bam"
    output:
        bam=temp("results/{species}/dedup/{sample}.bam"),
        metrics="results/picard_dedup/{species}/{sample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{species}/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "0.64.0/bio/picard/markduplicates"

rule bamtobed:
    input:
        "results/{species}/dedup/{sample}.bam",
    output:
        "results/{species}/dedup_bed/{sample}.bed.gz",
    conda:
        "../envs/bamtobed.yaml"
    shell:
        "samtools view -b -f 64 {input} | bedtools bamtobed -i - | gzip > {output}"

rule merge_reads:
    input:
        lambda w: expand_se_combo("results/{species}/dedup_bed/{sample}.bed.gz", w)
    output:
        "results/{species}/merged/{celltype}_{condition}.bed.gz",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "zcat {input} | gzip | bedsort > {output}"
