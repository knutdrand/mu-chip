rule cutadapt_se:
    input:
        "results/reads/{sesample}.fastq.gz"
    output:
        fastq="results/trimmed/{sesample}.fastq.gz",
        qc="results/trimmed/{sesample}.qc.txt"
    params:
        adapters='-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AATGATACGGCGACCACCGAGATCTACAC'
        others='--nextseq-trim=20 -m 10'
    log:
        "logs/cutadapt/{sample}.log"
    threads: 16
    script:
        "0.66.0/bio/cutadapt/se"

rule bwa_mem_se:
    input:
        reads="trimmed/{sample}.fastq.gz"
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

rule filter_se:
    input:
        "{folder}mapped/{sample}.bam"
    output:
        "{folder}mapped_filtered/{sample}.bam"
    params:
        "-Bb -q 30" 
    wrapper:
        "0.50.4/bio/samtools/view"

rule remove_duplicates:
    input:
        "results/{species}/mapped_filtered/{sesample}.bam"
    output:
        bam=temp("results/{species}/dedup/{sesample}.bam"),
        metrics="results/picard_dedup/{species}/{sesample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{species}/{sesample}.log"
    params:
        "REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "0.64.0/bio/picard/markduplicates"

rule bamtobedpe2se:
    input:
        "{species}/dedup/{sesample}.bam",
    output:
        "{species}/dedup_bed/{sesample}.bed.gz",
    shell:
        "samtools view -b -f 64 {input} | bedtools bamtobed -i - | gzip > {output}"

rule genomecov_se:
    input:
        bed="results/{species}/dedup_bed/{sesample}.bed.gz",
        ref="results/{species}/data/chrom.sizes.txt"
    output:
        "results/{species}/coverage/{sesample}.bdg"
    log:
        "logs/coverage/{species}/{sesample}.log"
    params: 
        "-bga"
    wrapper:
        "0.64.0/bio/bedtools/genomecov"
