#rule cutadapt:
#    input:
#        "results/reads/{sample}_R1_001.fastq.gz"
#    output:
#        fastq="results/trimmed/{sample}.fastq.gz",
#        qc="results/trimmed/{sample}.qc.txt"
#    params:
#        adapters='-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AATGATACGGCGACCACCGAGATCTACAC',
#        others='--nextseq-trim=20 -m 10'
#    log:
#        "logs/cutadapt/{sample}.log"
#    threads: 16
#    script:
#        "0.66.0/bio/cutadapt/se"

rule trimmomatic:
    input:
        "results/reads/{sample}_R1_001.fastq.gz"
    output:
        "results/trimmed/{sample}.fastq.gz",
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
    threads:
        32
    wrapper:
        "0.67.0/bio/trimmomatic/se"

rule bwa_mem:
    input:
        reads="results/trimmed/{sample}.fastq.gz",
        index=lambda w: config["index_path"].format(species=w.species)
    output:
        "results/{species}/bwa_mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{species}/{sample}.log"
    params:
        index=lambda w, input: input.index,
        sort="samtools",
    threads: 16
    wrapper:
        "0.49.0/bio/bwa/mem"

rule filter:
    input:
        "results/{species}/%s_mapped/{sample}.bam" % mapper
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
