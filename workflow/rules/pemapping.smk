rule cutadapt_pe:
    input:
        "results/reads/{pesample}_R1_001.fastq.gz",
        "results/reads/{pesample}_R2_001.fastq.gz",
    output:
        temp(fastq1="results/trimmed/{pesample}_R1.fastq.gz"),
        temp(fastq2="results/trimmed/{pesample}_R2.fastq.gz"),
        qc="results/trimmed/{pesample}.qc.txt"
    params:
        adapters = '-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AATGATACGGCGACCACCGAGATCTACAC',
        others = '--nextseq-trim=20 -m 10'
    log:
        "logs/cutadapt/{pesample}.log"
    threads: 16
    wrapper:
        "0.66.0/bio/cutadapt/pe"

rule bwa_mem_pe:
    input:
        reads=expand("results/trimmed/{{pesample}}_R{read}.fastq.gz", read=[1,2])
    output:
        temp("results/{species}/mapped_pe/{pesample}.bam")
    log:
        "logs/bwa_mem/{species}/{pesample}.log"
    params:
        index=config["data_dir"]+"{species}/{species}.fa.gz",
        sort="samtools",
        sort_order="coordinate"
    threads:
        16
    wrapper:
        "0.49.0/bio/bwa/mem"

rule remove_duplicates:
    input:
        "results/{species}/mapped_pe/{pesample}.bam"
    output:
        bam=temp("results/{species}/dedup_pe/{pesample}.bam"),
        metrics="results/picard_dedup/{species}/{pesample}.metrics.txt"
    log:
        "logs/picard/picard_dedup/{species}/{pesample}.log"
    params:
        "REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
    wrapper:
        "0.64.0/bio/picard/markduplicates"

rule filter_reads:
    input:
        "results/{species}/dedup_pe/{pesample}.bam"
    output:
        temp("results/{species}/filtered_dedup_pe/{pesample}.bam")
    params:
        "-Bb -q %s -F 1796 -f 2" % config.get("mapq", "30")
    wrapper:
        "0.50.4/bio/samtools/view"

rule bamtobedpe:
    input:
        "results/{species}/filtered_dedup_pe/{pesample}.bam",
    output:
        temp("results/{species}/all_dedup_bed/{pesample}.bed"),
    shell:
        "samtools collate -Ou {input} | bedtools bamtobed -i - | chiptools pairbed | bedtools sort > {output}"

rule remove_small_fragments:
    input:
        "results/{species}/all_dedup_bed/{pesample}.bed",
    output:
        "results/{species}/dedup_bed/{pesample}.bed",
    shell:
        "awk '{{if (($3-$2)>50) print}}' {input} > {output}"

rule genomecov:
    input:
        bed="results/{species}/dedup_bed/{pesample}.bed",
        ref="results/{species}/data/chrom.sizes.txt"
    output:
        "results/{species}/coverage/{pesample}.bdg"
    log:
        "logs/coverage/{species}/{pesample}.log"
    params: 
        "-bga"
    wrapper:
        "0.64.0/bio/bedtools/genomecov"

rule download_chrom_sizes:
    output:
        "results/{species}/data/chromInfo.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/chromInfo.txt.gz -O {output}"

rule clean_chrom_sizes:
    input:
        "results/{species}/data/chromInfo.txt.gz"
    output:
        "results/{species}/data/chrom.sizes.txt"
    shell:
        "zcat {input} > {output}"
