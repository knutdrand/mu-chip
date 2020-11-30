# rule cutadapt_pe:
#     input:
#         "results/reads/{pesample}_R1_001.fastq.gz",
#         "results/reads/{pesample}_R2_001.fastq.gz",
#     output:
#         fastq1=temp("results/trimmed/{pesample}_R1.fastq.gz"),
#         fastq2=temp("results/trimmed/{pesample}_R2.fastq.gz"),
#         qc="results/trimmed/{pesample}.qc.txt"
#     params:
#         adapters = '-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AATGATACGGCGACCACCGAGATCTACAC',
#         others = '--nextseq-trim=20 -m 10'
#     log:
#         "logs/cutadapt/{pesample}.log"
#     threads: 16
#     wrapper:
#         "0.66.0/bio/cutadapt/pe"

rule trimmomatic_pe:
    input:
        r1="results/reads/{pesample}_R1_001.fastq.gz",
        r2="results/reads/{pesample}_R2_001.fastq.gz"
    output:
        r1=temp("results/trimmed/{pesample}_R1.fastq.gz"),
        r2=temp("results/trimmed/{pesample}_R2.fastq.gz"),
        r1_unpaired=temp("results/trimmed/{pesample}_R1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/trimmed/{pesample}_R2.unpaired.fastq.gz"),
    log:
        "logs/trimmomatic/{pesample}.log"
    params:
        trimmer=["ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10", "SLIDINGWINDOW:4:20", "MINLEN:30"]
    threads:
        32
    wrapper:
        "0.67.0/bio/trimmomatic/pe"

rule bwa_mem_pe:
    input:
        reads=expand("results/trimmed/{{pesample}}_R{read}.fastq.gz", read=[1,2]),
        index="{params.index}"
    output:
        temp("results/{species}/bwa_mapped_pe/{pesample}.bam")
    log:
        "logs/bwa_mem/{species}/{pesample}.log"
    params:
        index=lambda w: config["index_path"].format(species=w.species),
        sort="samtools",
        sort_order="coordinate"
    threads:
        16
    wrapper:
        "0.49.0/bio/bwa/mem"

rule bowtie2_pe:
    input:
        sample=expand("results/trimmed/{{pesample}}_R{read}.fastq.gz", read=[1,2]),
        index=config["index_path"]+".1.bt2"
        # index="{params.index}.1.bt2"
    output:
        temp("results/{species}/bowtie2_unsorted_mapped_pe/{pesample}.bam")
    log:
        "logs/bowtie2/{species}/{pesample}.log"
    params:
        index=lambda w: config["index_path"].format(species=w.species),
        extra="-X 2000"  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "0.67.0/bio/bowtie2/align"

rule samtools_sort_pe:
    input:
        "results/{species}/bowtie2_unsorted_mapped_pe/{pesample}.bam"
    output:
        "results/{species}/bowtie2_mapped_pe/{pesample}.bam"
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.67.0/bio/samtools/sort"

rule remove_duplicates_pe:
    input:
        f"results/{{species}}/{mapper}_mapped_pe/{{pesample}}.bam"
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

rule remove_small_fragments_bam:
    input:
        "results/{species}/filtered_dedup_pe/{pesample}.bam"
    output:
        "results/{species}/size_filtered_dedup_pe/{pesample}.bam"        
    shell:
        """samtools view -h {input} | awk '($9>=0 ? $9 : -$9) > 50 || $1 ~ /^@/' | samtools view -bS - > {output}"""

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

