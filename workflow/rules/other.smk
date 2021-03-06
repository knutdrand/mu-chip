# An example collection of Snakemake rules imported in the main Snakefile.

chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"


rule download_chrom_sizes:
    output:
        temp("results/{species}/data/chromInfo.txt.gz")
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/chromInfo.txt.gz -O {output}"

rule clean_chrom_sizes:
    input:
        "results/{species}/data/chromInfo.txt.gz"
    output:
        "results/{species}/data/chrom.sizes.txt"
    shell:
        "z%s {input} > {output}" % chromosome_grep

rule download_genes:
    output:
        "results/{species}/data/refGene.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/refGene.txt.gz -O {output}"

rule get_genes_bed:
    input:
        "results/{species}/data/refGene.txt.gz"
    output:
        "results/{species}/data/genes.bed"
    shell:
        """z%s {input} | awk '{{OFS="\t"}}{{print $3, $5, $6, ".", ".", $4}}' | uniq > {output}""" % chromosome_grep


rule get_unique_tss:
    input:
        "results/{species}/data/genes.bed"
    output:
        "results/{species}/data/unique_tss.bed"
    shell:
        """awk '{{OFS="\t"}}{{if ($6=="+") {{print $1,$2,$2+1,".",".","+"}} else {{print $1, $3-1, $3,".",".","-"}}}}' {input} | uniq > {output}"""

rule download_bowtie_index:
    output:
        temp(config["index_path"]+".zip")
    log:
        "logs/wget/{species}.log"
    shell:
        "wget https://genome-idx.s3.amazonaws.com/bt/{wildcards.species}.zip -O {output} -o {log}"

rule unzip_bowtie_index:
    input:
        "{path}/{species}.zip"
    output:
        "{path}/{species}.1.bt2"
    shell:
        "unzip {input} -d {wildcards.path}"


rule fastqc:
    input:
        "results/{folder}/{filename}.fastq.gz"
    output:
        html="results/qc/fastqc/{folder}/{filename}.html",
        zip="results/qc/fastqc/{folder}/{filename}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{folder}/{filename}.log"
    threads: 1
    wrapper:
        "0.67.0/bio/fastqc"
