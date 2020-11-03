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
        """awk '{{OFS="\t"}}{{if ($6=="+") {{print $1,$2,$2+1}} else {{print $1, $3-1, $3}}}}' {input} | uniq > {output}"""

