# An example collection of Snakemake rules imported in the main Snakefile.

chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"

rule download_genes:
    output:
        "{species}/data/refGene.txt.gz"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.species}/database/refGene.txt.gz -O {output}"

rule get_genes_bed:
    input:
        lambda wildcards: "{species}/data/refGene.txt.gz"
    output:
        "{species}/data/genes.bed"
    shell:
        """z%s {input} | awk '{{OFS="\t"}}{{print $3, $5, $6, ".", ".", $4}}' | uniq > {output}""" % chromosome_grep
