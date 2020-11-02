genome_sizes = {"mm10": "mm", "hg38": "hs"}

rule call_broad_peak:
    input:
        treatment="results/{species}/dedup/{celltype}_{condition}.bed.gz",
        control="results/{species}/dedup/{celltype}_input.bed.gz"
    output:
        "results/{species}/broadpeakcalling/{celltype}_{condition}_peaks.broadPeak",
        "results/{species}/broadpeakcalling/{celltype}_{condition}_treat_pileup.bdg"
    conda:
        "envs/oldmacs.yaml"
    shell:
        "macs2 -t {input.treatment} -c {input.control} -g {genome_sizes[wildcards.species]} --bdg --broad --outdir results/{snakemake.wildcards.species}/broadpeakcalling -n {snakemake.wildcards.celltype}_{snakemake.wildcards.condition}"

rule merge_domains:
    input:
        "results/{species}/broadpeakcalling/{name}_peaks.broadPeak"
    output:
        "results/{species}/domains/{name}.bed"
    shell:
        "bedtools merge -d 5000 -i {input} > {output}"
