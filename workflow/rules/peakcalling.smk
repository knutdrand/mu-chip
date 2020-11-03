genome_sizes = {"mm10": "mm", "hg38": "hs"}

rule call_broad_peak:
    input:
        treatment="results/{species}/merged/{celltype}_{condition}.bed.gz",
        control="results/{species}/merged/{celltype}_input.bed.gz"
    output:
        "results/{species}/broadpeakcalling/{celltype}_{condition}_peaks.broadPeak",
        "results/{species}/broadpeakcalling/{celltype}_{condition}_treat_pileup.bdg",
        "results/{species}/broadpeakcalling/{celltype}_{condition}_control_lambda.bdg"
    conda:
        "envs/oldmacs.yaml"
    params:
        gs=lambda w: genome_sizes[w.species]
    shell:
        "macs2 -t {input.treatment} -c {input.control} -g {params.gs} --bdg --broad --outdir results/{wildcards.species}/broadpeakcalling -n {wildcards.celltype}_{wildcards.condition}"


rule move_coverage:
    input:
        "results/{species}/broadpeakcalling/{combo}_treat_pileup.bdg"
    output:
        "results/{species}/coverage/{combo}.bdg"
    shell:
        "cp {input} {output}"

rule merge_domains:
    input:
        "results/{species}/broadpeakcalling/{name}_peaks.broadPeak"
    output:
        "results/{species}/domains/{name}.bed"
    shell:
        "bedtools merge -d 5000 -i {input} > {output}"
