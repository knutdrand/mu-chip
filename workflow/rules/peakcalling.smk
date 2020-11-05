genome_sizes = {"mm10": "mm", "hg38": "hs"}

rule macs2:
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
        "macs2 callpeak -t {input.treatment} -c {input.control} -g {params.gs} --bdg --broad --outdir results/{wildcards.species}/broadpeakcalling -n {wildcards.celltype}_{wildcards.condition}"


rule move_coverage:
    input:
        "results/{species}/broadpeakcalling/{combo}_treat_pileup.bdg"
    output:
        "results/{species}/coverage/{combo}.bdg"
    shell:
        "cp {input} {output}"

rule merge_domains:
    input:
        "results/{species}/broadpeakcalling/{combo}_peaks.broadPeak"
    output:
        "results/{species}/domains/{combo}.bed"
    shell:
        "bedtools merge -d 5000 -i {input} > {output}"

rule clip_bed:
    input:
        bdg="results/{species}/{folder}/{combo}{filetype}.{suffix}",
        sizes="results/{species}/data/chrom.sizes.txt"
    output:
        "results/{species}/{folder}/{combo}{filetype}.clipped.{suffix}"
    wildcard_constraints:
        filetype=".*",
        suffix="bed|bdg|narrowPeak|broadPeak"

    shell:
        "%s {input.bdg} | bedtools slop -i - -g {input.sizes} -b 0 | bedClip stdin {input.sizes} {output}" % chromosome_grep
