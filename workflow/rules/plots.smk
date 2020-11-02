include: "rules/common.smk"


plottypes=["v", "heat", "tss", "average"]
wildcard_constraints:
    kind="treat_pileup|control_lambda|qvalues",
    plottype="|".join(plottypes),

rule heatplot:
    input:
        bedgraph="results/{species}/coverage/{celltype}_{condition}.bdg",
        regions="results/{species}/domains/{celltype}_{condition}.bed"
    output:
        "{species}/plots/{species}/{sample}_{kind}_average.png"
        "{species}/plots/{species}/{sample}_{kind}_average.pkl"
    shell:
        "bdgplot heat {input.bedgraph} {input.regions} -o {output[0]} -od {output[1]}"

rule tssplot:
    input:
        bedgraph="results/{species}/coverage/{celltype}_{condition}.bdg",
        regions="results/{species}/data/genes.bed"
    output:
        "{species}/plots/{species}/{sample}_{kind}_tss.png"
        "{species}/plots/{species}/{sample}_{kind}_tss.pkl"
    shell:
        "bdgplot tss {input.bedgraph} {input.regions} -o {output[0]} -od {output[1]}"

rule averageplot:
    input:
        bedgraph="results/{species}/coverage/{celltype}_{condition}.bdg",
        regions="results/{species}/domains/{celltype}_{condition}.bed"
    output:
        "{species}/plots/{species}/{sample}_{kind}_average.png"
        "{species}/plots/{species}/{sample}_{kind}_average.pkl"
    shell:
        "bdgplot average {input.bedgraph} {input.regions} -o {output[0]} -od {output[1]}"

rule joinplots:
    input:
        expand("results/plots/{species}/{combo}_{plottype}.pkl", combo=combos)
    output:
        report("results/plots/{species}/{comparisongroup}_{plottype}.png", category="BDGPlots")
    shell:
        "bdgtools joinfigs {wildcards.plottype} {input} -o {output} --name {wildcards.comparisongroup}"
