plottypes=["v", "heat", "tss", "average"]

wildcard_constraints:
    plottype="|".join(plottypes),

rule heatplot:
    input:
        bedgraph="results/{species}/broadpeakcalling/{combo}_treat_pileup.bdg",
        regions="results/{species}/domains/{combo}.bed"
    output:
        "results/plots/{species}/{combo}_heat.png",
        "results/plots/{species}/{combo}_heat.pkl"
    shell:
        "bdgplot heat {input.bedgraph} {input.regions} -o {output[0]} -od {output[1]}"

rule tssplot:
    input:
        bedgraph="results/{species}/broadpeakcalling/{combo}_treat_pileup.bdg",
        regions="results/{species}/data/unique_tss.bed"
    output:
        "results/plots/{species}/{combo}_tss.png",
        "results/plots/{species}/{combo}_tss.pkl"
    shell:
        "bdgplot tss {input.bedgraph} {input.regions} -o {output[0]} -od {output[1]}"

rule averageplot:
    input:
        bedgraph="results/{species}/broadpeakcalling/{combo}_treat_pileup.bdg",
        regions="results/{species}/domains/{combo}.bed"
    output:
        "results/plots/{species}/{combo}_average.png",
        "results/plots/{species}/{combo}_average.pkl"
    shell:
        "bdgplot average {input.bedgraph} {input.regions} -o {output[0]} -od {output[1]}"

rule joinplots:
    input:
        lambda w: expand("results/plots/{{species}}/{combo}_{{plottype}}.pkl", combo=get_combos(w.comparisongroup))
    output:
        report("results/plots/{species}/joined/{comparisongroup}_{plottype}.png", category="BDGPlots")
    shell:
        "bdgtools joinfigs {wildcards.plottype} {input} -o {output} --name {wildcards.comparisongroup}"
