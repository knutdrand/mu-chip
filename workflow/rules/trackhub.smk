track_types = ["domains.bb", "treat_pileup.bw", "control_lambda.bw"]

rule create_trackhub:
    input:
        lambda wildcards: expand("results/trackhub/{{species}}/{combo}_{filetype}",
                                 combo=get_combos_for_species(wildcards.species),
                                 filetype=track_types)
    output:
        "results/trackhub/{species}/trackDb.txt"
    script:
        "../scripts/trackhub.py"

rule ucsc_sort:
    input:
        "results/{species}/broadpeakcalling/{combo}{filetype}.clipped.bdg"
    output:
        bdg="results/{species}/broadpeakcalling/{combo}{filetype}.clipped.bdg.uscssort",
        tmp=directory("results/tmp/{species}-{combo}-{filetype}/")
    wildcard_constraints:
        filetype=".*"
    shell:
        "LC_COLLATE=C sort -k1,1 -k2,2n {input} -T {output.tmp} > {output.bdg}"

rule create_bw_track:
    input:
        bedGraph="results/{species}/broadpeakcalling/{combo}{filetype}.clipped.bdg.uscssort",
        chromsizes="results/{species}/data/chrom.sizes.txt"
    output:
        "results/trackhub/{species}/{combo}{filetype}.bw"
    wrapper:
        "0.50.3/bio/ucsc/bedGraphToBigWig"

rule trunctate_score:
    input:
        "results/{species}/broadpeakcalling/{combo}_{filetype}.clipped.broadPeak"
    output:
        temp("results/{species}/broadpeakcalling/{combo}_{filetype}.clipped.broadPeak.trunc")
    shell:
        """awk  '{{OFS="\t"}}{{$5=($5>1000?1000:$5);print}} {{input}} > {{output}}'"""

rule create_peak_track:
    input:
        peaks="results/{species}/broadpeakcalling/{combo}_peaks.clipped.broadPeak.trunc",
        sizes="results/{species}/data/chrom.sizes.txt"
    output:
        "results/trackhub/{species}/{combo}_peaks.bb"
    shell:
        "bedToBigBed -type=bed6+3 {input.peaks} {input.sizes} {output}"

rule create_domain_track:
    input:
        domains="results/{species}/domains/{combo}.clipped.bed",
        sizes="results/{species}/data/chrom.sizes.txt"
    output:
        "results/trackhub/{species}/{combo}_domains.bb"
    shell:
        "bedToBigBed {input.domains} {input.sizes} {output}"
