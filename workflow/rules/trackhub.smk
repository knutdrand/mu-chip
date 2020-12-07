track_types = ["domains.bb", "treat_pileup.bw", "control_lambda.bw"]

rule create_genomes_file:
    input:
        expand("results/trackhub/{species}/trackDb.txt", species=set(samples["species"]))
    output:
        "results/trackhub/genomes.txt"
    run:
        text = "\n".join(
            f"genome {species}\ntrackDb {species}/trackDb.txt\n" 
            for species in set(samples["species"]))
        open(output[0], "w").write(text)

rule create_trackhub:
    input:
        "results/trackhub/genomes.txt"
    output:
        "results/trackhub/hub.txt"
    run:
        name = config.get("name", "mu-chip")
        mail = config.get("mail", "example@mail.com")
        open(output[0], "w").write(f"""hub {name}
shortLabel {name}
longLabel {name}
genomesFile genomes.txt
email {mail}
descriptionUrl ucscHub.html
""")

rule create_trackdb:
    input:
        domains = lambda w: expand_species("results/trackhub/{{species}}/{endedness}_{celltype}_{condition}_domains.bb", species=w.species),
        all_tracks = lambda w: [fstr % tt for tt in track_types for fstr in expand_species("results/trackhub/{{species}}/{endedness}_{celltype}_{condition}_%s", species=w.species)]
    output:
        "results/trackhub/{species}/trackDb.txt"
    conda:
        "../envs/trackhub.yaml"
    script:
        "../scripts/trackhub.py"

rule ucsc_sort:
    input:
        "results/{species}/{endedness}_broadpeakcalling/{combo}{filetype}.clipped.bdg"
    output:
        bdg="results/{species}/{endedness}_broadpeakcalling/{combo}{filetype}.clipped.bdg.uscssort",
    params:
        tmp="results/tmp/ucsc_sort_{endedness}-{species}-{combo}-{filetype}"
    wildcard_constraints:
        filetype=".*",
    shell:
        "LC_COLLATE=C sort -k1,1 -k2,2n {input} -T {params.tmp} > {output.bdg}"

rule create_bw_track:
    input:
        bedGraph="results/{species}/{endedness}_broadpeakcalling/{combo}{filetype}.clipped.bdg.uscssort",
        chromsizes="results/{species}/data/chrom.sizes.txt"
    output:
        "results/trackhub/{species}/{endedness}_{combo}{filetype}.bw"
    wrapper:
        "0.50.3/bio/ucsc/bedGraphToBigWig"

rule trunctate_score:
    input:
        "results/{species}{endedness}_/broadpeakcalling/{combo}_{filetype}.clipped.broadPeak"
    output:
        temp("results/{species}/{endedness}_broadpeakcalling/{combo}_{filetype}.clipped.broadPeak.trunc")
    shell:
        """awk  '{{OFS="\t"}}{{$5=($5>1000?1000:$5);print}} {{input}} > {{output}}'"""

rule create_peak_track:
    input:
        peaks="results/{species}/{endedness}_broadpeakcalling/{combo}_peaks.clipped.broadPeak.trunc",
        sizes="results/{species}/data/chrom.sizes.txt"
    output:
        "results/trackhub/{species}/{endedness}_{combo}_peaks.bb"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "bedToBigBed -type=bed6+3 {input.peaks} {input.sizes} {output}"

rule create_domain_track:
    input:
        domains="results/{species}/{endedness}_domains/{combo}.clipped.bed",
        sizes="results/{species}/data/chrom.sizes.txt"
    output:
        "results/trackhub/{species}/{endedness}_{combo}_domains.bb"
    conda:
        "../envs/ucsc.yaml"
    shell:
        "bedToBigBed {input.domains} {input.sizes} {output}"
