wildcard_constraints:
    pe_sra="|".join(samples.index[(samples["endedness"]=="pe") & (samples["sra"]==1)]) or "None",
    se_sra="|".join(samples.index[(samples["endedness"]=="se") & (samples["sra"]==1)]) or "None",
    sra="|".join(samples.index[samples["sra"]==1]) or "None"

rule prefetch:
    output:
        temp("results/sra_data/{sra}.sra")
    shell:
        "prefetch {wildcards.sra} -o {output}"

rule fastq_dump_pe:
    input:
        "results/sra_data/{pe_sra}.sra"
    output:
        "results/dumped_fastq/{pe_sra}_1.fastq.gz",
        "results/dumped_fastq/{pe_sra}_2.fastq.gz"
    shell:
        "fastq-dump --split-files --gzip {input} -O results/dumped_fastq/"

rule fastq_dump_se:
    input:
        "results/sra_data/{se_sra}.sra"
    output:
        "results/dumped_fastq/{se_sra}.fastq.gz"
    shell:
        "fastq-dump --gzip {input} -O results/dumped_fastq/"

rule rename_dumped_sra_pe:
    input:
        "results/dumped_fastq/{pe_sra}_{read}.fastq.gz"
    output:
        "results/reads/{pe_sra}_R{read}_001.fastq.gz"
    shell:
        "mv {input} {output}"

rule rename_dumped_sra_se:
    input:
        "results/dumped_fastq/{se_sra}.fastq.gz"
    output:
        "results/reads/{se_sra}_R1_001.fastq.gz"
    shell:
        "mv {input} {output}"
