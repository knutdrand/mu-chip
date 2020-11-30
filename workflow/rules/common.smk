from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"
defaults = {"species": "mm10", 
            "endedness": "pe", 
            "sra": False, 
            "comparisongroup": "all"}

validate(config, schema="../schemas/config.schema.yaml")
chromosome_grep = "grep -Ew -e 'chr[0-9]{{1,2}}' -e chrX -e chrY"
samples = pd.read_table(config["samples"], sep="\t").set_index("sample", drop=False)
mapper = config.get("mapper", "bowtie2")
validate(samples, schema="../schemas/samples.schema.yaml")
for key, value in defaults.items():
    if key not in samples:
        samples[key] = config.get(key, value)
combos = {f"{celltype}_{condition}" for celltype, condition in zip(samples["celltype"], samples["condition"])}
combo_frame = samples.groupby(by=["species", "celltype", "condition"], as_index=False).last()
def get_combos(comparisongroup):
    return {f"{celltype}_{condition}" for celltype, condition, g in zip(samples["celltype"], samples["condition"], samples["comparisongroup"]) if g == comparisongroup and condition.lower() != "input"}

def expand_species(format_string, species):
    mask = combo_frame["species"] == species
    mask &= combo_frame["condition"].str.lower() != "input"
    res =  [format_string.format(endedness=row["endedness"], 
                                 celltype=row["celltype"],
                                 condition=row["condition"])
            for _, row in combo_frame[mask].iterrows()]
    return res

def expand_comparison_group(format_string, comparisongroup):
    mask = combo_frame["comparisongroup"]==comparisongroup
    mask &= combo_frame["condition"].str.lower() != "input"
    return [format_string.format(species=row["species"],
                                 endedness=row["endedness"],
                                 celltype=row["celltype"],
                                 condition=row["condition"])
            for _, row in combo_frame[mask].iterrows()]

def expand_se_combo(format_string, wildcards):
    mask = (samples["endedness"]=="se")
    mask &= (samples["celltype"].str.lower()==wildcards.celltype.lower()) 
    mask &= (samples["condition"].str.lower()==wildcards.condition.lower())
    mask &= (samples["species"].str.lower()==wildcards.species.lower())
    combos = list(samples.index[mask])
    return expand(format_string, species=wildcards.species, sample=combos)

def expand_pe_combo(format_string, wildcards):
    mask = (samples["endedness"]=="pe")
    mask &= (samples["celltype"].str.lower()==wildcards.celltype.lower()) 
    mask &= (samples["condition"].str.lower()==wildcards.condition.lower())
    mask &= (samples["species"].str.lower()==wildcards.species.lower())
    combos = list(samples.index[mask])
    return expand(format_string, species=wildcards.species, pesample=combos)

def expand_pe_input(format_string, wildcards):
    mask = (samples["endedness"]=="pe")
    mask &= (samples["celltype"].str.lower()==wildcards.celltype.lower())
    mask &= (samples["condition"].str.lower()=="input")
    mask &= (samples["species"].str.lower()==wildcards.species.lower())
    combos = list(samples.index[mask])
    return expand(format_string, species=wildcards.species, pesample=combos)

def expand_all_combos(format_string, wildcards):
    mask = (samples["endedness"]=="se")
    mask &= (samples["celltype"].str.lower()==wildcards.celltype.lower()) 
    mask &= (samples["condition"].str.lower()==wildcards.condition.lower())
    combos = list(samples.index[mask])
    return expand(format_string, species=wildcards.species, sample=combos)

wildcard_constraints:
    pesample="|".join(samples.index[samples["endedness"]=="pe"]),
    sample="|".join(samples.index[samples["endedness"]=="se"]),
    combo="|".join(combos),
    species="|".join(set(samples["species"]))
