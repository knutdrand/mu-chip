from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"], sep="\t").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")
if "endedness" not in samples:
    samples["endedness"] = config.get("endedness", "se")

if "species" not in samples:
    samples["species"] = config.get("species", "mm10")

if "comparisongroup" not in samples:
    samples["comparisongroup"] = "all"


combos = {f"{celltype}_{condition}" for celltype, condition in zip(samples["celltype"], samples["condition"])}

def get_combos(comparisongroup):
    return {f"{celltype}_{condition}" for celltype, condition, g in zip(samples["celltype"], samples["condition"], samples["comparisongroup"]) if g == comparisongroup and condition.lower() != "input"}

def get_combos_for_species(species):
    return {f"{celltype}_{condition}" for celltype, condition, s in zip(samples["celltype"], samples["condition"], samples["species"]) if s == species and condition.lower() != "input"}

def expand_se_combo(format_string, wildcards):
    mask = (samples["endedness"]=="se")
    mask &= (samples["celltype"].str.lower()==wildcards.celltype.lower()) 
    mask &= (samples["condition"].str.lower()==wildcards.condition.lower())
    combos = list(samples.index[mask])
    return expand(format_string, species=wildcards.species, sample=combos)

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
