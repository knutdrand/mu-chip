from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")
if "endedness" not in samples:
    samples["endedness"] = config["endedness"]

wildcard_constraints:
    pesamples="|".join(samples.index[samples["endedness"]=="pe"])
