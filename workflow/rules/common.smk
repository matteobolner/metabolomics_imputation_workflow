import pandas as pd
from metabotk import MetaboTK


configfile: "config/config.yaml"


def get_mice_covariates():
    return ",".join(config["mice_imputation_covariates"])


def setup_dataset(file_path):
    dataset = MetaboTK().io.from_excel(
        file_path,
        sample_id_column=config["sample_id_column"],
        metabolite_id_column=config["metabolite_id_column"],
        sample_metadata_sheet=config["sample_metadata_sheet"],
        chemical_annotation_sheet=config["chemical_annotation_sheet"],
        data_sheet=config["data_sheet"],
    )
    return dataset


mice_seeds = config["mice_seeds"]
imputation_cycles = config["imputation_cycles"]


wildcard_constraints:
    mice_seed="[^_]+",
    imputation_cycle="[^_/]+",

