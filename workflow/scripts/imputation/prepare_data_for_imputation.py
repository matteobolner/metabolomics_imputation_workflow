from metabotk import MetaboTK

ds = MetaboTK().io.from_excel(
    file_path=snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)
ds.ops.merge_sample_metadata_data().to_csv(
    snakemake.output.data_metadata, sep="\t", index=True
)
ds.chemical_annotation.to_csv(
    snakemake.output.chemical_annotation, sep="\t", index=True
)
ds.sample_metadata.to_csv(snakemake.output.sample_metadata, sep="\t", index=True)
