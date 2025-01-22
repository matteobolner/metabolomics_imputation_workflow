from metabotk import MetaboTK


ds = MetaboTK().io.from_excel(
    file_path=snakemake.input.dataset,
    sample_id_column=snakemake.config["sample_id_column"],
    metabolite_id_column=snakemake.config["metabolite_id_column"],
)


no_empty = ds.stats.remove_missing(on="metabolites", threshold=0.99999999)
ds = ds.ops.subset(what="metabolites", ids=list(no_empty.columns))
ds.data = ds.stats.remove_outliers(threshold=5)

low_missing = ds.stats.remove_missing(threshold=0.25)

removed = ds.ops.drop(
    what="metabolites", ids=list(low_missing.columns)
).chemical_annotation
removed.to_csv(snakemake.output.missing_removed, index=True, sep="\t")
ds = ds.ops.subset(what="metabolites", ids=list(low_missing.columns))
ds.io.save_excel(snakemake.output.dataset)
