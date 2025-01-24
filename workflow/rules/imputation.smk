rule split_subsets_for_imputation:
    input:
        dataset=config["raw_dataset"],
    output:
        group="data/imputation/{group}/input/raw_dataset.xlsx",
    run:
        dataset = setup_dataset(input.dataset)
        if config["impute_groups_separately"] == False:
            dataset.io.save_excel(output.group)
        else:
            dataset_subset = dataset.ops.split(
                by="samples", columns=config["imputation_group_column"]
            )[wildcards.group]
            print(dataset_subset.samples)
            dataset_subset.io.save_excel(output.group)


rule remove_outliers_and_missing:
    input:
        dataset=rules.split_subsets_for_imputation.output.group,
    output:
        missing_removed="data/imputation/{group}/missing_removed.tsv",
        dataset="data/imputation/{group}/input/clean_raw_dataset.xlsx",
    script:
        "../scripts/remove_outliers_and_missing.py"


rule prepare_dataset_for_imputation:
    input:
        dataset=rules.remove_outliers_and_missing.output.dataset,
    output:
        data_metadata="data/imputation/{group}/input/data_metadata.tsv",
        chemical_annotation="data/imputation/{group}/input/chemical_annotation.tsv",
        sample_metadata="data/imputation/{group}/input/samples.tsv",
    script:
        "../scripts/prepare_data_for_imputation.py"


rule download_imputation_script:
    output:
        script="workflow/scripts/impute.R",
    shell:
        "curl https://raw.githubusercontent.com/matteobolner/metabolomics_missing_data_imputation/refs/heads/main/impute.R -o {output.script}"


rule impute:
    input:
        data=rules.prepare_dataset_for_imputation.output.data_metadata,
        chemical_annotation=rules.prepare_dataset_for_imputation.output.chemical_annotation,
        script=rules.download_imputation_script.output.script,
    output:
        imputed="data/imputation/{group}/imputed/{mice_seed}.tsv",
    params:
        covariates=get_mice_covariates(),
        metabolite_id_column=config["metabolite_id_column"],
        super_pathway_column=config["super_pathway_column"],
    conda:
        "../envs/imputation.yaml"
    shell:
        """Rscript {input.script} -d {input.data} -c {input.chemical_annotation} -o {output.imputed} -s {wildcards.mice_seed} -m pmm -n 5 -r 0.25 -u 5 -a {params.covariates} --metabolite_id_column "{params.metabolite_id_column}" --super_pathway_column "{params.super_pathway_column}" """


rule split_imputed_datasets:
    input:
        imputed_data=rules.impute.output.imputed,
        samples=rules.prepare_dataset_for_imputation.output.sample_metadata,
        chemical_annotation=rules.prepare_dataset_for_imputation.output.chemical_annotation,
    output:
        imputations=expand(
            "data/imputation/{{group}}/imputed/seed_{{mice_seed}}/imputation_{imputation_cycle}.xlsx",
            imputation_cycle=[i for i in imputation_cycles],
        ),
    run:
        df = pd.read_table(input.imputed_data)
        df = df[df[".imp"] != 0]

        outputs = {
            i + 1: output.imputations[i] for i in range(0, len(config["mice_seeds"]))
        }

        for name, group in df.groupby(by=".imp"):
            group = group.drop(columns=[".id", ".imp"])
            dataset = MetaboTK().io.from_tables(
                sample_metadata=input.samples,
                chemical_annotation=input.chemical_annotation,
                data=group,
                sample_id_column=config["sample_id_column"],
                metabolite_id_column=config["metabolite_id_column"],
            )
            dataset.io.save_excel(outputs[name])


if config["impute_groups_separately"]:

    rule get_imputations:
        input:
            whole_raw_dataset=config["raw_dataset"],
            groups=expand(
                "data/imputation/{group}/imputed/seed_{{mice_seed}}/imputation_{{imputation_cycle}}.xlsx",
                group=config["imputation_groups"],
            ),
        output:
            dataset="data/imputation/imputed/seed_{mice_seed}/imputation_{imputation_cycle}.xlsx",
        run:
            whole_raw_dataset = setup_dataset(input.whole_raw_dataset)
            merged_data = []
            merged_sample_annotation = []
            merged_chemical_annotation = whole_raw_dataset.chemical_annotation
            for i in input.groups:
                tempsubset = setup_dataset(i)
                merged_data.append(tempsubset.data)
                merged_sample_annotation.append(tempsubset.sample_metadata)
            merged_data = pd.concat(merged_data)
            merged_sample_metadata = pd.concat(merged_sample_annotation)
            merged_dataset = MetaboTK().io.from_tables(
                data=merged_data,
                sample_metadata=merged_sample_metadata,
                chemical_annotation=merged_chemical_annotation,
                sample_id_column=whole_raw_dataset._sample_id_column,
                metabolite_id_column=whole_raw_dataset._metabolite_id_column,
            )
            merged_dataset.io.save_excel(output.dataset)


else:

    rule get_imputations:
        input:
            dataset=expand(
                "data/imputation/{group}/imputed/seed_{{mice_seed}}/imputation_{{imputation_cycle}}.xlsx",
                group="whole",
            ),
        output:
            dataset="data/imputation/imputed/seed_{mice_seed}/imputation_{imputation_cycle}.xlsx",
        shell:
            "cp {input.dataset} {output.dataset}"
