rule classify_leukemia:
    input:
        # Fusions use RNA sample ID (DNA_TO_RNA mapping)
        fusions = lambda wildcards: config["output"] + f"output_vkaryo/arriba_annotated/{DNA_TO_RNA.get(wildcards.sample, wildcards.sample)}_fusions_annotated.tsv",
        cnvs = config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_validated.tsv",
        # NOTA: usa o multianno do ANNOVAR, não merged_variants/*_driver_variants_final.tsv
        # (Mutect2+VarDict sem filtro de frequência populacional -- maioria
        # dos "achados" seria polimorfismo germinativo comum). O multianno
        # já tem gnomAD AF e ClinVar para filtragem adequada.
        snvs = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_mutect2.hg38_multianno.txt",
        aneuploidy = config["output"] + "output_vkaryo/plots/{sample}_aneuploidy_summary.txt",
        drivers = config["variant_calling"]["driver_genes"]
    output:
        classification = config["output"] + "output_vkaryo/classification/{sample}_classification.json",
        report = config["output"] + "output_vkaryo/classification/{sample}_classification_report.md"
    params:
        sample = "{sample}"
    log:
        config["output"] + "logs/classify_leukemia/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/classify_leukemia.py \
            --fusions {input.fusions} \
            --cnvs {input.cnvs} \
            --snvs {input.snvs} \
            --driver-genes {input.drivers} \
            --aneuploidy-summary {input.aneuploidy} \
            --sample {params.sample} \
            --output {output.classification} \
            --report {output.report} \
            2>&1 | tee {log}
        """
