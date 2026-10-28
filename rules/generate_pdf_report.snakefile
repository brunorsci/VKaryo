rule generate_pdf_report:
    input:
        classification = config["output"] + "output_vkaryo/classification/{sample}_classification.json",
        # Fusions use RNA sample ID (DNA_TO_RNA mapping)
        fusions = lambda wildcards: config["output"] + f"output_vkaryo/arriba_annotated/{DNA_TO_RNA.get(wildcards.sample, wildcards.sample)}_fusions_annotated.tsv",
        cnvs = config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_validated.tsv",
        ideogram = config["output"] + "output_vkaryo/plots/{sample}_aneuploidy_ideogram.pdf"
    output:
        report = config["output"] + "output_vkaryo/reports/{sample}_vkaryo_report.pdf"
    params:
        sample = "{sample}"
    log:
        config["output"] + "logs/pdf_report/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/generate_pdf_report.py \
            --sample {params.sample} \
            --classification {input.classification} \
            --fusions {input.fusions} \
            --cnvs {input.cnvs} \
            --ideogram {input.ideogram} \
            --output {output.report} \
            2>&1 | tee {log}
        """
