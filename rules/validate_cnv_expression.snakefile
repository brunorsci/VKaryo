rule validate_cnv_expression:
    input:
        cnv = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns",
        # Use RNA sample ID for expression data (DNA_TO_RNA mapping)
        expression = lambda wildcards: config["output"] + f"output_vkaryo/star_aligned/{DNA_TO_RNA.get(wildcards.sample, wildcards.sample)}_ReadsPerGene.out.tab",
        gtf = config["refgenome"]["gtf"]
    output:
        validated = config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_validated.tsv",
        high_conf = config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_high_confidence.tsv"
    log:
        config["output"] + "logs/validate_cnv/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/validate_cnv_expression.py \
            --cnv {input.cnv} \
            --expression {input.expression} \
            --gtf {input.gtf} \
            --output {output.validated} \
            --output-high-conf {output.high_conf} \
            2>&1 | tee {log}
        """
