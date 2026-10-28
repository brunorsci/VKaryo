rule plot_aneuploidy:
    input:
        cnv = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns",
        validated = config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_validated.tsv"
    output:
        ideogram = config["output"] + "output_vkaryo/plots/{sample}_aneuploidy_ideogram.pdf",
        circular = config["output"] + "output_vkaryo/plots/{sample}_cnv_circular.pdf",
        summary = config["output"] + "output_vkaryo/plots/{sample}_aneuploidy_summary.txt"
    params:
        sample = "{sample}"
    log:
        config["output"] + "logs/plot_aneuploidy/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/plot_aneuploidy_ideogram.py \
            --cnv {input.cnv} \
            --validated {input.validated} \
            --sample {params.sample} \
            --output {output.ideogram} \
            --output-circular {output.circular} \
            --output-summary {output.summary} \
            2>&1 | tee {log}
        """
