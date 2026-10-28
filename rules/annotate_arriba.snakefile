rule annotate_arriba:
    input:
        fusions = config["output"] + "output_vkaryo/arriba/{sample}_fusions.tsv"
    output:
        annotated = config["output"] + "output_vkaryo/arriba_annotated/{sample}_fusions_annotated.tsv"
    log:
        config["output"] + "logs/annotate_arriba/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/annotate_arriba_leukemia.py \
            --input {input.fusions} \
            --output {output.annotated} \
            2>&1 | tee {log}
        """
