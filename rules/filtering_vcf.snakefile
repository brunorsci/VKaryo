rule filter_leukemia_svs:
    input:
        vcf=config["output"] + "output_vkaryo/annotsv/{sample}_annotated.tsv"
    params:
        annotation_file="reference/related_genes.csv"  # Caminho do arquivo CSV com as anotações
    threads: 5
    output:
        filtered_vcf=config["output"] + "output_vkaryo/vcf_results/{sample}_related.tsv"
    shell:
        "python scripts/filtering_vcf.py -i {input.vcf} -o {output.filtered_vcf} -csv {params.annotation_file}"
