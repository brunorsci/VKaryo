"""
Snakemake rule para consolidar fusões gênicas de múltiplos callers
"""

rule consolidate_fusions:
    input:
        genefuse = config["output"] + "output_vkaryo/genefuse/{sample}_fusions.json",
        arriba = config["output"] + "output_vkaryo/arriba/{sample}_fusions.tsv",
        manta = config["output"] + "output_vkaryo/manta_rna/{sample}/results/variants/candidateSV.vcf.gz"
    output:
        consolidated = config["output"] + "output_vkaryo/consolidated_fusions/{sample}_consolidated_fusions.tsv",
        high_conf = config["output"] + "output_vkaryo/consolidated_fusions/{sample}_high_confidence_fusions.tsv"
    params:
        min_support = 1  # Mínimo de callers suportando a fusão
    log:
        config["output"] + "logs/consolidate_fusions/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    script:
        "../scripts/consolidate_fusions.py --genefuse {input.genefuse} --arriba {input.arriba} --manta {input.manta} --output {output.consolidated} --min-support {params.min_support} 2>&1 | tee {log}"


rule filter_high_confidence_fusions:
    """
    Filtra apenas fusões de alta confiança (suportadas por múltiplos callers)
    """
    input:
        config["output"] + "output_vkaryo/consolidated_fusions/{sample}_consolidated_fusions.tsv"
    output:
        config["output"] + "output_vkaryo/consolidated_fusions/{sample}_high_confidence_fusions.tsv"
    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep='\t')

        # Filtra alta confiança
        high_conf = df[
            (df['support_count'] >= 2) |  # Múltiplos callers
            (df['consolidated_confidence'].isin(['HIGH', 'VERY_HIGH']))
        ]

        high_conf.to_csv(output[0], sep='\t', index=False)
