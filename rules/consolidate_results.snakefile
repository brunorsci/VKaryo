import os  # defensivo -- import idempotente, seguro mesmo se já importado alhures


def _consolidate_extra_flag(flag_name, config_key):
    """Monta a flag opcional só se o config apontar pro arquivo e ele existir
    -- painel de genes e tabela de verdade conhecida são específicos de uma coorte e não se aplicam a outra
    projeto, então ficam de fora silenciosamente quando ausentes."""
    path = config.get("consolidation", {}).get(config_key)
    if path and os.path.exists(path):
        return f"--{flag_name} {path}"
    return ""


rule consolidate_arriba:
    input:
        # arriba_annotated/ é nomeado pelo ID de RNA, não de DNA -- SAMPLES_BOTH
        # é a lista de IDs de DNA com par de RNA, então precisa mapear.
        expand(config["output"] + "output_vkaryo/arriba_annotated/{sample}_fusions_annotated.tsv",
               sample=[DNA_TO_RNA[d] for d in SAMPLES_BOTH]) if config.get("RNA") else [],
        driver_genes=config["variant_calling"]["driver_genes"],
        sample_mapping="sample_mapping.tsv",
    output:
        config["output"] + "output_vkaryo/consolidated/arriba_all_samples.tsv",
    params:
        extra=_consolidate_extra_flag("extra-gene-panel", "extra_gene_panel"),
        truth=_consolidate_extra_flag("ground-truth", "ground_truth"),
        outdir=config["output"] + "output_vkaryo",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/consolidate_program_results.py --program arriba \
            --output-dir {params.outdir} \
            --sample-mapping {input.sample_mapping} \
            --driver-genes {input.driver_genes} \
            {params.extra} {params.truth} \
            --out {output}
        """


rule consolidate_manta:
    input:
        expand(config["output"] + "output_vkaryo/annotsv/manta/{sample}_annotsv_manta.tsv", sample=SAMPLES) if config.get("DNA") else [],
        driver_genes=config["variant_calling"]["driver_genes"],
        sample_mapping="sample_mapping.tsv",
    output:
        config["output"] + "output_vkaryo/consolidated/manta_all_samples.tsv",
    params:
        extra=_consolidate_extra_flag("extra-gene-panel", "extra_gene_panel"),
        truth=_consolidate_extra_flag("ground-truth", "ground_truth"),
        outdir=config["output"] + "output_vkaryo",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/consolidate_program_results.py --program manta \
            --output-dir {params.outdir} \
            --sample-mapping {input.sample_mapping} \
            --driver-genes {input.driver_genes} \
            {params.extra} {params.truth} \
            --out {output}
        """


rule consolidate_delly:
    input:
        expand(config["output"] + "output_vkaryo/annotsv/delly/{sample}_annotsv_delly.tsv", sample=SAMPLES) if config.get("DNA") else [],
        driver_genes=config["variant_calling"]["driver_genes"],
        sample_mapping="sample_mapping.tsv",
    output:
        config["output"] + "output_vkaryo/consolidated/delly_all_samples.tsv",
    params:
        extra=_consolidate_extra_flag("extra-gene-panel", "extra_gene_panel"),
        truth=_consolidate_extra_flag("ground-truth", "ground_truth"),
        outdir=config["output"] + "output_vkaryo",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/consolidate_program_results.py --program delly \
            --output-dir {params.outdir} \
            --sample-mapping {input.sample_mapping} \
            --driver-genes {input.driver_genes} \
            {params.extra} {params.truth} \
            --out {output}
        """


rule consolidate_cnvkit:
    input:
        expand(config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.call.cns", sample=SAMPLES) if config.get("DNA") else [],
        driver_genes=config["variant_calling"]["driver_genes"],
        sample_mapping="sample_mapping.tsv",
    output:
        config["output"] + "output_vkaryo/consolidated/cnvkit_all_samples.tsv",
    params:
        extra=_consolidate_extra_flag("extra-gene-panel", "extra_gene_panel"),
        truth=_consolidate_extra_flag("ground-truth", "ground_truth"),
        outdir=config["output"] + "output_vkaryo",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/consolidate_program_results.py --program cnvkit \
            --output-dir {params.outdir} \
            --sample-mapping {input.sample_mapping} \
            --driver-genes {input.driver_genes} \
            {params.extra} {params.truth} \
            --out {output}
        """


rule consolidate_mutect2:
    input:
        expand(config["output"] + "output_vkaryo/mutect2/annovar/{sample}_mutect2.hg38_multianno.txt", sample=SAMPLES) if config.get("DNA") else [],
        driver_genes=config["variant_calling"]["driver_genes"],
        sample_mapping="sample_mapping.tsv",
    output:
        config["output"] + "output_vkaryo/consolidated/mutect2_all_samples.tsv",
    params:
        extra=_consolidate_extra_flag("extra-gene-panel", "extra_gene_panel"),
        truth=_consolidate_extra_flag("ground-truth", "ground_truth"),
        outdir=config["output"] + "output_vkaryo",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/consolidate_program_results.py --program mutect2 \
            --output-dir {params.outdir} \
            --sample-mapping {input.sample_mapping} \
            --driver-genes {input.driver_genes} \
            {params.extra} {params.truth} \
            --out {output}
        """
