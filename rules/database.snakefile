# Conditional input functions for optional tools
def get_purecn_db_input(wildcards):
    """Return PureCN file if enabled, else empty list"""
    if config.get("tools", {}).get("purecn", False):
        return config["output"] + f"output_vkaryo/purecn/{wildcards.sample}_purity_ploidy.txt"
    return []

def get_optitype_db_input(wildcards):
    """Return OptiType file if enabled, else empty list"""
    if config.get("tools", {}).get("optitype", False):
        return config["output"] + f"output_vkaryo/optitype/{wildcards.sample}_hla_typing.txt"
    return []

def get_vardict_mrd_summary_input(wildcards):
    """Return VarDict MRD summary if enabled, else empty list"""
    if config.get("tools", {}).get("vardict", False):
        return config["output"] + f"output_vkaryo/vardict/{wildcards.sample}_mrd_summary.txt"
    return []

def get_vardict_mrd_tracking_input(wildcards):
    """
    Return VarDict MRD tracking if enabled, else empty list
    NOTE: vardict_track_mrd rule is disabled (requires ANNOVAR), so always return empty
    """
    # Always return empty list since vardict_track_mrd is disabled
    return []


rule populate_sample_database:
    """
    Popula banco de dados SQLite com todos os resultados da amostra
    """
    input:
        # PureCN (conditional)
        purecn_summary = get_purecn_db_input,
        # OptiType (conditional)
        hla_typing = get_optitype_db_input,
        # VarDict MRD (conditional)
        mrd_summary = get_vardict_mrd_summary_input,
        mrd_tracking = get_vardict_mrd_tracking_input,
        # Mutect2 driver mutations
        driver_mutations = config["output"] + "output_vkaryo/mutect2/{sample}_driver_mutations.tsv",
        # CNVs
        cnv_calls = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns",
        # All somatic SNVs (Mutect2 filtered - replaced Strelka)
        snvs = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_filtered.vcf.gz",
        # Consensus SNVs (Mutect2 + VarDict merged)
        consensus_snvs = config["output"] + "output_vkaryo/merged_variants/{sample}_driver_variants_final.tsv",
    output:
        db = config["output"] + "output_vkaryo/database/{sample}_results.db"
    log:
        config["output"] + "logs/database/{sample}.log"
    conda:
        "../envs/reporting.yaml"
    script:
        "../scripts/populate_database.py"


rule populate_integrated_database:
    """
    Popula database com resultados integrados DNA+RNA
    Adiciona CNVs validados, fusões, classificação
    """
    input:
        db = config["output"] + "output_vkaryo/database/{sample}_results.db",
        # CNV validado (DNA+RNA)
        cnv_validated = config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_validated.tsv",
        # Fusões anotadas (use RNA sample ID via DNA_TO_RNA mapping)
        fusions_annotated = lambda wildcards: config["output"] + f"output_vkaryo/arriba_annotated/{DNA_TO_RNA.get(wildcards.sample, wildcards.sample)}_fusions_annotated.tsv",
        # Classificação
        classification = config["output"] + "output_vkaryo/classification/{sample}_classification.json",
    output:
        db_integrated = config["output"] + "output_vkaryo/database/{sample}_results_integrated.db"
    log:
        config["output"] + "logs/database_integrated/{sample}.log"
    conda:
        "../envs/reporting.yaml"
    shell:
        """
        # Copia database base
        cp {input.db} {output.db_integrated}

        # Adiciona dados integrados ao database existente
        python scripts/populate_database.py \
            --db {output.db_integrated} \
            --sample {wildcards.sample} \
            --cnv-validated {input.cnv_validated} \
            --fusions {input.fusions_annotated} \
            --classification {input.classification} \
            2>&1 | tee {log}
        """


rule merge_all_databases:
    """
    Consolida databases de todas as amostras em um único database central
    """
    input:
        # Se tiver DNA+RNA integrado, usa integrado; senão usa básico
        dbs = lambda wildcards: (
            expand(config["output"] + "output_vkaryo/database/{sample}_results_integrated.db",
                   sample=list(set(SAMPLES) & set(SAMPLES_RNA)))
            if config["RNA"] and config["DNA"]
            else expand(config["output"] + "output_vkaryo/database/{sample}_results.db", sample=SAMPLES)
        )
    output:
        merged_db = config["output"] + "output_vkaryo/database/vkaryo_results.db"
    log:
        config["output"] + "logs/merge_databases/merge.log"
    conda:
        "../envs/reporting.yaml"
    script:
        "../scripts/merge_databases.py"
