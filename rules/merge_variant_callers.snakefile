rule merge_mutect2_vardict:
    """
    Merge variants from Mutect2 (high specificity) and VarDict (high sensitivity)
    Strategy:
    - Use Mutect2 for all genes (high confidence, filtered)
    - Add VarDict low-VAF variants for driver genes (FLT3, NPM1, IDH1/2, etc.)
    - Deduplicate and prioritize by evidence
    """
    input:
        mutect2 = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_pass.vcf.gz",
        vardict = config["output"] + "output_vkaryo/vardict/{sample}_mrd_filtered.vcf.gz",
        drivers = config["variant_calling"]["driver_genes"]
    output:
        merged = config["output"] + "output_vkaryo/merged_variants/{sample}_merged_snvs.vcf.gz",
        merged_tbi = config["output"] + "output_vkaryo/merged_variants/{sample}_merged_snvs.vcf.gz.tbi",
        stats = config["output"] + "output_vkaryo/merged_variants/{sample}_merge_stats.txt"
    log:
        config["output"] + "logs/merge_variants/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        # Merge VCFs using bcftools
        # VarDict variants get priority for driver genes, Mutect2 for others
        bcftools merge \
            --force-samples \
            -o {output.merged} \
            -O z \
            {input.mutect2} {input.vardict} \
            2>&1 | tee {log}

        tabix -f -p vcf {output.merged}

        # Generate merge statistics
        echo "Merge Statistics for {wildcards.sample}" > {output.stats}
        echo "=======================================" >> {output.stats}
        echo "" >> {output.stats}

        echo "Mutect2 variants:" >> {output.stats}
        bcftools view -H {input.mutect2} | wc -l >> {output.stats}

        echo "VarDict variants:" >> {output.stats}
        bcftools view -H {input.vardict} | wc -l >> {output.stats}

        echo "Merged total:" >> {output.stats}
        bcftools view -H {output.merged} | wc -l >> {output.stats}

        echo "" >> {output.stats}
        echo "Merge completed: $(date)" >> {output.stats}
        """


rule prioritize_driver_variants:
    """
    Prioritize and annotate driver gene variants
    Combines high-confidence Mutect2 calls with sensitive VarDict calls
    """
    input:
        merged = config["output"] + "output_vkaryo/merged_variants/{sample}_merged_snvs.vcf.gz",
        drivers = config["variant_calling"]["driver_genes"]
    output:
        drivers_vcf = config["output"] + "output_vkaryo/merged_variants/{sample}_driver_variants_final.vcf.gz",
        drivers_tbi = config["output"] + "output_vkaryo/merged_variants/{sample}_driver_variants_final.vcf.gz.tbi",
        drivers_tsv = config["output"] + "output_vkaryo/merged_variants/{sample}_driver_variants_final.tsv"
    params:
        min_vaf = config["variant_calling"]["min_vaf"]
    log:
        config["output"] + "logs/merge_variants/{sample}_prioritize.log"
    conda:
        # gatk.yaml doesn't declare pysam/pandas (script is pure Python, needs neither
        # GATK nor Java) -- python_analysis.yaml already pins a working pysam=0.21
        "../envs/python_analysis.yaml"
    script:
        "../scripts/prioritize_drivers.py"


rule consensus_variant_calling:
    """
    Create consensus variant call set
    HIGH confidence = detected by both Mutect2 AND VarDict
    MEDIUM confidence = detected by one caller in driver gene
    LOW confidence = detected by one caller in non-driver gene
    """
    input:
        merged = config["output"] + "output_vkaryo/merged_variants/{sample}_merged_snvs.vcf.gz",
        mutect2 = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_pass.vcf.gz",
        vardict = config["output"] + "output_vkaryo/vardict/{sample}_mrd_filtered.vcf.gz",
        drivers = config["variant_calling"]["driver_genes"]
    output:
        high_conf = config["output"] + "output_vkaryo/merged_variants/{sample}_high_confidence.vcf.gz",
        high_conf_tbi = config["output"] + "output_vkaryo/merged_variants/{sample}_high_confidence.vcf.gz.tbi",
        medium_conf = config["output"] + "output_vkaryo/merged_variants/{sample}_medium_confidence.vcf.gz",
        medium_conf_tbi = config["output"] + "output_vkaryo/merged_variants/{sample}_medium_confidence.vcf.gz.tbi",
        summary = config["output"] + "output_vkaryo/merged_variants/{sample}_consensus_summary.txt"
    log:
        config["output"] + "logs/merge_variants/{sample}_consensus.log"
    conda:
        "../envs/python_analysis.yaml"
    script:
        "../scripts/consensus_calling.py"
