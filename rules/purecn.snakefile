rule purecn_coverage:
    """
    Calcula cobertura para PureCN
    """
    input:
        bam = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        intervals = config["refgenome"]["exome_bed"]
    output:
        coverage = config["output"] + "output_vkaryo/purecn/{sample}_coverage_loess.txt.gz"
    log:
        config["output"] + "logs/purecn_coverage/{sample}.log"
    threads: 4
    conda:
        "../envs/purecn.yaml"
    shell:
        """
        # Find Coverage.R script from PureCN installation
        PURECN_EXTDATA=$(Rscript -e "cat(system.file('extdata', package='PureCN'))")

        Rscript "$PURECN_EXTDATA/Coverage.R" \
            --out-dir $(dirname {output.coverage}) \
            --bam {input.bam} \
            --intervals {input.intervals} \
            --cores {threads} \
            2>&1 | tee {log}
        """


rule purecn_run:
    """
    Executa PureCN para estimar pureza, ploidy e CNVs corrigidos
    SEM controle normal (usa PoN)
    Usa VCF do Mutect2 (somatic filtered)
    """
    input:
        coverage = config["output"] + "output_vkaryo/purecn/{sample}_coverage_loess.txt.gz",
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_filtered.vcf.gz",
        intervals = config["refgenome"]["exome_bed"]
    output:
        rds = config["output"] + "output_vkaryo/purecn/{sample}.rds",
        csv = config["output"] + "output_vkaryo/purecn/{sample}.csv",
        pdf = config["output"] + "output_vkaryo/purecn/{sample}.pdf",
        genes = config["output"] + "output_vkaryo/purecn/{sample}_genes.csv",
        loh = config["output"] + "output_vkaryo/purecn/{sample}_loh.csv"
    params:
        out_dir = config["output"] + "output_vkaryo/purecn/",
        sample_id = "{sample}",
        genome = "hg38"
    log:
        config["output"] + "logs/purecn/{sample}.log"
    threads: 8
    conda:
        "../envs/purecn.yaml"
    shell:
        """
        # Find PureCN.R script from PureCN installation
        PURECN_EXTDATA=$(Rscript -e "cat(system.file('extdata', package='PureCN'))")

        Rscript "$PURECN_EXTDATA/PureCN.R" \
            --out {params.out_dir} \
            --tumor {input.coverage} \
            --sampleid {params.sample_id} \
            --vcf {input.vcf} \
            --intervals {input.intervals} \
            --genome {params.genome} \
            --fun-segmentation PSCBS \
            --post-optimize \
            --seed 123 \
            2>&1 | tee {log}
        """


rule purecn_summary:
    """
    Extrai sumário de pureza/ploidy do PureCN
    """
    input:
        rds = config["output"] + "output_vkaryo/purecn/{sample}.rds"
    output:
        summary = config["output"] + "output_vkaryo/purecn/{sample}_purity_ploidy.txt"
    log:
        config["output"] + "logs/purecn_summary/{sample}.log"
    conda:
        "../envs/purecn.yaml"
    script:
        "../scripts/extract_purecn_summary.R"
