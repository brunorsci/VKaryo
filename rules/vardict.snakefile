rule vardict_call:
    """
    VarDict para detecção de variantes de baixa frequência (MRD)
    Modo tumor-only com filtros rigorosos
    Usa BED de genes drivers (75 genes, ~2Mb total vs 3Gb do genoma)
    """
    input:
        bam = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        bai = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bai",
        bed = config["variant_calling"]["driver_genes_bed"],  # BED de genes drivers
        ref = config["refgenome"]["ref"]
    output:
        vcf = config["output"] + "output_vkaryo/vardict/{sample}_vardict.vcf.gz",
        tbi = config["output"] + "output_vkaryo/vardict/{sample}_vardict.vcf.gz.tbi"
    params:
        sample = "{sample}",
        # Parâmetros para detecção sensível de baixa VAF
        min_af = 0.01,  # VAF mínimo 1% para MRD
        min_qual = 20
    log:
        config["output"] + "logs/vardict/{sample}.log"
    threads: 8
    conda:
        "../envs/vardict.yaml"
    shell:
        """
        # Fix locale issue: VarDict expects decimal point not comma
        export LC_NUMERIC=C
        export LC_ALL=C

        # VarDict tumor-only mode
        # BED já tem 4 colunas: chr, start, end, gene_name
        (vardict-java \
            -G {input.ref} \
            -f {params.min_af} \
            -N {params.sample} \
            -b {input.bam} \
            -c 1 -S 2 -E 3 -g 4 \
            -Q {params.min_qual} \
            -th {threads} \
            {input.bed} \
        | teststrandbias.R \
        | var2vcf_valid.pl -N {params.sample} -E -f {params.min_af} \
        | bgzip -c > {output.vcf}) 2>&1 | tee {log}

        tabix -f -p vcf {output.vcf}
        """


rule vardict_filter_mrd:
    """
    Filtra variantes VarDict para MRD clínico
    Foco em genes recorrentes em leucemia
    """
    input:
        vcf = config["output"] + "output_vkaryo/vardict/{sample}_vardict.vcf.gz"
    output:
        filtered = config["output"] + "output_vkaryo/vardict/{sample}_mrd_filtered.vcf.gz",
        mrd_report = config["output"] + "output_vkaryo/vardict/{sample}_mrd_summary.txt"
    params:
        # Genes recorrentes em leucemia para MRD
        mrd_genes = "FLT3,NPM1,DNMT3A,TET2,IDH1,IDH2,TP53,RUNX1,CEBPA,KIT,NRAS,KRAS,PTPN11,JAK2,ASXL1,EZH2,CBL,SF3B1,SRSF2,U2AF1,ZRSR2",
        min_vaf_mrd = 0.005,  # 0.5% para MRD profundo
        max_vaf_germline = 0.35  # Evita germline (VAF ~50%)
    log:
        config["output"] + "logs/vardict_filter/{sample}.log"
    conda:
        "../envs/vardict.yaml"
    script:
        "../scripts/filter_vardict_mrd.py"


rule vardict_annovar:
    """
    Anota VarDict com ANNOVAR para contexto clínico (instalação local)
    """
    input:
        vcf = config["output"] + "output_vkaryo/vardict/{sample}_mrd_filtered.vcf.gz"
    output:
        annotated = config["output"] + "output_vkaryo/vardict/annovar/{sample}_vardict.hg38_multianno.txt",
        vcf_annotated = config["output"] + "output_vkaryo/vardict/annovar/{sample}_vardict_annotated.vcf"
    params:
        out_prefix = config["output"] + "output_vkaryo/vardict/annovar/{sample}_vardict",
        annovar_dir = config["refgenome"]["annovar_dir"],
        annovar_db = config["refgenome"]["annovar_db"],
        buildver = "hg38"
    log:
        config["output"] + "logs/vardict_annovar/{sample}.log"
    threads: 4
    conda:
        "../envs/annovar.yaml"
    shell:
        """
        {params.annovar_dir}/table_annovar.pl {input.vcf} {params.annovar_db} \
            -buildver {params.buildver} \
            -out {params.out_prefix} \
            -remove \
            -protocol refGene,cytoBand,clinvar_20240917,cosmic101,gnomad41_genome,dbnsfp47a,avsnp154,exac03,dbscsnv11 \
            -operation g,r,f,f,f,f,f,f,f \
            -nastring . \
            -vcfinput \
            -polish \
            2>&1 | tee {log}
        """


rule vardict_track_mrd:
    """
    Rastreia variantes de MRD ao longo do tempo
    Gera relatório de monitoramento
    """
    input:
        mrd_report = config["output"] + "output_vkaryo/vardict/{sample}_mrd_summary.txt",
        annotated = config["output"] + "output_vkaryo/vardict/annovar/{sample}_vardict.hg38_multianno.txt"
    output:
        tracking = config["output"] + "output_vkaryo/vardict/{sample}_mrd_tracking.tsv",
        plot = config["output"] + "output_vkaryo/vardict/{sample}_mrd_plot.pdf"
    log:
        config["output"] + "logs/vardict_tracking/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    script:
        "../scripts/track_mrd_variants.py"
