rule annotate_gnomad:
    """
    Anota variantes Mutect2 com frequências alélicas do gnomAD
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_pass.vcf.gz",
        gnomad = config["refgenome"]["gnomad_genome"]
    output:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_gnomad_annotated.vcf.gz",
        tbi = config["output"] + "output_vkaryo/mutect2/{sample}_gnomad_annotated.vcf.gz.tbi"
    log:
        config["output"] + "logs/filter_germline/{sample}_annotate.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        bcftools annotate \
            -a {input.gnomad} \
            -c INFO/AF \
            -o {output.vcf} \
            -O z \
            {input.vcf} \
            2>&1 | tee {log}

        tabix -f -p vcf {output.vcf}
        """


# COSMIC annotation is optional - check if file exists in config
if config["refgenome"].get("cosmic_heme") and os.path.exists(config["refgenome"].get("cosmic_heme", "")):
    rule annotate_cosmic:
        """
        Anota variantes com COSMIC hematológico (OPCIONAL)
        """
        input:
            vcf = config["output"] + "output_vkaryo/mutect2/{sample}_gnomad_annotated.vcf.gz",
            cosmic = config["refgenome"]["cosmic_heme"]
        output:
            vcf = config["output"] + "output_vkaryo/mutect2/{sample}_cosmic_annotated.vcf.gz",
            tbi = config["output"] + "output_vkaryo/mutect2/{sample}_cosmic_annotated.vcf.gz.tbi"
        log:
            config["output"] + "logs/filter_germline/{sample}_cosmic.log"
        conda:
            "../envs/gatk.yaml"
        shell:
            """
            bcftools annotate \
                -a {input.cosmic} \
                -c ID \
                -m COSMIC \
                -o {output.vcf} \
                -O z \
                {input.vcf} \
                2>&1 | tee {log}

            tabix -f -p vcf {output.vcf}
            """

    # Input para filtro germinal vem do COSMIC annotado
    _filter_input = config["output"] + "output_vkaryo/mutect2/{sample}_cosmic_annotated.vcf.gz"
    _cosmic_filter = " | INFO/COSMIC=1"
else:
    # Sem COSMIC - input vem direto do gnomAD annotado
    _filter_input = config["output"] + "output_vkaryo/mutect2/{sample}_gnomad_annotated.vcf.gz"
    _cosmic_filter = ""


rule filter_germline_variants:
    """
    Filtra variantes germinais usando gnomAD
    Mantém variantes com AF < 0.001 (0.1%) OU presentes no COSMIC (se disponível)
    """
    input:
        vcf = _filter_input
    output:
        somatic = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_filtered.vcf.gz",
        somatic_tbi = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_filtered.vcf.gz.tbi",
        germline = config["output"] + "output_vkaryo/mutect2/{sample}_germline_filtered.vcf.gz",
        germline_tbi = config["output"] + "output_vkaryo/mutect2/{sample}_germline_filtered.vcf.gz.tbi"
    params:
        af_threshold = config["variant_calling"]["min_vaf"],
        cosmic_filter = _cosmic_filter
    log:
        config["output"] + "logs/filter_germline/{sample}_filter.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        # Filtrar variantes somáticas (AF < 0.001 ou COSMIC presente se disponível)
        bcftools view \
            -i 'INFO/AF < 0.001{params.cosmic_filter}' \
            -O z \
            -o {output.somatic} \
            {input.vcf} \
            2>&1 | tee {log}

        tabix -f -p vcf {output.somatic}

        # Filtrar variantes germinais (AF >= 0.001)
        bcftools view \
            -i 'INFO/AF >= 0.001' \
            -O z \
            -o {output.germline} \
            {input.vcf} \
            2>&1 >> {log}

        tabix -f -p vcf {output.germline}
        """


rule extract_driver_mutations:
    """
    Extrai mutações em genes drivers conhecidos de leucemia
    Aplicando filtros mais relaxados para VAF baixo
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_filtered.vcf.gz",
        drivers = config["variant_calling"]["driver_genes"]
    output:
        drivers = config["output"] + "output_vkaryo/mutect2/{sample}_driver_mutations.vcf.gz",
        drivers_tbi = config["output"] + "output_vkaryo/mutect2/{sample}_driver_mutations.vcf.gz.tbi",
        other = config["output"] + "output_vkaryo/mutect2/{sample}_other_mutations.vcf.gz",
        other_tbi = config["output"] + "output_vkaryo/mutect2/{sample}_other_mutations.vcf.gz.tbi"
    params:
        min_vaf = config["variant_calling"]["min_vaf"],
        min_depth = config["variant_calling"]["min_depth"]
    log:
        config["output"] + "logs/filter_germline/{sample}_drivers.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        # Criar arquivo BED no formato correto: chr start end gene_name
        # Input format: gene_name chr start end
        grep -v "^#" {input.drivers} | awk '{{print $2, $3, $4, $1}}' OFS='\t' > {output.drivers}.bed

        # Extrair variantes em genes drivers (VAF >= 5%, DP >= 20)
        bcftools view \
            -R {output.drivers}.bed \
            -i "FORMAT/AF >= {params.min_vaf} & FORMAT/DP >= {params.min_depth}" \
            -O z \
            -o {output.drivers} \
            {input.vcf} \
            2>&1 | tee {log}

        tabix -f -p vcf {output.drivers}

        # Extrair outras variantes (não-drivers)
        bcftools view \
            -T ^{output.drivers}.bed \
            -O z \
            -o {output.other} \
            {input.vcf} \
            2>&1 >> {log}

        tabix -f -p vcf {output.other}

        # Limpar arquivo temporário
        rm -f {output.drivers}.bed
        """
