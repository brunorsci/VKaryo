rule mutect2_call:
    """
    Mutect2 para chamada de variantes somáticas (padrão GATK)
    Modo tumor-only com Panel of Normals e germline resource
    Otimizado para detecção de SNVs em leucemias
    """
    input:
        bam = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        bai = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bai",
        ref = config["refgenome"]["ref"],
        # NOTA: exome_bed cru referencia contigs alternativos ausentes do
        # FASTA de referencia (mesmo problema do CNVkit, corrigido do mesmo
        # jeito) -- GATK falha com "Badly formed genome unclippedLoc" se
        # usado diretamente em -L. Usa a versao ja filtrada
        # (scripts/filter_alt_contig_probes.py / rules/cnvkitTargets.snakefile).
        intervals = config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.noalt.bed",
        gnomad = config["refgenome"]["gnomad_genome"],
        pon = config["refgenome"]["panel_of_normals"]
    output:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_unfiltered.vcf.gz",
        tbi = config["output"] + "output_vkaryo/mutect2/{sample}_unfiltered.vcf.gz.tbi",
        stats = config["output"] + "output_vkaryo/mutect2/{sample}_unfiltered.vcf.gz.stats",
        f1r2 = config["output"] + "output_vkaryo/mutect2/{sample}_f1r2.tar.gz"
    params:
        sample = "{sample}",
        extra = "--max-mnp-distance 0 " +
                "--dont-use-soft-clipped-bases " +
                "--pcr-indel-qual 40"
    log:
        config["output"] + "logs/mutect2/{sample}.log"
    threads: 4
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.intervals} \
            --tumor-sample {params.sample} \
            --germline-resource {input.gnomad} \
            --panel-of-normals {input.pon} \
            {params.extra} \
            -O {output.vcf} \
            --f1r2-tar-gz {output.f1r2} \
            --native-pair-hmm-threads {threads} \
            2>&1 | tee {log}
        """


rule mutect2_learn_orientation:
    """
    Aprende modelo de orientação de artifacts (OxoG, FFPE)
    """
    input:
        f1r2 = config["output"] + "output_vkaryo/mutect2/{sample}_f1r2.tar.gz"
    output:
        model = config["output"] + "output_vkaryo/mutect2/{sample}_read-orientation-model.tar.gz"
    log:
        config["output"] + "logs/mutect2_orientation/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk LearnReadOrientationModel \
            -I {input.f1r2} \
            -O {output.model} \
            2>&1 | tee {log}
        """


rule mutect2_pileup_summaries:
    """
    Calcula pileup summaries para estimativa de contaminação
    Usa gnomAD common variants (AF > 0.05)
    """
    input:
        bam = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        bai = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bai",
        intervals = config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.noalt.bed",
        gnomad = config["refgenome"]["gnomad_genome"]
    output:
        pileup = config["output"] + "output_vkaryo/mutect2/{sample}_pileup.table"
    params:
        sample = "{sample}"
    log:
        config["output"] + "logs/mutect2_pileup/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.gnomad} \
            -L {input.intervals} \
            -O {output.pileup} \
            2>&1 | tee {log}
        """


rule mutect2_calculate_contamination:
    """
    Estima contaminação da amostra usando pileup summaries
    """
    input:
        pileup = config["output"] + "output_vkaryo/mutect2/{sample}_pileup.table"
    output:
        contamination = config["output"] + "output_vkaryo/mutect2/{sample}_contamination.table",
        segments = config["output"] + "output_vkaryo/mutect2/{sample}_segments.table"
    log:
        config["output"] + "logs/mutect2_contamination/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk CalculateContamination \
            -I {input.pileup} \
            -O {output.contamination} \
            --tumor-segmentation {output.segments} \
            2>&1 | tee {log}
        """


rule mutect2_filter:
    """
    Filtra variantes Mutect2 aplicando modelos de orientação e contaminação
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_unfiltered.vcf.gz",
        tbi = config["output"] + "output_vkaryo/mutect2/{sample}_unfiltered.vcf.gz.tbi",
        stats = config["output"] + "output_vkaryo/mutect2/{sample}_unfiltered.vcf.gz.stats",
        ref = config["refgenome"]["ref"],
        contamination = config["output"] + "output_vkaryo/mutect2/{sample}_contamination.table",
        segments = config["output"] + "output_vkaryo/mutect2/{sample}_segments.table",
        orientation = config["output"] + "output_vkaryo/mutect2/{sample}_read-orientation-model.tar.gz"
    output:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_filtered.vcf.gz",
        tbi = config["output"] + "output_vkaryo/mutect2/{sample}_filtered.vcf.gz.tbi"
    log:
        config["output"] + "logs/mutect2_filter/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            --stats {input.stats} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segments} \
            --ob-priors {input.orientation} \
            -O {output.vcf} \
            2>&1 | tee {log}

        # Create tabix index for filtered VCF (force overwrite)
        tabix -f -p vcf {output.vcf}
        """


rule mutect2_pass_only:
    """
    Extrai apenas variantes PASS
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_filtered.vcf.gz"
    output:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_pass.vcf.gz",
        tbi = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_pass.vcf.gz.tbi"
    log:
        config["output"] + "logs/mutect2_pass/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        bcftools view -f PASS {input.vcf} | bgzip -c > {output.vcf}
        tabix -f -p vcf {output.vcf}
        """


rule mutect2_funcregion_scan:
    """
    Estágio 1 (rápido) da anotação em 2 estágios -- ver
    PLANO_correcao_artigo_e_banco_dados.md. Roda SÓ o refGene (anotação
    baseada em gene, rápida mesmo com 800k+ variantes -- ao contrário dos
    bancos de filtro tipo gnomAD/dbNSFP, que dominam o tempo do ANNOVAR) no
    VCF COMPLETO pós-FilterMutectCalls (_filtered.vcf.gz, PASS e não-PASS).
    Objetivo único: descobrir qual posição é exônica/splicing, pra alimentar
    só essas no estágio 2 (caro). Testado: ~1 min mesmo com 818 mil
    variantes de entrada.
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_filtered.vcf.gz"
    output:
        multianno = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_funcregion.hg38_multianno.txt"
    params:
        out_prefix = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_funcregion",
        annovar_dir = config["refgenome"]["annovar_dir"],
        annovar_db = config["refgenome"]["annovar_db"],
        buildver = "hg38"
    log:
        config["output"] + "logs/mutect2_funcregion_scan/{sample}.log"
    threads: 2
    conda:
        "../envs/annovar.yaml"
    shell:
        """
        {params.annovar_dir}/table_annovar.pl {input.vcf} {params.annovar_db} \
            -buildver {params.buildver} \
            -out {params.out_prefix} \
            -remove \
            -protocol refGene \
            -operation g \
            -nastring . \
            -vcfinput \
            -polish \
            2>&1 | tee {log}
        """


rule mutect2_exonic_subset_vcf:
    """
    Estágio 1b -- extrai do VCF completo (PASS e não-PASS) só as posições
    classificadas como exônicas/splicing pelo estágio 1 acima. Filtro por
    REGIÃO FUNCIONAL, não por nome de gene -- substitui o antigo filtro por
    painel de genes curado (driver_genes_leukemia.txt), que escondia
    variantes reais em genes fora da lista: um painel curado nunca cobre todos
    os genes relevantes, e usá-lo como critério de inclusão descarta achados
    antes de qualquer revisão. Mantém o VCF (não avinput puro)
    pra preservar FILTER/FORMAT/SAMPLE que scripts/consolidate_program_results.py
    já usa pra classificar camada_confianca e extrair AF.
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_filtered.vcf.gz",
        multianno = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_funcregion.hg38_multianno.txt"
    output:
        bed = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_exonic_positions.bed",
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_exonic_subset.vcf.gz",
        tbi = config["output"] + "output_vkaryo/mutect2/{sample}_exonic_subset.vcf.gz.tbi"
    log:
        config["output"] + "logs/mutect2_exonic_subset_vcf/{sample}.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        awk -F'\t' '
            $6=="exonic" || $6=="splicing" || $6=="exonic;splicing" ||
            $6=="ncRNA_exonic" || $6=="ncRNA_splicing" || $6=="ncRNA_exonic;splicing" {{
                print $1"\t"($2-1)"\t"$3
            }}' {input.multianno} | sort -k1,1 -k2,2n -u > {output.bed}
        bcftools view -R {output.bed} {input.vcf} -Oz -o {output.vcf} 2>&1 | tee {log}
        tabix -f -p vcf {output.vcf}
        """


rule mutect2_annovar:
    """
    Estágio 2 (caro) -- anotação completa (ClinVar/COSMIC/gnomAD/dbNSFP/
    avSNP/ExAC/dbscSNV), mas só no subconjunto exônico/splicing do estágio 1
    (PASS e não-PASS) -- não mais em _somatic_pass.vcf.gz. Testado: ~9 min
    pro subconjunto (32.687 variantes) vs 50+ min pro VCF completo (818 mil)
    na mesma amostra -- reduz o custo em ~5-6x e ainda cobre variante
    não-PASS em região codificante (o que o filtro PASS-only antigo escondia).
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_exonic_subset.vcf.gz"
    output:
        annotated = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_mutect2.hg38_multianno.txt",
        vcf_annotated = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_mutect2.hg38_multianno.vcf"
    params:
        out_prefix = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_mutect2",
        annovar_dir = config["refgenome"]["annovar_dir"],
        annovar_db = config["refgenome"]["annovar_db"],
        buildver = "hg38"
    log:
        config["output"] + "logs/mutect2_annovar/{sample}.log"
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


rule mutect2_summary:
    """
    Gera sumário de variantes somáticas Mutect2
    Usa anotação ANNOVAR para análise de driver mutations
    """
    input:
        vcf = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_pass.vcf.gz",
        annotated = config["output"] + "output_vkaryo/mutect2/annovar/{sample}_mutect2.hg38_multianno.txt"
    output:
        summary = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_summary.txt",
        drivers = config["output"] + "output_vkaryo/mutect2/{sample}_driver_mutations.tsv"
    log:
        config["output"] + "logs/mutect2_summary/{sample}.log"
    resources:
        mem_mb = 6000
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/summarize_mutect2.py \
            --vcf {input.vcf} \
            --annotated {input.annotated} \
            --summary {output.summary} \
            --drivers {output.drivers} \
            2>&1 | tee {log}
        """
