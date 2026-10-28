rule cnvkit_reference:
    input:
        target_bed=config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.target.bed",
        antitarget_bed=config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.antitarget.bed",
        fasta=config["refgenome"]["ref"]
    output:
        reference=config["output"] + "output_vkaryo/cnvkit/flat_reference.cnn"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_reference.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

        # Cria a referência "flat" a partir dos BEDs de target/antitarget,
        # COM correção de GC/RepeatMasker (Seção 6.6 da tese) -- não usar
        # --no-gc --no-rmask aqui.
        #
        # NOTA: o FASTA de referência (GRCh38_character_genome.fasta, 195
        # contigs) não tem os contigs alternativos ("_alt"/"_random"/chrUn_...)
        # que 151 das 244655 sondas do kit Agilent SureSelect V6 referenciam
        # (o FASTA completo hg38.fa, usado para desenhar o kit, tem 455
        # contigs) -- extrair sequência desses contigs para calcular
        # GC/RepeatMasker lança KeyError se o BED não estiver filtrado antes.
        # Sem essa correção, a referência "flat" resultante produz um desvio
        # sistemático de log2 GLOBAL por amostra (quase todo o genoma "ganho"
        # ou "perdido" de uma vez -- ver karyotype_quality_flag em
        # classify_leukemia.py), inutilizando o componente de CNV das camadas
        # de estratificação de risco.
        #
        # target_bed/antitarget_bed acima já são derivados de
        # reference/Agilent_SureSelect_V6_hg38.noalt.bed (as 151 sondas de
        # contig alternativo removidas -- 0,06% do total, nenhum gene perdido
        # exclusivamente por elas; ver scripts/filter_alt_contig_probes.py),
        # então a extração de FASTA abaixo não encontra contigs ausentes
        # e a correção de GC/RepeatMasker roda normalmente.
        cnvkit.py reference \
            -t {input.target_bed} \
            -a {input.antitarget_bed} \
            -f {input.fasta} \
            -o {output.reference} \
            2>&1 | tee {log}
        """
