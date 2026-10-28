rule annotsvDUO:
    """
    Manta já entra completo (manta/{sample}/{sample}_fixed.vcf não é
    filtrado por PASS em nenhuma etapa anterior -- confirmado). DELLY agora
    usa {sample}_sv_full.vcf.gz (rule delly_full_vcf, PASS e não-PASS) em
    vez do antigo {sample}_sv.vcf.gz (PASS-only) -- ver
    rules/delly_filter.snakefile para o porquê.
    """
    input:
        manta=config["output"] + "output_vkaryo/manta/{sample}/{sample}_fixed.vcf",
        delly=config["output"] + "output_vkaryo/delly/svs/{sample}_sv_full.vcf.gz",
    output:
        manta=config["output"] + "output_vkaryo/annotsv/manta/{sample}_annotsv_manta.tsv",
        delly=config["output"] + "output_vkaryo/annotsv/delly/{sample}_annotsv_delly.tsv",
    params:
        output_dir=config["output"] + "output_vkaryo/annotsv/",
        db="reference/AnnotSV_DB",  # Diretório contendo bases de dados necessárias
    # AnnotSV (tclsh) é efetivamente monothread: medido em 2026-10-21 durante
    # a re-anotação do painel, usa ~0,5 core por job. O valor anterior
    # (threads: 60) não acelerava nada e ainda forçava o Snakemake a rodar UMA
    # amostra por vez em máquina de 48 cores -- ~3 min/amostra em série, ~5 h
    # para 100 amostras. Com 2, roda ~24 em paralelo sem oversubscrição real.
    threads: 2
    conda:
        "../envs/annotsv.yaml"
    shell:
        """
        echo "##############################################"
        echo "-----------    Running AnnotSV    ------------"
        echo "##############################################"

        # NOTA sobre o diretório temporário por amostra: o AnnotSV grava seus
        # arquivos intermediários NO DIRETÓRIO DO ARQUIVO DE SAÍDA, com nome
        # derivado do timestamp em segundos
        # (ex.: 20261021-014722_AnnotSV_inputSVfile.formatted.sorted.bed).
        # Com várias amostras em paralelo, duas que iniciem no mesmo segundo
        # colidem nesse nome e uma delas morre com "Unable to open file ...
        # Exiting" (observado em 2026-10-21, 22 de 24 jobs perdidos). Dar a
        # cada amostra um diretório próprio elimina a colisão sem precisar
        # serializar a execução.
        TMP_MANTA=$(dirname {output.manta})/_tmp_{wildcards.sample}_manta
        TMP_DELLY=$(dirname {output.delly})/_tmp_{wildcards.sample}_delly
        rm -rf "$TMP_MANTA" "$TMP_DELLY"
        mkdir -p "$TMP_MANTA" "$TMP_DELLY"

        ## MANTA results
        AnnotSV \
            -SVinputFile {input.manta} \
            -outputDir "$TMP_MANTA" \
            -outputFile $(basename {output.manta}) \
            -annotationsDir {params.db}/AnnotSV_annotations \
            -genomeBuild GRCh38

        ## DELLY results
        AnnotSV \
            -SVinputFile {input.delly} \
            -outputDir "$TMP_DELLY" \
            -outputFile $(basename {output.delly}) \
            -annotationsDir {params.db}/AnnotSV_annotations \
            -genomeBuild GRCh38

        mv "$TMP_MANTA"/$(basename {output.manta}) {output.manta}
        mv "$TMP_DELLY"/$(basename {output.delly}) {output.delly}
        # mantém os .unannotated.tsv ao lado da saída principal, como antes
        for extra in "$TMP_MANTA"/*.unannotated.tsv; do
            [ -e "$extra" ] && mv "$extra" $(dirname {output.manta})/ || true
        done
        for extra in "$TMP_DELLY"/*.unannotated.tsv; do
            [ -e "$extra" ] && mv "$extra" $(dirname {output.delly})/ || true
        done
        rm -rf "$TMP_MANTA" "$TMP_DELLY"
        """
