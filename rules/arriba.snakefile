rule arriba:
    input:
        bam=config["output"] + "output_vkaryo/star_aligned/{sample}_Aligned.sortedByCoord.out.bam",
        gtf=config["refgenome"]["gtf"],
        ref=config["refgenome"]["ref"],
        # Os tres arquivos de recurso do Arriba sao instalados pelo proprio
        # pacote, dentro do ambiente conda da regra. O prefixo do ambiente
        # varia por maquina, entao e resolvido em tempo de execucao a partir
        # de $CONDA_PREFIX, no shell abaixo -- nao como input, que o Snakemake
        # exigiria existir antes de criar o ambiente.
    output:
        fusions_tsv=config["output"] + "output_vkaryo/arriba/{sample}_fusions.tsv",
        fusions_rm_tsv=config["output"] + "output_vkaryo/arriba/{sample}_fusions_discarded.tsv",
    params:
        dir=config["output"] + "output_vkaryo/arriba/",
        chimeric=config["output"] + "output_vkaryo/star_aligned/{sample}_Aligned.Chimeric.out.sam",
    benchmark: config["output"] + "output_vkaryo/benchmarks/arriba/arriba_{sample}_benchmark.txt",

    conda:
        "../envs/arriba.yaml"
    threads: 20
    shell:
        """
        echo "##############################################" 
        echo "-----------    Running ARRIBA    ------------" 
        echo "##############################################"
        mkdir -p {params.dir}

        # Recursos do Arriba, dentro do ambiente conda desta regra.
        ARRIBA_DB="$CONDA_PREFIX/var/lib/arriba"

        arriba -x {input.bam} -g "{input.gtf}" -a "{input.ref}" \
        -b "$ARRIBA_DB/blacklist_hg38_GRCh38_v2.5.0.tsv.gz" \
        -k "$ARRIBA_DB/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz" \
        -p "$ARRIBA_DB/protein_domains_hg38_GRCh38_v2.5.0.gff3" \
        -o "{output.fusions_tsv}" \
        -O "{output.fusions_rm_tsv}" -@ {threads}
        """
