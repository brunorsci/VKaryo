rule filter_alt_contig_probes:
    """
    Filtra sondas do BED de captura (Agilent SureSelect V6) que referenciam
    contigs alternativos ("_alt"/"_random"/chrUn_...) ausentes do FASTA de
    referência do pipeline (GRCh38_character_genome.fasta, 195 contigs, vs.
    os 455 do hg38.fa completo usado para desenhar o kit) -- ver Seção 6.6
    da tese / scripts/filter_alt_contig_probes.py.
    """
    input:
        bed = config["refgenome"]["exome_bed"],
        fai = config["refgenome"]["ref"] + ".fai"
    output:
        noalt_bed = config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.noalt.bed"
    log:
        config["output"] + "logs/cnvkit_reference/filter_alt_contig_probes.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/filter_alt_contig_probes.py \
            --bed {input.bed} \
            --fai {input.fai} \
            --output {output.noalt_bed} \
            2>&1 | tee {log}
        """


rule cnvkit_target:
    """
    Gera o BED de alvos (target.bed) do CNVkit a partir do BED de captura
    já filtrado (sem sondas de contig alternativo), anotado com nomes de
    genes via refFlat.
    """
    input:
        bed = config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.noalt.bed",
        refflat = config["refgenome"]["ref_dir"] + "refFlat_hg38.txt"
    output:
        target_bed = config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.target.bed"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_target.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py target {input.bed} \
            --annotate {input.refflat} \
            --split \
            -o {output.target_bed} \
            2>&1 | tee {log}
        """


rule cnvkit_antitarget:
    """
    Gera o BED de antialvos (antitarget.bed) do CNVkit a partir do target.bed
    -- regiões fora do capture kit usadas como sinal de fundo/off-target.
    """
    input:
        target_bed = config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.target.bed"
    output:
        antitarget_bed = config["output"] + "output_vkaryo/cnvkit/Agilent_SureSelect_V6_hg38.antitarget.bed"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_antitarget.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py antitarget {input.target_bed} \
            -o {output.antitarget_bed} \
            2>&1 | tee {log}
        """
