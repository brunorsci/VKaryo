# ============================================================================
# Desequilibrio de numero de copias em nivel de braco cromossomico
# ============================================================================
# POR QUE ESTA REGRA EXISTE
#
# A regra cnv_focal_events (rules/cnvFocalEvents.snakefile) so reporta perda
# homozigotica e amplificacao de alto nivel -- os regimes em que a evidencia de
# cobertura e inequivoca. Eventos de braco sao, em sua maioria, HEMIZIGOTICOS:
# perder uma de duas copias de um braco nao zera a cobertura, so a reduz pela
# metade. Esses eventos apareciam na nuvem de janelas do cariograma mas nao
# entravam em tabela alguma.
#
# A verificacao direta da profundidade confirma a escala: um braco com uma
# copia contra duas fica em torno de 0,5 da cobertura de um braco equilibrado
# da mesma amostra, e uma monossomia completa se aproxima de -1,0 em log2. Um
# braco perdido pode conter, na copia remanescente, uma delecao homozigotica
# focal -- as duas camadas sao complementares e nao redundantes.
#
# CRITERIO alinhado a pratica publicada: amplitude mais extensao, como em
# Shukla et al., Nat Commun 2020 (que inclui linhagens celulares). Ver o
# cabecalho de scripts/cnv_detect_arm_level.py para a discussao e as fontes.
#
# Generico: nao usa lista de genes nem pressuposto de neoplasia.
# ============================================================================


rule cnv_arm_level:
    """Desequilibrio por braco, por amostra."""
    input:
        cnr=config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr",
        stats=bin_statistics_for,
    output:
        tsv=config["output"] + "output_vkaryo/cnv_arm/{sample}_arm_level.tsv",
    params:
        cytobands=config["cnv"].get(
            "cytobands", "reference/arriba/cytobands_hg38_GRCh38_v2.5.0.tsv"),
        min_janelas=config["cnv"].get("arm_min_bins", 50),
        cob_min=config["cnv"].get("focal_min_expected_coverage", 5.0),
        perda=config["cnv"].get("arm_loss_threshold", -0.30),
        ganho=config["cnv"].get("arm_gain_threshold", 0.25),
        amplitude=config["cnv"].get("arm_bin_amplitude", 0.20),
        fracao=config["cnv"].get("arm_min_fraction", 0.50),
    log:
        config["output"] + "logs/cnv_arm/{sample}.log",
    threads: 1
    shell:
        """
        python3 scripts/cnv_detect_arm_level.py \
            --cnr {input.cnr} \
            --bin-statistics {input.stats} \
            --cytobands {params.cytobands} \
            --sample {wildcards.sample} \
            --min-janelas {params.min_janelas} \
            --cobertura-esperada-minima {params.cob_min} \
            --limiar-perda {params.perda} \
            --limiar-ganho {params.ganho} \
            --limiar-amplitude-janela {params.amplitude} \
            --fracao-minima-do-braco {params.fracao} \
            --out {output.tsv} \
            2>&1 | tee {log}
        """


rule cnv_arm_level_consolidated:
    """Uma tabela com os bracos de todas as amostras da execucao."""
    input:
        tsv=expand(config["output"] + "output_vkaryo/cnv_arm/{sample}_arm_level.tsv",
                   sample=SAMPLES),
    output:
        tsv=config["output"] + "output_vkaryo/cnv_arm/all_samples_arm_level.tsv",
    threads: 1
    run:
        import os
        with open(output.tsv, 'w') as saida:
            cabecalho_escrito = False
            for caminho in sorted(input.tsv):
                if not os.path.exists(caminho):
                    continue
                with open(caminho) as entrada:
                    cabecalho = entrada.readline()
                    if not cabecalho_escrito:
                        saida.write(cabecalho)
                        cabecalho_escrito = True
                    for linha in entrada:
                        saida.write(linha)
