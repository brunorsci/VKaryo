# ============================================================================
# Eventos focais de número de cópias apoiados em evidência de cobertura
# ============================================================================
# POR QUÊ ESTAS REGRAS EXISTEM
#
# Com referência plana (flat), a razão log2 do CNVkit perde validade
# quantitativa quando a cobertura do segmento tende a zero: o valor reportado
# é o piso numérico do método, não a magnitude da perda. Um segmento cuja
# cobertura observada corresponde a log2 em torno de -5 pode ser reportado
# como -23, erro de quase vinte unidades de log2. No lado dos ganhos não há
# esse problema: a razão nativa e a razão de cobertura concordam.
#
# Consequência: a deleção homozigótica -- que é o achado de maior peso clínico
# -- não podia ser chamada nem graduada pelo workflow. Estas regras suprem
# isso, avaliando cada bin pela razão entre cobertura observada e esperada.
#
# A cobertura esperada por bin vem de `cnv_bin_statistics`, derivada da PRÓPRIA
# coorte em execução, com o kit de captura dessa execução. Não é arquivo fixo:
# um kit diferente tem bins diferentes, e cobertura esperada de outro kit
# estaria errada. Quem analisa uma amostra isolada aponta
# `cnv.bin_statistics_file` no config para um arquivo gerado antes com o mesmo
# kit -- é o caminho previsto para uso clínico amostra a amostra.
#
# Genérico por construção: nenhuma destas regras usa lista de genes, painel de
# doença ou pressuposto de neoplasia. Serve a câncer e a não câncer.
# ============================================================================


rule cnv_bin_statistics:
    """Cobertura esperada e frequência de queda de cada bin, derivadas da
    coorte em execução.

    Requer um número mínimo de amostras para a mediana ser estável
    (`cnv.bin_statistics_min_samples`, padrão 20). Abaixo disso, o caminho
    correto é apontar `cnv.bin_statistics_file` para um arquivo gerado
    previamente com o mesmo kit.

    RESSALVA A DECLARAR EM QUALQUER PUBLICAÇÃO: a mediana é da coorte em
    execução. Numa coorte em que uma região é recorrentemente alterada, a
    cobertura esperada dessa região carrega essa alteração. A mediana de uma
    coorte heterogênea atenua muito o efeito, não o elimina. Por isso estas
    estatísticas servem para julgar CONFIABILIDADE do bin, e não como
    referência de número de cópias.
    """
    input:
        cnr=expand(config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr",
                   sample=SAMPLES),
    output:
        stats=config["output"] + "output_vkaryo/cnvkit/cnv_bin_statistics.tsv",
    params:
        cnr_dir=config["output"] + "output_vkaryo/cnvkit/",
        min_samples=config["cnv"].get("bin_statistics_min_samples", 20),
    log:
        config["output"] + "logs/cnv_focal/bin_statistics.log",
    threads: 1
    shell:
        """
        python3 scripts/cnv_bin_statistics.py \
            --cnr-dir {params.cnr_dir} \
            --min-amostras {params.min_samples} \
            --out {output.stats} \
            2>&1 | tee {log}
        """


def bin_statistics_for(wildcards):
    """Arquivo fixo declarado no config, quando houver; senão, o derivado da
    coorte em execução."""
    declarado = config["cnv"].get("bin_statistics_file")
    if declarado:
        return declarado
    return config["output"] + "output_vkaryo/cnvkit/cnv_bin_statistics.tsv"


rule cnv_focal_events:
    """Eventos focais por amostra, com a evidência de cobertura ao lado."""
    input:
        cnr=config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr",
        stats=bin_statistics_for,
    output:
        tsv=config["output"] + "output_vkaryo/cnv_focal/{sample}_focal_events.tsv",
    params:
        homo=config["cnv"].get("focal_homozygous_ratio", 0.05),
        perda=config["cnv"].get("focal_loss_ratio", 0.70),
        ganho=config["cnv"].get("focal_gain_ratio", 1.50),
        cob_min=config["cnv"].get("focal_min_expected_coverage", 5.0),
        min_bins=config["cnv"].get("focal_min_bins", 3),
        max_dist=config["cnv"].get("focal_max_bin_gap", 100000),
        max_freq=config["cnv"].get("focal_max_cohort_frequency", 0.05),
        qc_min=config["cnv"].get("focal_min_sample_global_ratio", 0.30),
    log:
        config["output"] + "logs/cnv_focal/{sample}.log",
    threads: 1
    shell:
        """
        python3 scripts/cnv_detect_focal_events.py \
            --cnr {input.cnr} \
            --bin-statistics {input.stats} \
            --sample {wildcards.sample} \
            --limiar-homozigotico {params.homo} \
            --limiar-perda {params.perda} \
            --limiar-ganho {params.ganho} \
            --cobertura-esperada-minima {params.cob_min} \
            --min-bins {params.min_bins} \
            --max-distancia {params.max_dist} \
            --max-frequencia-coorte {params.max_freq} \
            --razao-global-minima {params.qc_min} \
            --out {output.tsv} \
            2>&1 | tee {log}
        """


rule cnv_focal_events_consolidated:
    """Uma tabela com os eventos focais de todas as amostras da execução."""
    input:
        tsv=expand(config["output"] + "output_vkaryo/cnv_focal/{sample}_focal_events.tsv",
                   sample=SAMPLES),
    output:
        tsv=config["output"] + "output_vkaryo/cnv_focal/all_samples_focal_events.tsv",
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
