# ============================================================================
# Cariotipo virtual -- a figura
# ============================================================================
# POR QUE ESTA REGRA EXISTE
#
# As duas etapas de CNV deste workflow entregam TABELAS: cnv_focal_events
# (perda homozigotica e ganho de alto nivel) e cnv_arm_level (desequilibrio
# hemizigotico por braco). As duas sao chamadas -- conclusoes ja tomadas.
# Faltava a EVIDENCIA na forma em que a citogenetica a le: a nuvem de janelas
# ao longo do genoma, em que uma trissomia ou uma perda de braco aparece
# porque a nuvem inteira sai do zero.
#
# O ideograma de barras (rules/plot_aneuploidy.snakefile) nao cumpre esse
# papel: ele pinta a conclusao, nao mostra o dado. O scatter de genoma inteiro
# (rules/cnv_scatter.snakefile) usa a coluna log2 do .cnr, vinda de referencia
# PLANA -- a mesma razao que config.yaml ja documenta como quantitativamente
# invalida no regime de cobertura baixa (um segmento cuja cobertura
# reporta -23,4).
#
# Esta regra desenha a MESMA grandeza em que as tabelas foram chamadas:
# log2((observado + 1)/(esperado + 1)) contra a mediana da coorte por janela,
# menos a mediana autossomica da amostra -- expressao identica a de
# scripts/cnv_detect_arm_level.py. O leitor pode conferir a figura contra a
# tabela sem trocar de metrica.
#
# QUANTO SE GANHA, MEDIDO. Nas MESMAS janelas de 1 Mb e nos mesmos
# autossomos, pelo DLRS (espalhamento da razao derivada, metrica padrao de
# array-CGH -- ver a funcao dlrs() no script e o porque de nao ser o MAD
# simples), a razao contra a estatistica da coorte apresenta cerca de metade
# da dispersao ponto a ponto da razao de referencia plana: valores tipicos em
# torno de 0,055 contra 0,115. Isso pela mesma leitura, sem realinhar nada.
# Com DLRS de 0,055, um ganho de uma copia (+0,585) fica a cerca de dez vezes
# a dispersao da figura.
#
# Os bracos que a tabela de cnv_arm chama sao os deslocamentos EXTENSOS da
# figura; os picos estreitos restantes sao eventos focais, que pertencem a
# outra tabela. O painel por cromossomo acrescenta o que a tabela de braco nao
# consegue dizer: uma perda pode nao tomar o braco inteiro, sendo
# intersticial e assentando sobre a linha de uma copia, com a porcao proximal
# preservada. A tabela chama o braco porque a mediana passa no limiar; a
# figura mostra a extensao real.
#
# NAO ENTRA NA rule all DE PROPOSITO. E figura de inspecao, gerada sob
# demanda: incluir no alvo padrao faria toda execucao reprocessar amostras ja
# analisadas. Pedir
# esta figura pela rule all arrastaria o DAG inteiro. Gere por alvo explicito:
#
#   snakemake --rerun-triggers mtime \
#     results/output_vkaryo/cariotipo_virtual/<amostra>/<amostra>_genoma.png
#
# conferindo antes com -n que so esta regra roda.
# ============================================================================


# Nome de linhagem para o titulo. NAO usa SAMPLE_TO_CELLLINE do Snakefile: no
# ponto em que ele chega ate aqui ele ja foi reatribuido para conter so as 70
# amostras com alteracao-chave documentada (ver Snakefile, bloco do
# cnv_scatter), e as outras 30 sairiam com o ERR repetido no titulo. Este dict
# e local a esta regra justamente para nao alterar aquele.
import csv as _csv
import os as _os

CARIOTIPO_LINHAGEM = {}
_mapa = sample_mapping_file  # definido no Snakefile, linha 54
if _os.path.exists(_mapa):
    with open(_mapa) as _f:
        for _r in _csv.DictReader(_f, delimiter="\t"):
            CARIOTIPO_LINHAGEM[_r["DNA_sample"]] = _r["Cell_line"]


rule cariotipo_virtual:
    """Traco de genoma inteiro + um painel por cromossomo, por amostra."""
    input:
        cnr=config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr",
        stats=bin_statistics_for,
        script="scripts/plot_cariotipo_virtual.py",
    output:
        png=config["output"] + "output_vkaryo/cariotipo_virtual/{sample}/{sample}_genoma.png",
        dir=directory(config["output"] + "output_vkaryo/cariotipo_virtual/{sample}/cromossomos"),
    params:
        cytobands=config["cnv"].get(
            "cytobands", "reference/arriba/cytobands_hg38_GRCh38_v2.5.0.tsv"),
        cob_min=config["cnv"].get("focal_min_expected_coverage", 5.0),
        # Mesma janela usada para medir o ruido declarado no titulo. 1 Mb no
        # genoma inteiro e 0,25 Mb por cromossomo: o eixo por cromossomo e
        # ~50x mais largo por base, entao resolve evento mais curto antes de o
        # ruido dominar.
        janela=config["cnv"].get("cariotipo_janela_mb", 1.0),
        janela_chrom=config["cnv"].get("cariotipo_janela_mb_cromossomo", 0.25),
        titulo=lambda wc: (
            f"{CARIOTIPO_LINHAGEM.get(wc.sample, wc.sample)} ({wc.sample})"
            " -- cariotipo virtual"),
    log:
        config["output"] + "logs/cariotipo_virtual/{sample}.log",
    threads: 1
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        # Interpretador resolvido por $CONDA_PREFIX, NAO por "python3" no PATH.
        # Um perfil de login pode prepender ao PATH o bin de outro ambiente
        # conda, e o bash que o Snakemake abre le
        # esse perfil DEPOIS de ativar o ambiente da regra -- entao "python3"
        # cai no interpretador errado, que nao tem matplotlib, embora o
        # ambiente certo esteja ativo e $CONDA_PREFIX aponte para ele. Chamar
        # o binario pelo prefixo e imune a ordem do PATH.
        #
        # Mesma nota de rules/cnv_scatter.snakefile quanto a libstdc++: o
        # matplotlib deste ambiente foi compilado contra uma mais nova
        # (CXXABI_1.3.15) do que a do sistema, e o linker dinamico acha a do
        # sistema primeiro.
        export LD_LIBRARY_PATH="${{CONDA_PREFIX}}/lib:${{LD_LIBRARY_PATH:-}}"
        "${{CONDA_PREFIX}}/bin/python" scripts/plot_cariotipo_virtual.py \
            --cnr {input.cnr} \
            --bin-statistics {input.stats} \
            --cytobands {params.cytobands} \
            --sample {wildcards.sample} \
            --titulo "{params.titulo}" \
            --cobertura-esperada-minima {params.cob_min} \
            --janela-mb {params.janela} \
            --janela-mb-cromossomo {params.janela_chrom} \
            --out {output.png} \
            --dir-por-cromossomo {output.dir} \
            2>&1 | tee {log}
        """
