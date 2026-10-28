#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cariotipo virtual: traco de numero de copias de genoma inteiro, no layout
usado pelas figuras publicadas de array-CGH e de WGS.

POR QUE ESTA FIGURA EXISTE
--------------------------
O workflow afirma produzir um cariotipo virtual e entrega duas tabelas: a de
eventos focais (rules/cnvFocalEvents.snakefile) e a de desequilibrio por braco
(rules/cnvArmLevel.snakefile). As duas sao chamadas -- conclusoes. Falta a
EVIDENCIA na forma em que a citogenetica a le: a nuvem de janelas ao longo do
genoma, onde uma trissomia ou uma perda de braco se ve porque a nuvem inteira
sai do zero. Um ideograma com barras coloridas nao e isso: ele desenha a
conclusao, nao o dado.

A METRICA E EXATAMENTE A DO RESTO DO PIPELINE
---------------------------------------------
Por janela: log2((observado + 1) / (esperado + 1)), com o esperado vindo da
mediana da coorte em scripts/cnv_bin_statistics.py, menos a MEDIANA
AUTOSSOMICA da propria amostra. E a mesma expressao de
scripts/cnv_detect_arm_level.py, linha por linha, de proposito: a figura
mostra a grandeza em que as tabelas foram chamadas, portanto o leitor pode
conferir uma contra a outra. Nao e um segundo metodo nem uma segunda
normalizacao.

NAO se usa a coluna log2 do .cnr. Aquela razao vem de referencia PLANA e, como
ja esta documentado em config.yaml, perde validade quantitativa quando a
cobertura tende a zero: um log2 real em torno de -5 pode ser reportado como -23. A
mediana da coorte por janela e, na pratica, um painel de normais em espaco de
cobertura, e e o que torna o traco legivel.

O QUE A FIGURA NAO TRAZ
-----------------------
Nada derivado de chamada: sem sombreado, sem seta, sem rotulo ISCN. Desenhar o
veredito sobre a imagem a torna circular como evidencia -- ela passaria a
exibir uma conclusao em vez de permitir que o leitor chegue a uma. As chamadas
estao nas tabelas de cnv_arm/ e cnv_focal/; aqui ficam o dado e o eixo.

EXCLUSOES, DECLARADAS E NAO SILENCIOSAS
---------------------------------------
- Bracos curtos acrocentricos (13p, 14p, 15p, 21p, 22p): satelite e rDNA, sem
  sequencia exonica unica. As poucas janelas que mapeiam ali sao artefato de
  repeticao e dominam o eixo y, mostrando uma delecao dramatica que nao
  existe. A citogenetica convencional tambem nao os pontua. A exclusao sai
  escrita no titulo da figura.
- chrX e chrY por padrao: a cobertura esperada por janela vem da mediana da
  coorte, que MISTURA amostras com uma e com duas copias de X -- a mesma
  ressalva ja registrada em cnv_detect_arm_level.py. Com --incluir-sexuais a
  figura os inclui, e ai o desvio deles nao distingue sexo de alteracao.
- Janelas com cobertura esperada abaixo de --cobertura-esperada-minima:
  ausencia de leitura ali nao informa nada.
"""
import argparse
import csv
import math
import os
import statistics
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ACROCENTRICOS = {"chr13", "chr14", "chr15", "chr21", "chr22"}
ORDEM = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def carregar_bandas(caminho):
    """chrom -> lista ordenada de (inicio, fim, nome, giemsa).

    Aceita tanto o cytoBand do UCSC quanto a tabela de bandas do Arriba
    (reference/arriba/cytobands_hg38_*.tsv), que vem com cabecalho e sem o
    prefixo 'chr' -- as duas estao em uso no repositorio.
    """
    bandas = defaultdict(list)
    with open(caminho) as fh:
        primeira = fh.readline()
        if not primeira.lower().startswith(('contig', 'chrom', '#')):
            fh.seek(0)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            if len(c) < 5:
                continue
            chrom = c[0] if c[0].startswith('chr') else f'chr{c[0]}'
            try:
                bandas[chrom].append((int(c[1]), int(c[2]), c[3], c[4]))
            except ValueError:
                continue
    for chrom in bandas:
        bandas[chrom].sort()
    return bandas


def centromero(bandas_chrom):
    acen = [b for b in bandas_chrom if b[3] == 'acen']
    if not acen:
        return None
    return (min(b[0] for b in acen), max(b[1] for b in acen))


def comprimentos(bandas):
    return {c: max(b[1] for b in bs) for c, bs in bandas.items() if bs}


def trechos_descartados(bandas):
    """(inicio, fim) do braco curto de cada acrocentrico -- ver cabecalho."""
    fora = {}
    for chrom in ACROCENTRICOS:
        cen = centromero(bandas.get(chrom, []))
        if cen:
            fora[chrom] = (0, cen[0])
    return fora


def carregar_estatisticas(caminho, cobertura_minima):
    """(chrom, start) -> cobertura esperada. Mesmo recorte de
    cnv_detect_arm_level.py: so bins de alvo, so acima da cobertura minima."""
    est = {}
    with open(caminho) as fh:
        cab = fh.readline().rstrip('\n').split('\t')
        i = {c: n for n, c in enumerate(cab)}
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            if c[i['tipo_de_bin']] != 'alvo':
                continue
            esperada = float(c[i['profundidade_mediana_absoluta']])
            if esperada < cobertura_minima:
                continue
            est[(c[i['chrom']], int(c[i['start']]))] = esperada
    return est


def razoes(cnr, est, cromossomos, fora):
    """chrom -> lista de (posicao_central, razao_log2), e a mediana
    autossomica que serve de base, como em cnv_detect_arm_level.py."""
    pontos = {c: [] for c in cromossomos}
    autossomicas = []
    with open(cnr) as fh:
        next(fh)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            chrom = c[0]
            try:
                a = int(c[1]); b = int(c[2]); obs = float(c[4])
            except (ValueError, IndexError):
                continue
            esperada = est.get((chrom, a))
            if esperada is None:
                continue
            r = math.log2((obs + 1.0) / (esperada + 1.0))
            if chrom not in ('chrX', 'chrY'):
                autossomicas.append(r)
            if chrom not in pontos:
                continue
            meio = (a + b) / 2
            trecho = fora.get(chrom)
            if trecho and trecho[0] <= meio < trecho[1]:
                continue
            pontos[chrom].append((meio, r))
    return pontos, autossomicas


def dlrs(med):
    """Dispersao TECNICA da figura, pelo espalhamento da razao log2 derivada
    (DLRS: derivative log ratio spread), metrica padrao de controle de
    qualidade em array-CGH.

    Por que nao o MAD simples das medianas de janela: aquele mede o quanto a
    figura inteira se espalha, e num genoma muito rearranjado quase todo esse
    espalhamento e BIOLOGIA. Em uma amostra muito rearranjada
    (mieloma hipotriploide) da MAD 0,599 com um traco visivelmente limpo --
    platos nitidos e pontos de quebra abruptos. Chamar aquilo de ruido seria
    mentira, e o numero no titulo faria o leitor desconfiar justamente da
    figura mais informativa do painel.

    O DLRS olha a diferenca entre janelas VIZINHAS. Um segmento alterado, por
    mais extenso que seja, nao contribui: dentro dele as vizinhas valem o
    mesmo. So os pontos de quebra contribuem, e eles sao poucos, o que a
    mediana absorve. O que sobra e a dispersao ponto a ponto, que e o que o
    olho precisa vencer para enxergar um deslocamento.
    """
    difs = []
    for (chrom, k), v in med.items():
        viz = med.get((chrom, k + 1))
        if viz is not None:
            difs.append(abs(viz - v))
    if not difs:
        return float('nan')
    return 1.4826 * statistics.median(difs) / math.sqrt(2)


def medianas_por_janela(pontos, janela_mb, minimo):
    med = {}
    for chrom, ps in pontos.items():
        baldes = defaultdict(list)
        for meio, r in ps:
            baldes[int(meio // (janela_mb * 1e6))].append(r)
        for k, v in baldes.items():
            if len(v) >= minimo:
                med[(chrom, k)] = statistics.median(v)
    return med


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cnr', required=True)
    ap.add_argument('--bin-statistics', required=True)
    ap.add_argument('--cytobands', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--out', required=True, help='PNG do traco de genoma inteiro')
    ap.add_argument('--dir-por-cromossomo', default=None,
                    help='Tambem escreve um painel por cromossomo aqui. O traco de '
                         'genoma inteiro comprime um cromossomo em algumas dezenas de '
                         'pixels -- suficiente para ver deslocamento do cromossomo '
                         'inteiro, insuficiente para ver um braco ou um segmento '
                         'sub-braco.')
    ap.add_argument('--janela-mb', type=float, default=1.0)
    ap.add_argument('--janela-mb-cromossomo', type=float, default=0.25,
                    help='Menor que a do genoma inteiro porque o eixo fica ~50x mais '
                         'largo por base, entao o traco resolve um evento mais curto '
                         'antes de o ruido dominar.')
    ap.add_argument('--min-janelas-balde', type=int, default=10)
    ap.add_argument('--cobertura-esperada-minima', type=float, default=5.0)
    ap.add_argument('--ylim', type=float, default=2.0)
    ap.add_argument('--incluir-sexuais', action='store_true')
    ap.add_argument('--titulo', default=None)
    args = ap.parse_args()

    bandas = carregar_bandas(args.cytobands)
    if not bandas:
        sys.exit('ERRO: arquivo de bandas vazio')
    tam = comprimentos(bandas)
    fora = trechos_descartados(bandas)
    est = carregar_estatisticas(args.bin_statistics, args.cobertura_esperada_minima)
    if not est:
        sys.exit('ERRO: estatisticas de bin vazias')

    cromossomos = [c for c in ORDEM if c in tam
                   and (args.incluir_sexuais or c not in ('chrX', 'chrY'))]
    pontos, autossomicas = razoes(args.cnr, est, cromossomos, fora)
    if not autossomicas:
        sys.exit('ERROR: no usable autosomal window')
    base = statistics.median(autossomicas)

    med = medianas_por_janela(pontos, args.janela_mb, args.min_janelas_balde)
    # Medido entre medianas de janela, nao entre bins crus: um punhado de bins
    # saturados nao pode definir a escala de dispersao da figura. E pelo DLRS,
    # nao pelo espalhamento total -- ver a funcao.
    disp = dlrs(med)

    deslocamento, x0, marcas = {}, 0.0, []
    for c in cromossomos:
        deslocamento[c] = x0
        marcas.append((x0 + tam[c] / 2e6, c.replace('chr', '')))
        x0 += tam[c] / 1e6
    total = x0

    fig, ax = plt.subplots(figsize=(22, 4.2))
    for c in cromossomos:
        ax.scatter([deslocamento[c] + m / 1e6 for m, _ in pontos[c]],
                   [r - base for _, r in pontos[c]],
                   s=0.7, c='#333333', alpha=0.18, linewidths=0, rasterized=True)
    for c in cromossomos:
        ks = sorted(k for k in med if k[0] == c)
        if not ks:
            continue
        ax.plot([deslocamento[c] + (k[1] + 0.5) * args.janela_mb for k in ks],
                [med[k] - base for k in ks],
                color='#c0392b', linewidth=0.9, solid_joinstyle='round')
    for c in cromossomos:
        ax.axvline(deslocamento[c], color='#bbbbbb', linewidth=0.6, zorder=0)
    ax.axhline(0, color='black', linewidth=0.7)
    for nivel, rotulo in ((0.585, '3/2'), (-1.0, '1/2')):
        ax.axhline(nivel, color='#7f8c8d', linewidth=0.5, linestyle=':')
        ax.text(total * 1.002, nivel, rotulo, fontsize=7, va='center', color='#7f8c8d')

    ax.set_xlim(0, total)
    ax.set_ylim(-args.ylim, args.ylim)
    ax.set_xticks([m[0] for m in marcas])
    ax.set_xticklabels([m[1] for m in marcas], fontsize=8)
    ax.set_ylabel('Log$_2$ (observed / expected)', fontsize=10)
    ax.set_title(f"{args.titulo or args.sample}    "
                 f"[{args.janela_mb:g} Mb windows, DLRS {disp:.3f}, "
                 f"one-copy gain = {0.585 / disp:.1f}x the spread"
                 + ('; acrocentric short arms excluded]' if fora else ']'),
                 fontsize=10, fontweight='bold')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=160)
    plt.close(fig)
    print(f"{args.sample}: DLRS={disp:.4f} janelas={len(med)} -> {args.out}")

    if not args.dir_por_cromossomo:
        return
    os.makedirs(args.dir_por_cromossomo, exist_ok=True)
    jw = args.janela_mb_cromossomo
    med_c = medianas_por_janela(pontos, jw, max(4, args.min_janelas_balde // 3))
    for c in cromossomos:
        if not pontos[c]:
            continue
        f2, a2 = plt.subplots(figsize=(16, 4.0))
        a2.scatter([m / 1e6 for m, _ in pontos[c]], [r - base for _, r in pontos[c]],
                   s=2.5, c='#333333', alpha=0.22, linewidths=0, rasterized=True)
        ks = sorted(k[1] for k in med_c if k[0] == c)
        if ks:
            a2.plot([(k + 0.5) * jw for k in ks], [med_c[(c, k)] - base for k in ks],
                    color='#c0392b', linewidth=1.3)
        # O centromero entra porque um evento de braco e lido contra ele: um
        # deslocamento que para no centromero e achado diferente de um que o
        # atravessa.
        cen = centromero(bandas.get(c, []))
        if cen:
            a2.axvspan(cen[0] / 1e6, cen[1] / 1e6, color='#95a5a6', alpha=0.25, zorder=0)
        a2.axhline(0, color='black', linewidth=0.7)
        for nivel, rotulo in ((0.585, '3/2'), (-1.0, '1/2')):
            a2.axhline(nivel, color='#7f8c8d', linewidth=0.5, linestyle=':')
            a2.text(tam[c] / 1e6 * 1.003, nivel, rotulo, fontsize=8,
                    va='center', color='#7f8c8d')
        a2.set_xlim(0, tam[c] / 1e6)
        a2.set_ylim(-args.ylim, args.ylim)
        a2.set_xlabel('Position (Mb)', fontsize=10)
        a2.set_ylabel('Log$_2$ (observed / expected)', fontsize=10)
        a2.set_title(f"{args.titulo or args.sample} — {c}   "
                     f"[{jw:g} Mb windows, genome DLRS {disp:.3f}]",
                     fontsize=11, fontweight='bold')
        for sp in ('top', 'right'):
            a2.spines[sp].set_visible(False)
        f2.tight_layout()
        f2.savefig(os.path.join(args.dir_por_cromossomo, f"{c}.png"), dpi=170)
        plt.close(f2)
    print(f"{args.sample}: {len(cromossomos)} paineis por cromossomo "
          f"-> {args.dir_por_cromossomo}")


if __name__ == '__main__':
    main()
