#!/usr/bin/env python3
"""Cariograma virtual de um cromossomo: nuvem de bins mais segmentos
coloridos, numa escala em que delecao homozigotica e amplificacao aparecem.

Estilo
------
Um cromossomo por figura, no formato de `cnvkit.py scatter` e do
um painel por cromossomo: posicao em Mb no eixo X, razao de copias
no eixo Y, nuvem cinza dos bins e segmentos como linhas horizontais grossas
coloridas por estado.

Por que a razao log2 do CNVkit NAO serve aqui
---------------------------------------------
O plot anterior
(`plot_cnv_scatter_colored.py`) fixa o eixo em -5..+5, seguindo o corte
automatico do cnvkit, e com isso **11 das 13 delecoes homozigoticas validadas
ficam invisiveis**: seus segmentos tem log2 entre -16 e -23 e sao desenhados
fora da tela: uma delecao homozigotica podia aparecer CINZA, neutra
-- o oposto do resultado, que e a delecao bialelica documentada no painel.

Aumentar o eixo nao resolve. Com referencia plana, o log2 de uma regiao sem
cobertura nao e uma medida: zero copias da log2 = -infinito, e o valor
reportado (-23, -25) e so o piso numerico do metodo. As duas delecoes que
"apareciam" estavam igualmente erradas, como perda leve, porque o segmento do
.cns e muito maior que o evento focal e dilui o sinal.

A escala usada aqui
-------------------
`log2((observado + 1) / (esperado + 1))`, por bin, onde o esperado vem das
estatisticas da coorte (scripts/cnv_bin_statistics.py). A pseudocontagem e a
correcao padrao para razao com denominador ou numerador perto de zero: torna a
metrica definida em zero e limitada por baixo, sem alterar o regime em que o
log2 ja era valido.

Comportamento medido:


Ou seja: do lado do ganho a escala concorda com a familiar, e do lado da perda
profunda ela traz o sinal de volta para dentro do grafico, onde a nuvem de
bins aparece deslocada para baixo, agrupada e separada do ruido.

Classificacao dos segmentos
---------------------------
Vem de `cnv_focal/{amostra}_focal_events.tsv`, isto e, da razao de cobertura
-- nao de um corte de log2. Perda homozigotica e amplificacao de alto nivel
sao desenhadas; perda heterozigotica e ganho de baixo nivel nao, porque seus
limiares nao foram calibrados contra verdade independente.
"""
import argparse
import csv
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

GAIN_COLOR = '#1f77b4'
LOSS_COLOR = '#d62728'
NEUTRAL_COLOR = '#888888'
PSEUDO = 1.0


def carregar_estatisticas(caminho, chrom):
    est = {}
    with open(caminho) as fh:
        cab = fh.readline().rstrip('\n').split('\t')
        i = {c: n for n, c in enumerate(cab)}
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            if c[i['chrom']] != chrom:
                continue
            est[(int(c[i['start']]), int(c[i['end']]))] = (
                c[i['tipo_de_bin']], float(c[i['profundidade_mediana_absoluta']]))
    return est


def carregar_bins(cnr, chrom, est, cobertura_minima):
    pontos = []
    with open(cnr) as fh:
        next(fh)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            if c[0] != chrom:
                continue
            a, b = int(c[1]), int(c[2])
            e = est.get((a, b))
            if not e or e[0] != 'alvo' or e[1] < cobertura_minima:
                continue
            try:
                obs = float(c[4])
            except ValueError:
                continue
            pontos.append(((a + b) / 2 / 1e6,
                           math.log2((obs + PSEUDO) / (e[1] + PSEUDO))))
    return pontos


def carregar_eventos(caminho, chrom, min_janelas, razao_ganho_minima):
    eventos = []
    with open(caminho, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['chrom'] != chrom or r.get('qc_da_amostra') == 'insuficiente':
                continue
            n = int(r['n_bins'])
            if n < min_janelas or int(r['n_bins_raros_na_coorte']) != n:
                continue
            try:
                razao = float(r['razao_de_cobertura'])
            except ValueError:
                continue
            classe = r['classe_do_evento']
            if classe == 'perda_homozigotica':
                tipo = 'loss'
            elif classe == 'ganho' and razao >= razao_ganho_minima:
                tipo = 'gain'
            else:
                continue
            try:
                esperada = float(r['cobertura_esperada_media_x'])
                observada = float(r['cobertura_observada_media_x'])
            except ValueError:
                continue
            eventos.append({
                'start': int(r['start']) / 1e6, 'end': int(r['end']) / 1e6,
                'tipo': tipo, 'n_bins': n,
                'y': math.log2((observada + PSEUDO) / (esperada + PSEUDO)),
                'razao': razao, 'obs': observada, 'esp': esperada,
                'genes': [g for g in r['genes'].split(';') if g and g != '-'],
            })
    return eventos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cnr', required=True)
    ap.add_argument('--eventos', required=True, help='{amostra}_focal_events.tsv')
    ap.add_argument('--bin-statistics', required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--titulo', default=None)
    ap.add_argument('--saida', required=True)
    ap.add_argument('--destacar-gene', action='append', default=[],
                    help='nome do gene a rotular; pode repetir')
    ap.add_argument('--min-janelas', type=int, default=5)
    ap.add_argument('--razao-ganho-minima', type=float, default=4.0)
    ap.add_argument('--cobertura-esperada-minima', type=float, default=5.0)
    args = ap.parse_args()

    est = carregar_estatisticas(args.bin_statistics, args.chrom)
    if not est:
        sys.exit(f'ERRO: nenhum bin de {args.chrom} nas estatísticas')
    pontos = carregar_bins(args.cnr, args.chrom, est, args.cobertura_esperada_minima)
    if not pontos:
        sys.exit(f'ERRO: nenhum bin utilizável em {args.chrom}')
    eventos = carregar_eventos(args.eventos, args.chrom, args.min_janelas,
                               args.razao_ganho_minima)

    fig, ax = plt.subplots(figsize=(16, 6))
    xs = [p[0] for p in pontos]
    ys = [p[1] for p in pontos]
    ax.scatter(xs, ys, s=3, c='#555555', alpha=0.28, linewidths=0, zorder=2)
    ax.axhline(0, color='black', linewidth=0.8, zorder=3)

    # Desenha todos os segmentos, mas ROTULA de forma esparsa: quando varios
    # eventos caem proximos -- dezenas num trecho de poucos Mb
    # -- os rotulos se empilham e a figura fica ilegivel. Rotula-se o evento de
    # maior numero de janelas em cada agrupamento, e o agrupamento ganha a
    # contagem dos demais. Nenhum evento deixa de ser DESENHADO.
    extensao = max(xs) - min(xs)
    for e in eventos:
        cor = GAIN_COLOR if e['tipo'] == 'gain' else LOSS_COLOR
        largura_desenho = max(e['end'] - e['start'], extensao * 0.004)
        ax.plot([e['start'], e['start'] + largura_desenho], [e['y'], e['y']],
                color=cor, linewidth=5, solid_capstyle='butt', zorder=5)

    destaque = {x.upper() for x in args.destacar_gene}
    candidatos = [e for e in eventos
                  if not destaque or any(g.upper() in destaque for g in e['genes'])]
    candidatos.sort(key=lambda e: -e['n_bins'])
    rotulados, ocupado = [], []
    for e in candidatos:
        centro = (e['start'] + e['end']) / 2
        if any(abs(centro - c) < extensao * 0.07 for c in ocupado):
            vizinho = min(rotulados, key=lambda r: abs(
                (r['start'] + r['end']) / 2 - centro))
            vizinho['_agrupados'] = vizinho.get('_agrupados', 0) + 1
            continue
        ocupado.append(centro)
        rotulados.append(e)

    for e in rotulados:
        cor = GAIN_COLOR if e['tipo'] == 'gain' else LOSS_COLOR
        largura_desenho = max(e['end'] - e['start'], extensao * 0.004)
        alvo = [g for g in e['genes'] if not destaque or g.upper() in destaque]
        rotulo = alvo[0] if alvo else (e['genes'][0] if e['genes'] else None)
        if not rotulo:
            continue
        if len(e['genes']) > 1 and not destaque:
            rotulo += f" (+{len(e['genes']) - 1})"
        desc = (f"{rotulo}\n{e['obs']:.2f}× obs / {e['esp']:.1f}× esp"
                f"  ·  {e['n_bins']} janelas")
        if e.get('_agrupados'):
            desc += f"\n(+{e['_agrupados']} evento(s) vizinho(s))"
        ax.annotate(desc, xy=(e['start'] + largura_desenho / 2, e['y']),
                    xytext=(0, -34 if e['tipo'] == 'loss' else 26),
                    textcoords='offset points', ha='center',
                    va='top' if e['tipo'] == 'loss' else 'bottom',
                    fontsize=9, fontweight='bold', color=cor, zorder=6,
                    arrowprops=dict(arrowstyle='-', color=cor, linewidth=1.1))

    inferior = min(min(ys), min([e['y'] for e in eventos], default=0)) - 1.0
    superior = max(max(ys), max([e['y'] for e in eventos], default=0)) + 1.0
    ax.set_ylim(min(inferior, -3.0), max(superior, 3.0))
    ax.set_xlim(min(xs) - 1, max(xs) + 1)
    ax.set_xlabel('Posição (Mb)', fontsize=12)
    ax.set_ylabel('Razão de cópias — log2((obs+1)/(esp+1))', fontsize=12)
    ax.set_title(args.titulo or f'{args.chrom}', fontsize=14, fontweight='bold')
    ax.legend(handles=[
        mpatches.Patch(color=LOSS_COLOR, label='Perda homozigótica'),
        mpatches.Patch(color=GAIN_COLOR, label='Amplificação de alto nível'),
        mpatches.Patch(color=NEUTRAL_COLOR, label='Bins (janelas de medição)'),
    ], loc='lower right', fontsize=9, framealpha=0.95)
    ax.grid(axis='y', alpha=0.18, zorder=0)
    for lado in ('top', 'right'):
        ax.spines[lado].set_visible(False)
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.saida)) or '.', exist_ok=True)
    plt.savefig(args.saida, dpi=200)
    plt.close(fig)
    print(f'{args.chrom}: {len(pontos)} bins, {len(eventos)} evento(s) -> {args.saida}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
