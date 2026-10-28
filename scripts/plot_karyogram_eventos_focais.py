#!/usr/bin/env python3
"""Desenha o cariograma dos eventos focais de numero de copias validados,
a partir dos dados que o workflow ja produziu.

Origem da logica
----------------
Um painel
por cromossomo, barra cinza representando o cromossomo inteiro em escala, e
retangulos coloridos sobre as regioes alteradas -- azul para ganho, vermelho
para perda. Nao se reexecuta o CNVkit:
apenas a logica do plot e reaproveitada, sobre a tabela de eventos focais
`cnv_focal/{amostra}_focal_events.tsv` ja existente.

Criterio do que e desenhado
---------------------------
Por padrao, **so entram os eventos para os quais ha evidencia direta de
cobertura**, com no minimo `--min-janelas` janelas, **todas** com queda
infrequente no conjunto de referencia, e em amostra de controle de qualidade
suficiente:

- **perda homozigotica** -- razao de cobertura abaixo de 0,05, isto e, ausencia
  praticamente completa de leituras onde se esperava cobertura normal;
- **ganho com razao de cobertura >= `--razao-ganho-minima`** (padrao 4,0), ou
  seja, amplificacao de alto nivel.

Perda heterozigotica NAO e desenhada, e ganho de baixo nivel tambem nao. Os
limiares dessas duas classes (0,70 e 1,50) nao foram calibrados contra verdade
independente -- ver as Limitacoes do artigo. Desenha-los daria ao leitor a
impressao de validacao que nao existe: na primeira versao deste script, aceitar
todo ganho acima de 1,50 produzia de 389 a 1.258 retangulos por figura, quase
todos de baixa magnitude, soterrando os achados reais. Com o criterio atual, o
limiar de 4,0 recupera as amplificacoes de BCL2 verificadas contra a coorte em
apenas as amplificacoes de alto nivel, de razao proxima de dez ou mais.

Verificacao contra falso resultado
----------------------------------
O script grava, ao lado da figura, um `.tsv` com exatamente as regioes
desenhadas e a evidencia de cada uma. A figura e conferivel linha a linha
contra essa tabela -- se um retangulo aparece na imagem, ele esta na tabela, e
vice-versa. Nenhuma classificacao e recalculada aqui.
"""
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Constantes de plotagem
CHR_ORDER = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
GAIN_COLOR = '#1f77b4'
LOSS_COLOR = '#d62728'
NEUTRAL_COLOR = '#ecf0f1'

CHR_LENGTHS = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
}

# Realce minimo para que um evento focal de poucos kb seja visivel numa barra
# de centenas de Mb. Sem isso, uma delecao de 6,4 kb no chr6 teria largura de
# 0,004% do painel e simplesmente nao apareceria -- a figura mentiria por
# omissao. A tabela de acompanhamento preserva as coordenadas reais.
LARGURA_MINIMA_MB = 1.5


def carregar(caminho, min_janelas, exigir_raras, incluir_heterozigotica,
             razao_ganho_minima):
    eventos, descartados = [], 0
    with open(caminho, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            classe = r['classe_do_evento']
            try:
                razao = float(r['razao_de_cobertura'])
            except (ValueError, KeyError):
                razao = None
            if classe == 'perda_homozigotica':
                tipo = 'loss'
            elif classe == 'ganho':
                if razao is None or razao < razao_ganho_minima:
                    descartados += 1
                    continue
                tipo = 'gain'
            elif classe == 'perda_heterozigotica' and incluir_heterozigotica:
                tipo = 'loss'
            else:
                descartados += 1
                continue
            if r.get('qc_da_amostra') == 'insuficiente':
                descartados += 1
                continue
            n = int(r['n_bins'])
            if n < min_janelas:
                descartados += 1
                continue
            if exigir_raras and int(r['n_bins_raros_na_coorte']) != n:
                descartados += 1
                continue
            eventos.append({
                'chrom': r['chrom'], 'start': int(r['start']), 'end': int(r['end']),
                'tipo': tipo, 'classe': classe, 'n_bins': n,
                'n_raros': int(r['n_bins_raros_na_coorte']),
                'razao': r['razao_de_cobertura'],
                'log2': r['log2_por_razao_de_cobertura'],
                'cob_obs': r['cobertura_observada_media_x'],
                'cob_esp': r['cobertura_esperada_media_x'],
                'genes': r['genes'],
            })
    return eventos, descartados


def desenhar(eventos, amostra, saida, rotulo=None):
    """Um painel por cromossomo, todos na MESMA escala.

    Normalizar a escala e essencial para a figura nao mentir: se cada painel
    usasse a propria extensao, o eixo desenhado no rodape valeria so para o
    ultimo cromossomo, e o leitor leria todos na escala dele. Com escala comum,
    o comprimento de cada barra e proporcional ao comprimento real do
    cromossomo e o eixo vale para a figura inteira.
    """
    escala_max = max(CHR_LENGTHS.values()) / 1e6

    fig, axes = plt.subplots(nrows=24, ncols=1, figsize=(15, 18), sharex=True,
                             gridspec_kw={'hspace': 0.55})
    for idx, chrom in enumerate(CHR_ORDER):
        ax = axes[idx]
        comprimento = CHR_LENGTHS[chrom] / 1e6
        ax.add_patch(Rectangle((0, 0), comprimento, 1, facecolor=NEUTRAL_COLOR,
                               edgecolor='#555555', linewidth=0.6))

        do_crom = sorted([e for e in eventos if e['chrom'] == chrom],
                         key=lambda e: e['start'])
        ultimo_rotulo_fim = -1e9
        for e in do_crom:
            inicio = e['start'] / 1e6
            largura = max((e['end'] - e['start']) / 1e6, LARGURA_MINIMA_MB)
            if inicio + largura > comprimento:
                inicio = max(0.0, comprimento - largura)
            ax.add_patch(Rectangle((inicio, 0), largura, 1,
                                   facecolor=GAIN_COLOR if e['tipo'] == 'gain' else LOSS_COLOR,
                                   alpha=0.9, edgecolor='none'))
            genes = [g for g in e['genes'].split(';') if g and g != '-']
            if not genes:
                continue
            texto = genes[0] if len(genes) == 1 else f'{genes[0]} (+{len(genes) - 1})'
            centro = inicio + largura / 2
            # evita sobreposicao de rotulos vizinhos no mesmo cromossomo
            if centro - ultimo_rotulo_fim < escala_max * 0.055:
                continue
            ultimo_rotulo_fim = centro
            ax.annotate(texto, xy=(centro, 1), xytext=(0, 2),
                        textcoords='offset points', ha='center', va='bottom',
                        fontsize=6.5, color='#222222')

        ax.set_xlim(0, escala_max)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        for lado in ('top', 'right', 'left'):
            ax.spines[lado].set_visible(False)
        if idx < len(CHR_ORDER) - 1:
            ax.spines['bottom'].set_visible(False)
        ax.set_ylabel(chrom.replace('chr', ''), fontsize=9, fontweight='bold',
                      rotation=0, labelpad=16, va='center')

    axes[-1].set_xlabel('Posição no cromossomo (Mb) — escala comum a todos os painéis',
                        fontsize=11)

    n_perda = sum(1 for e in eventos if e['tipo'] == 'loss')
    n_ganho = sum(1 for e in eventos if e['tipo'] == 'gain')
    fig.suptitle(f'Cariograma de eventos focais — {rotulo or amostra}',
                 fontsize=15, fontweight='bold', y=0.985)
    fig.text(0.5, 0.963,
             f'{n_perda} perda(s) homozigótica(s) e {n_ganho} amplificação(ões) de alto '
             f'nível, com evidência direta de cobertura',
             ha='center', fontsize=10, color='#333333')
    fig.text(0.5, 0.948,
             f'Eventos menores que {LARGURA_MINIMA_MB:g} Mb são desenhados com essa largura '
             f'mínima para ficarem visíveis; as coordenadas reais estão na tabela que '
             f'acompanha a figura.',
             ha='center', fontsize=8, color='#666666', style='italic')
    fig.legend(handles=[mpatches.Patch(facecolor=LOSS_COLOR, label='Perda homozigótica'),
                        mpatches.Patch(facecolor=GAIN_COLOR, label='Amplificação de alto nível')],
               loc='lower center', bbox_to_anchor=(0.5, 0.005), ncol=2, fontsize=10,
               frameon=False)
    fig.subplots_adjust(top=0.935, bottom=0.055, left=0.07, right=0.985)
    plt.savefig(saida, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--eventos', required=True, help='{amostra}_focal_events.tsv')
    ap.add_argument('--amostra', required=True)
    ap.add_argument('--rotulo', default=None, help='nome para o titulo (ex.: linhagem)')
    ap.add_argument('--saida', required=True, help='arquivo da figura (.png ou .pdf)')
    ap.add_argument('--min-janelas', type=int, default=5)
    ap.add_argument('--sem-exigir-raras', action='store_true')
    ap.add_argument('--incluir-heterozigotica', action='store_true')
    ap.add_argument('--razao-ganho-minima', type=float, default=4.0,
                    help='ganho so e desenhado acima desta razao de cobertura; '
                         'o limiar de 1,50 que define a classe nao foi calibrado')
    args = ap.parse_args()

    eventos, descartados = carregar(args.eventos, args.min_janelas,
                                    not args.sem_exigir_raras,
                                    args.incluir_heterozigotica,
                                    args.razao_ganho_minima)
    if not eventos:
        print(f'{args.amostra}: nenhum evento satisfaz o critério '
              f'({descartados} descartados) -- figura não gerada', file=sys.stderr)
        return

    os.makedirs(os.path.dirname(os.path.abspath(args.saida)) or '.', exist_ok=True)
    desenhar(eventos, args.amostra, args.saida, args.rotulo)

    tabela = os.path.splitext(args.saida)[0] + '_regioes_desenhadas.tsv'
    with open(tabela, 'w', newline='') as fo:
        w = csv.writer(fo, delimiter='\t', lineterminator='\n')
        w.writerow(['amostra', 'chrom', 'start', 'end', 'tamanho_pb', 'classe_do_evento',
                    'n_janelas', 'n_janelas_raras', 'cobertura_observada_x',
                    'cobertura_esperada_x', 'razao_de_cobertura',
                    'log2_por_razao_de_cobertura', 'genes'])
        for e in sorted(eventos, key=lambda x: (CHR_ORDER.index(x['chrom'])
                                                if x['chrom'] in CHR_ORDER else 99,
                                                x['start'])):
            w.writerow([args.amostra, e['chrom'], e['start'], e['end'],
                        e['end'] - e['start'], e['classe'], e['n_bins'], e['n_raros'],
                        e['cob_obs'], e['cob_esp'], e['razao'], e['log2'], e['genes']])

    print(f'{args.amostra}: {len(eventos)} evento(s) desenhado(s) '
          f'({descartados} fora do critério)', file=sys.stderr)
    print(f'  figura  -> {args.saida}', file=sys.stderr)
    print(f'  regiões -> {tabela}', file=sys.stderr)


if __name__ == '__main__':
    main()
