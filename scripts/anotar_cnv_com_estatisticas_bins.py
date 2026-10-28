#!/usr/bin/env python3
"""Acrescenta à tabela consolidada de CNV as colunas que permitem julgar cada
segmento, usando as estatísticas por bin da coorte de referência
(scripts/gerar_estatisticas_bins_referencia.py).

Funciona para UMA amostra isolada: as estatísticas vêm do arquivo de
referência, não da coorte em análise.

Colunas acrescentadas
---------------------
cobertura_esperada_x
    Soma da profundidade mediana absoluta dos bins de referência que o
    segmento cobre, ponderada pelo número de bins -- a profundidade que se
    esperaria no segmento numa amostra típica.
razao_de_cobertura
    profundidade observada no segmento (coluna `depth` do CNVkit) dividida
    pela esperada acima.
log2_por_razao_de_cobertura
    log2 dessa razão. É a estimativa a usar no lado das deleções, onde o log2
    nativo do CNVkit com referência flat é quantitativamente inválido (medido:
    erro de quase vinte unidades de log2). No lado das amplificações os dois
    concordam.
n_bins_de_referencia
    Quantos bins de referência o segmento cobre.
fracao_bins_mal_capturados
    Fração desses bins cuja mediana normalizada na coorte é < 0,05, isto é,
    bins que quase nunca têm cobertura em amostra nenhuma. Alto => a ausência
    de leituras no segmento não informa sobre a amostra.
mediana_amostras_com_queda
    Mediana, entre os bins do segmento, do número de amostras da coorte de
    referência em que aquele bin caiu abaixo de 10% da mediana. Alto => o bin
    falha com frequência.
    ATENÇÃO: valor alto NÃO prova artefato. Deleção recorrente verdadeira
    (CDKN2A em leucemia, lócus TRA/TRD por recombinação V(D)J) também produz
    queda recorrente. A coluna informa; não decide.
segmento_em_contig_nao_cromossomico
    sim/não. Scaffolds não localizados (`_random`, `_alt`, `chrUn_`) não devem
    constar de resultado clínico.
"""
import argparse
import bisect
import csv
import math
import os
import statistics
import sys
from collections import defaultdict

csv.field_size_limit(10 ** 9)
NAO_CROMOSSOMICO = ('_random', '_alt', 'chrUn_', '_fix', '_decoy')
LIMIAR_MAL_CAPTURADO = 0.05

NOVAS = ['cobertura_esperada_x', 'razao_de_cobertura', 'log2_por_razao_de_cobertura',
         'n_bins_de_referencia', 'fracao_bins_mal_capturados',
         'mediana_amostras_com_queda', 'segmento_em_contig_nao_cromossomico']


def carregar(caminho):
    """Bins por cromossomo, ordenados por início, para busca binária."""
    por_crom = defaultdict(list)
    with open(caminho) as fh:
        next(fh)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            por_crom[c[0]].append((int(c[1]), int(c[2]), float(c[4]),
                                   float(c[5]), int(c[6])))
    for ch in por_crom:
        por_crom[ch].sort()
    return por_crom, {ch: [b[0] for b in v] for ch, v in por_crom.items()}


def bins_do_segmento(por_crom, inicios, chrom, inicio, fim):
    lista = por_crom.get(chrom)
    if not lista:
        return []
    # primeiro bin cujo início pode sobrepor: recua um para pegar bins que
    # começam antes e terminam dentro
    i = max(0, bisect.bisect_left(inicios[chrom], inicio) - 1)
    out = []
    for a, b, norm, absd, quedas in lista[i:]:
        if a >= fim:
            break
        if b > inicio:
            out.append((norm, absd, quedas))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--entrada', required=True)
    ap.add_argument('--estatisticas', required=True)
    ap.add_argument('--saida', required=True)
    args = ap.parse_args()

    por_crom, inicios = carregar(args.estatisticas)
    print(f'estatísticas de {sum(len(v) for v in por_crom.values())} bins', file=sys.stderr)

    n = sem_bins = 0
    with open(args.entrada, newline='') as fi, open(args.saida, 'w', newline='') as fo:
        leitor = csv.reader(fi, delimiter='\t')
        header = next(leitor)
        idx = {c: i for i, c in enumerate(header)}
        for c in ('chromosome', 'start', 'end', 'depth'):
            if c not in idx:
                sys.exit(f'ERRO: coluna {c} ausente na entrada')
        escritor = csv.writer(fo, delimiter='\t', lineterminator='\n')
        escritor.writerow(header + NOVAS)

        for linha in leitor:
            chrom = linha[idx['chromosome']]
            try:
                inicio, fim = int(linha[idx['start']]), int(linha[idx['end']])
                obs = float(linha[idx['depth']])
            except (ValueError, IndexError):
                escritor.writerow(linha + ['NA'] * len(NOVAS)); n += 1; continue

            nao_crom = 'sim' if any(m in chrom for m in NAO_CROMOSSOMICO) else 'nao'
            bins = bins_do_segmento(por_crom, inicios, chrom, inicio, fim)
            if not bins:
                escritor.writerow(linha + ['NA', 'NA', 'NA', '0', 'NA', 'NA', nao_crom])
                n += 1; sem_bins += 1; continue

            esperada = statistics.median(b[1] for b in bins)
            razao = obs / esperada if esperada > 0 else None
            log2r = math.log2(razao) if (razao and razao > 0) else None
            mal = sum(1 for b in bins if b[0] < LIMIAR_MAL_CAPTURADO) / len(bins)
            med_quedas = statistics.median(b[2] for b in bins)

            escritor.writerow(linha + [
                f'{esperada:.4f}',
                f'{razao:.6f}' if razao is not None else 'NA',
                f'{log2r:.4f}' if log2r is not None else 'sem_cobertura',
                str(len(bins)),
                f'{mal:.4f}',
                f'{med_quedas:.1f}',
                nao_crom,
            ])
            n += 1

    print(f'{n} segmentos anotados ({sem_bins} sem bin de referência) -> {args.saida}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
