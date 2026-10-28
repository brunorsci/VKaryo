#!/usr/bin/env python3
"""Detecta eventos focais de número de cópias apoiados em evidência de
cobertura, usando as estatísticas por bin da coorte (cnv_bin_statistics.py).

Problema que resolve
--------------------
Com referência plana, a razão log2 do CNVkit perde validade quantitativa
quando a cobertura do segmento tende a zero: o valor reportado é o piso
numérico do método, não a magnitude da perda. Um lócus com cobertura
residual de menos de 1x contra dezenas de x
esperado -- razão que corresponde a log2 = -5,0 --, e o CNVkit reporta -23,4.
Erro de 18 unidades de log2. No lado dos ganhos não há esse problema.

Consequência prática: deleções homozigóticas, que são justamente o achado de
maior peso clínico, não podem ser chamadas nem graduadas pelo log2 nativo.

O que este script faz
---------------------
Avalia cada bin pela razão entre a cobertura OBSERVADA na amostra e a
cobertura ESPERADA para aquele bin, vinda das estatísticas da coorte. Agrupa
bins contíguos que satisfazem o critério em eventos, e reporta cada evento com
a evidência que o sustenta.

Classificação (descritiva, derivada da razão de cobertura -- não é rótulo de
confiança inventado):
  perda_homozigotica   razão < limiar_homozigotico (padrão 0,05)
  perda_heterozigotica razão entre o limiar homozigótico e limiar_perda (0,70)
  ganho                razão > limiar_ganho (1,50)

Colunas de evidência por evento, para que a avaliação não dependa de um número
só: número de bins, cobertura observada e esperada, razão, log2 derivado da
razão, e `n_bins_raros` -- quantos dos bins do evento caem em menos de
`max_frequencia_coorte` das amostras da coorte de referência. Um evento cujos
bins caem em muitas amostras pode ser artefato de captura OU alteração
recorrente de verdade (deleção de CDKN2A, recombinação V(D)J no lócus
TRA/TRD): a coluna informa, não decide.

Guarda de qualidade da amostra
------------------------------
A raridade do bin na coorte NÃO protege contra uma amostra globalmente mal
sequenciada: se só ela tem cobertura baixa, suas quedas são raras por
construção. Uma amostra assim pode ter profundidade
mediana de 2,8x contra ~27x na coorte -- 51.436 bins em queda, 125 vezes o
segundo colocado -- e produzia 5 "deleções homozigóticas" em genes-driver
(MECOM, IKZF1, TET2, CUX1, BTK) que são artefato de cobertura.

Por isso cada execução calcula a razão global da amostra (cobertura observada
sobre esperada, em toda a extensão avaliável) e a reporta na coluna
`razao_global_da_amostra`. Abaixo de `--razao-global-minima` (padrão 0,30) a
amostra recebe `qc_da_amostra = insuficiente` em todas as suas linhas, e seus
eventos não devem ser lidos como biologia sem reavaliação.

Genérico por construção: não usa lista de genes, painel de doença nem
pressuposto de neoplasia. Serve a câncer e a não câncer.
"""
import argparse
import csv
import math
import os
import sys

# Regioes em que a ausencia de cobertura e biologia esperada e NAO achado de
# numero de copias: os locus de imunoglobulina e de receptor de celula T
# sofrem recombinacao V(D)J, que excisa fisicamente segmentos em qualquer
# amostra linfoide. Toda linhagem B deste painel mostra IGH em 14q32 vazio --
# isso e linfocito normal, nao delecao. Marcadas e MANTIDAS na saida, para
# que a chamada possa ser auditada, em vez de descartadas em silencio.
# Coordenadas GRCh38.
LOCI_VDJ = [
    ('chr14', 105_550_000, 106_900_000, 'IGH'),
    ('chr22', 22_020_000, 22_930_000, 'IGL'),
    ('chr2', 88_800_000, 90_300_000, 'IGK'),
    ('chr14', 21_620_000, 22_560_000, 'TRA/TRD'),
    ('chr7', 142_290_000, 142_820_000, 'TRB'),
    ('chr7', 38_240_000, 38_380_000, 'TRG'),
]


def locus_vdj(chrom, inicio, fim):
    for c, a, b, nome in LOCI_VDJ:
        if chrom == c and inicio < b and fim > a:
            return f'recombinacao_VDJ_{nome}'
    return '-'


def carregar_estatisticas(caminho):
    est = {}
    with open(caminho) as fh:
        cab = fh.readline().rstrip('\n').split('\t')
        i = {c: n for n, c in enumerate(cab)}
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            est[(c[i['chrom']], c[i['start']])] = (
                c[i['end']],
                c[i['tipo_de_bin']],
                float(c[i['profundidade_mediana_absoluta']]),
                int(c[i['n_amostras_com_queda']]),
                int(c[i['n_amostras']]),
            )
    return est


def classificar(razao, lim_homo, lim_perda, lim_ganho):
    if razao < lim_homo:
        return 'perda_homozigotica'
    if razao < lim_perda:
        return 'perda_heterozigotica'
    if razao > lim_ganho:
        return 'ganho'
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cnr', required=True, help='.cnr da amostra')
    ap.add_argument('--bin-statistics', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--limiar-homozigotico', type=float, default=0.05)
    ap.add_argument('--limiar-perda', type=float, default=0.70)
    ap.add_argument('--limiar-ganho', type=float, default=1.50)
    ap.add_argument('--cobertura-esperada-minima', type=float, default=5.0,
                    help='bins com cobertura esperada abaixo disto são ignorados: '
                         'a ausência de leituras ali não é informativa')
    ap.add_argument('--min-bins', type=int, default=3,
                    help='mínimo de bins contíguos para constituir um evento')
    ap.add_argument('--max-distancia', type=int, default=100_000,
                    help='distancia maxima entre bins do mesmo evento, em pb. '
                         'NAO se aplica quando o intervalo nao contem janela '
                         'avaliavel alguma -- ver --nao-quebrar-em-vao-sem-sondas')
    ap.add_argument('--quebrar-em-vao-sem-sondas', action='store_true',
                    help='volta ao comportamento antigo, quebrando o evento em '
                         'qualquer intervalo maior que --max-distancia')
    ap.add_argument('--max-frequencia-coorte', type=float, default=0.05,
                    help='fração da coorte abaixo da qual a queda do bin é tida por rara')
    ap.add_argument('--razao-global-minima', type=float, default=0.30,
                    help='abaixo disto a amostra e marcada qc_da_amostra=insuficiente')
    ap.add_argument('--incluir-antialvo', action='store_true',
                    help='por padrão só bins de alvo entram: em exoma os bins de '
                         'antialvo têm cobertura esparsa por construção')
    args = ap.parse_args()

    est = carregar_estatisticas(args.bin_statistics)
    if not est:
        sys.exit('ERRO: estatísticas de bin vazias')
    n_coorte = next(iter(est.values()))[4]

    eventos = []
    atual = None
    descartados_esperada = 0
    soma_obs = soma_esp = 0.0
    n_avaliaveis = 0

    with open(args.cnr) as fh:
        next(fh)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            chrom, inicio, fim, gene = c[0], c[1], c[2], c[3]
            try:
                obs = float(c[4])
            except (ValueError, IndexError):
                continue
            e = est.get((chrom, inicio))
            if e is None:
                continue
            _fim, tipo, esperada, quedas, _n = e
            if tipo == 'contig_nao_cromossomico':
                continue
            if tipo == 'antialvo' and not args.incluir_antialvo:
                continue
            if esperada < args.cobertura_esperada_minima:
                descartados_esperada += 1
                continue

            soma_obs += obs
            soma_esp += esperada
            n_avaliaveis += 1

            razao = obs / esperada
            classe = classificar(razao, args.limiar_homozigotico,
                                 args.limiar_perda, args.limiar_ganho)
            if classe is None:
                # janela avaliavel que NAO pertence a classe: se um evento
                # estava aberto, o vao a partir daqui contem sonda, logo uma
                # eventual retomada adiante e um evento NOVO
                if atual is not None:
                    atual['_janela_no_vao'] = True
                continue

            raro = quedas < args.max_frequencia_coorte * n_coorte
            bin_dado = (chrom, int(inicio), int(fim), gene, obs, esperada, raro)
            # Um intervalo grande entre duas janelas do mesmo evento NAO e
            # prova de interrupcao quando nao ha sonda nenhuma no meio:
            # ausencia de sonda nao e ausencia de delecao. Medido em
            # Caso tipico: uma delecao homozigotica que cobre
            # 1.286 kb continuos de janelas zeradas, mas dois vaos sem sonda
            # (417 kb e 303 kb) a fragmentavam em tres eventos de 346, 214 e
            # 5 kb -- subestimando a extensao em 5x e separando CDKN2A de
            # MTAP e do agrupamento de IFNA. Agora so se quebra num vao que
            # de fato contenha janela avaliavel.
            contiguo = False
            if atual and atual['chrom'] == chrom and atual['classe'] == classe:
                distancia = int(inicio) - atual['fim']
                if distancia <= args.max_distancia:
                    contiguo = True
                elif not args.quebrar_em_vao_sem_sondas:
                    contiguo = not atual.get('_janela_no_vao', False)
            if contiguo:
                atual['fim'] = int(fim)
                atual['bins'].append(bin_dado)
                atual['_janela_no_vao'] = False
            else:
                atual = {'chrom': chrom, 'inicio': int(inicio), 'fim': int(fim),
                         'classe': classe, 'bins': [bin_dado],
                         '_janela_no_vao': False}
                eventos.append(atual)

    razao_global = (soma_obs / soma_esp) if soma_esp > 0 else 0.0
    qc_amostra = ('suficiente' if razao_global >= args.razao_global_minima
                  else 'insuficiente')

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
    with open(args.out, 'w', newline='') as fo:
        w = csv.writer(fo, delimiter='\t', lineterminator='\n')
        w.writerow(['sample', 'chrom', 'start', 'end', 'tamanho_pb', 'classe_do_evento',
                    'n_bins', 'cobertura_observada_media_x', 'cobertura_esperada_media_x',
                    'razao_de_cobertura', 'log2_por_razao_de_cobertura',
                    'n_bins_raros_na_coorte', 'n_amostras_da_coorte',
                    'razao_global_da_amostra', 'qc_da_amostra',
                    'fracao_janelas_com_cobertura_na_coorte', 'marcacao', 'genes'])
        n_escritos = 0
        for ev in eventos:
            bins = ev['bins']
            if len(bins) < args.min_bins:
                continue
            obs = sum(b[4] for b in bins) / len(bins)
            esp = sum(b[5] for b in bins) / len(bins)
            razao = obs / esp if esp > 0 else 0.0
            log2 = f'{math.log2(razao):.4f}' if razao > 0 else 'sem_cobertura'
            genes = sorted({g for b in bins for g in b[3].split(',')
                            if g and g not in ('Antitarget', '-', '.')})
            w.writerow([args.sample, ev['chrom'], ev['inicio'], ev['fim'],
                        ev['fim'] - ev['inicio'], ev['classe'], len(bins),
                        f'{obs:.4f}', f'{esp:.4f}', f'{razao:.6f}', log2,
                        sum(1 for b in bins if b[6]), n_coorte,
                        f'{razao_global:.4f}', qc_amostra,
                        f'{sum(1 for b in bins if b[5] >= 10.0) / len(bins):.4f}',
                        locus_vdj(ev['chrom'], ev['inicio'], ev['fim']),
                        ';'.join(genes) if genes else '-'])
            n_escritos += 1

    print(f'{args.sample}: {n_escritos} eventos focais '
          f'(de {len(eventos)} grupos; mínimo {args.min_bins} bins) -> {args.out}',
          file=sys.stderr)
    print(f'  {descartados_esperada} bins ignorados por cobertura esperada '
          f'< {args.cobertura_esperada_minima:g}x', file=sys.stderr)
    print(f'  razao global da amostra: {razao_global:.4f} '
          f'({n_avaliaveis} bins avaliaveis) -> qc_da_amostra={qc_amostra}',
          file=sys.stderr)
    if qc_amostra == 'insuficiente':
        print(f'  AVISO: cobertura global abaixo de {args.razao_global_minima:g} da '
              'esperada -- os eventos desta amostra podem ser artefato de cobertura '
              'e nao devem ser lidos como biologia sem reavaliacao', file=sys.stderr)


if __name__ == '__main__':
    main()
