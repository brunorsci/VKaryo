#!/usr/bin/env python3
"""Detecta desequilibrio de numero de copias em nivel de BRACO cromossomico,
a partir da razao de cobertura por janela de medicao.

Por que esta etapa existe
-------------------------
A deteccao de eventos focais (scripts/cnv_detect_focal_events.py) so reporta
perda homozigotica e ganho de alto nivel -- os regimes em que a evidencia de
cobertura e inequivoca. Eventos de braco sao, em sua maioria, HEMIZIGOTICOS:
perder uma de duas copias de um braco inteiro nao zera a cobertura, so a
reduz pela metade. Esses eventos apareciam na nuvem de bins do cariograma mas
NAO entravam em tabela nenhuma.

Exemplo tipico: uma amostra com Xq em torno de 29x e
Xp com 15,50x, contra 30,45x no chr2 -- ou seja, Xq com duas copias (0,96 do
autossomo) e Xp com uma (0,51). E uma perda de braco, del(Xp), invisivel para
a etapa focal. A comparacao e interna a amostra, portanto sem efeito de sexo.
Na mesma medida, CMK mostra desequilibrio no sentido oposto (Xq/Xp = 0,51),
enquanto amostras equilibradas ficam proximas de zero.

Metrica
-------
Por janela: `log2((observado + 1) / (esperado + 1))`, com o esperado vindo das
estatisticas da coorte. A pseudocontagem mantem a metrica definida quando a
cobertura tende a zero. Por braco: a MEDIANA dessas razoes, menos a mediana
autossomica da propria amostra -- essa subtracao remove deslocamento global da
amostra (profundidade, conteudo GC, pureza) e deixa so o desequilibrio
relativo entre bracos.

Criterio, alinhado a literatura
-------------------------------
Uma perda de uma copia em fundo diploide daria -1,0 e um ganho de uma copia
+0,58, em teoria. Na pratica o valor e atenuado: linhagens nao sao puras nem
necessariamente diploides, e a medida em exoma sofre compressao. O proprio
na pratica um chrX haploide mede entre -0,79 e -0,95, e nao -1,0.

A pratica publicada para chamada em nivel de braco combina DOIS criterios --
uma amplitude e uma fracao do braco:

- Shukla et al., Nat Commun 2020 (aneuploidias de braco em evolucao tumoral e
  resposta a droga), que inclui linhagens do GDSC e e portanto o desenho mais
  proximo deste: segmentos com log2 > 0,2 ou < -0,2, e o braco e chamado
  quando >= 0,9 da sua extensao esta alterada.
- Revisao de limiares em meningioma, Nat Commun 2025: a fracao exigida varia
  de 5% a 80% entre estudos publicados -- ou seja, nao ha consenso, e a
  escolha precisa ser declarada.
- Estudos de CNA em exoma costumam exigir ainda um numero minimo de exons ou
  segmentos por evento (p.ex. >= 6 exons com log2 medio alem de +-0,3).

Aqui os dois criterios sao aplicados e ambos saem na tabela: a MEDIANA do
desvio do braco (amplitude robusta, sem depender de segmentacao) e a FRACAO
de janelas do braco alem do limiar de amplitude. O braco so e chamado quando
os dois passam. A fracao padrao e 0,50, valor intermediario na faixa
publicada, e e ajustavel -- quem quiser o criterio estrito de Shukla et al.
passa 0,90.

RESSALVA DOS CROMOSSOMOS SEXUAIS: a cobertura esperada por janela vem da
mediana da coorte, que MISTURA amostras com uma e com duas copias de X. Para
chrX e chrY, portanto, o desvio em relacao a essa expectativa nao distingue
sexo de alteracao, e as linhas correspondentes saem marcadas na coluna
`ressalva`. O que permanece valido nesses cromossomos e a comparacao INTERNA
entre os bracos da mesma amostra -- uma razao Xq/Xp proxima de 1,9 indica duas copias
de Xq contra uma de Xp independentemente do sexo --, e por isso a saida traz
tambem a razao entre os dois bracos do mesmo cromossomo.

LIMITACAO A DECLARAR: estes limiares NAO foram calibrados contra um conjunto
de cariotipos anotados. Um painel de referencia costuma documentar citogenetica para parte das
linhagens, mas quase toda ela e de translocacao BALANCEADA -- que por
definicao nao tem assinatura de profundidade de leitura. Apenas uma linhagem
(RPMI-8402, del(1)(p32)) traz alteracao com assinatura de numero de copias.
Logo esta saida e descritiva e exige validacao ortogonal antes de uso clinico.
"""
import argparse
import csv
import math
import os
import statistics
import sys
from collections import defaultdict


def carregar_bracos(caminho):
    """Limites p/q por cromossomo, a partir do arquivo de bandas. O centromero
    e a banda com giemsa 'acen'."""
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
    bracos = {}
    for chrom, bs in bandas.items():
        bs.sort()
        acen = [b for b in bs if b[3] == 'acen']
        if acen:
            centro_ini = min(b[0] for b in acen)
            centro_fim = max(b[1] for b in acen)
        else:
            ps = [b for b in bs if b[2].startswith('p')]
            if not ps:
                continue
            centro_ini = centro_fim = max(b[1] for b in ps)
        inicio = min(b[0] for b in bs)
        fim = max(b[1] for b in bs)
        if centro_ini > inicio:
            bracos[f'{chrom}p'] = (chrom, inicio, centro_ini)
        bracos[f'{chrom}q'] = (chrom, centro_fim, fim)
    return bracos


def carregar_estatisticas(caminho, cobertura_minima):
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


def razoes_por_braco(cnr, est, bracos):
    por_braco = defaultdict(list)
    autossomicas = []
    indice = defaultdict(list)
    for nome, (chrom, a, b) in bracos.items():
        indice[chrom].append((a, b, nome))
    with open(cnr) as fh:
        next(fh)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            chrom = c[0]
            try:
                a = int(c[1]); obs = float(c[4])
            except (ValueError, IndexError):
                continue
            esperada = est.get((chrom, a))
            if esperada is None:
                continue
            r = math.log2((obs + 1.0) / (esperada + 1.0))
            for ini, fim, nome in indice.get(chrom, ()):
                if ini <= a < fim:
                    por_braco[nome].append(r)
                    if chrom not in ('chrX', 'chrY'):
                        autossomicas.append(r)
                    break
    return por_braco, autossomicas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cnr', required=True)
    ap.add_argument('--bin-statistics', required=True)
    ap.add_argument('--cytobands', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min-janelas', type=int, default=50)
    ap.add_argument('--cobertura-esperada-minima', type=float, default=5.0)
    ap.add_argument('--limiar-perda', type=float, default=-0.30,
                    help='amplitude mediana do braco para perda')
    ap.add_argument('--limiar-ganho', type=float, default=0.25,
                    help='amplitude mediana do braco para ganho')
    ap.add_argument('--limiar-amplitude-janela', type=float, default=0.20,
                    help='alem deste valor a janela conta como alterada '
                         '(Shukla et al., Nat Commun 2020, usa 0,20)')
    ap.add_argument('--fracao-minima-do-braco', type=float, default=0.50,
                    help='fracao das janelas do braco que precisa estar alterada; '
                         'a literatura varia de 0,05 a 0,90 -- 0,90 reproduz o '
                         'criterio estrito de Shukla et al.')
    args = ap.parse_args()

    bracos = carregar_bracos(args.cytobands)
    est = carregar_estatisticas(args.bin_statistics, args.cobertura_esperada_minima)
    if not bracos or not est:
        sys.exit('ERRO: bandas ou estatísticas vazias')
    por_braco, autossomicas = razoes_por_braco(args.cnr, est, bracos)
    if not autossomicas:
        sys.exit('ERRO: nenhuma janela autossômica utilizável')
    base = statistics.median(autossomicas)

    linhas = []
    for nome in sorted(por_braco, key=lambda n: (len(n), n)):
        rs = por_braco[nome]
        if len(rs) < args.min_janelas:
            continue
        mediana = statistics.median(rs)
        ajustada = mediana - base
        lim = args.limiar_amplitude_janela
        frac_perda = sum(1 for r in rs if (r - base) <= -lim) / len(rs)
        frac_ganho = sum(1 for r in rs if (r - base) >= lim) / len(rs)
        if ajustada <= args.limiar_perda and frac_perda >= args.fracao_minima_do_braco:
            classe, frac = 'perda_de_braco', frac_perda
        elif ajustada >= args.limiar_ganho and frac_ganho >= args.fracao_minima_do_braco:
            classe, frac = 'ganho_de_braco', frac_ganho
        elif ajustada <= args.limiar_perda or ajustada >= args.limiar_ganho:
            classe = 'amplitude_sem_extensao'
            frac = max(frac_perda, frac_ganho)
        else:
            classe, frac = 'equilibrado', max(frac_perda, frac_ganho)
        q = statistics.quantiles(rs, n=4) if len(rs) >= 4 else [mediana] * 3
        chrom = bracos[nome][0]
        sexual = chrom in ('chrX', 'chrY')
        ressalva = ('expectativa_da_coorte_mistura_sexos; usar a razao entre bracos'
                    if sexual else '-')
        if sexual:
            classe = f'{classe}_sob_ressalva' if classe != 'equilibrado' else classe
        linhas.append({
            'sample': args.sample, 'braco': nome,
            'chrom': chrom, 'start': bracos[nome][1], 'end': bracos[nome][2],
            'n_janelas': len(rs),
            'mediana_log2_cobertura': f'{mediana:.4f}',
            'mediana_autossomica_da_amostra': f'{base:.4f}',
            'desvio_ajustado': f'{ajustada:.4f}',
            'q1': f'{q[0]:.4f}', 'q3': f'{q[2]:.4f}',
            'fracao_janelas_alteradas': f'{frac:.4f}',
            'fracao_janelas_em_perda': f'{frac_perda:.4f}',
            'fracao_janelas_em_ganho': f'{frac_ganho:.4f}',
            'classe_do_braco': classe,
            'ressalva': ressalva,
        })

    # razao entre os dois bracos do mesmo cromossomo -- comparacao interna,
    # imune a mistura de sexos na coorte
    por_nome = {l['braco']: l for l in linhas}
    for l in linhas:
        irmao = l['braco'][:-1] + ('q' if l['braco'].endswith('p') else 'p')
        outro = por_nome.get(irmao)
        l['desvio_do_braco_irmao'] = (
            f"{float(l['mediana_log2_cobertura']) - float(outro['mediana_log2_cobertura']):.4f}"
            if outro else 'NA')

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
    with open(args.out, 'w', newline='') as fo:
        w = csv.DictWriter(fo, fieldnames=list(linhas[0].keys()), delimiter='\t',
                           lineterminator='\n')
        w.writeheader()
        w.writerows(linhas)
    n_alt = sum(1 for l in linhas if l['classe_do_braco'] in
                ('perda_de_braco', 'ganho_de_braco'))
    print(f'{args.sample}: {len(linhas)} braços avaliados, {n_alt} com desequilíbrio '
          f'-> {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
