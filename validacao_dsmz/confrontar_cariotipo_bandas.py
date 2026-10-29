"""Confronta um cariotipo em bandas com o dado de cobertura, alteracao a alteracao.

Entrada: a notacao do cariotipo (ISCN) e o arquivo .cnr do CNVkit da amostra.
Saida: para cada alteracao de dosagem relatada, um veredito em tres camadas com
o intervalo que o dado delimita.

Calibrado contra leitura manual independente: reproduziu 11 de 11 casos numa
linhagem de referencia, incluindo quatro controles negativos. A taxa de bloco
espurio -- a mesma regra aplicada a bracos que o cariotipo nao declara alterados
-- e medida por scripts/controle_bloco_espurio.py, e sem ela a taxa de
recuperacao nao e interpretavel.

O CRITERIO, E POR QUE NAO E O DA VERSAO ANTERIOR
Pontuar "o braco inteiro acompanha a notacao" mede quanto o bandeamento erra,
nao quanto o workflow acerta: quando o cariotipo escreve -9 e o material do 9
esta retido num derivativo -- coisa que o proprio cariotipo registra com der(9)
--, exigir que os dois bracos caiam e exigir que o dado reproduza uma imprecisao.

Cada alteracao de dosagem relatada e classificada em tres camadas:
  concordante               bloco no sentido esperado cobrindo a maior parte do braco
  concordante_com_refino    bloco contiguo com borda definida, em parte do braco
  nao_reproduzido           nenhum bloco no sentido esperado

O bloco e procurado em janelas de 0,25 Mb, porque a mediana do braco dilui
delecao intersticial: na KASUMI-1 o -13 e uma delecao de ~19 a ~41 Mb cuja
mediana de braco fica em -0,24, abaixo de qualquer limiar razoavel de braco.
"""
import collections
import io
import os
import re
import statistics as st

REF = os.environ.get('VKARYO_REFERENCE', 'reference')
OUT = os.environ.get('VKARYO_OUTPUT', 'results/output_vkaryo')


JANELA = 250_000          # resolucao da busca de blocos
LIMIAR = 0.20             # desvio minimo por janela para contar no bloco
BLOCO_MIN = 3_000_000     # extensao minima de um bloco, em pb
FRACAO_BRACO = 0.80       # a partir daqui o bloco e tido por "o braco inteiro"
FRACAO_MIN = 0.10         # abaixo disso o bloco e pequeno demais para o braco
MAGNITUDE_MIN = 0.30      # mediana do proprio bloco; separa evento real de ruido
DENSIDADE_MIN = 0.50      # fracao do bloco que precisa ter sonda: um bloco feito
                          # sobretudo de vao e interpolacao, nao medicao


def carregar_bandas():
    bracos = {}
    with io.open(f'{REF}/arriba/cytobands_hg38_GRCh38_v2.5.0.tsv', encoding='utf-8') as fh:
        next(fh)
        for linha in fh:
            c = linha.rstrip('\n').split('\t')
            if len(c) < 4:
                continue
            cr, ini, fim, nome = c[0].replace('chr', ''), int(c[1]), int(c[2]), c[3]
            if not nome:
                continue
            b = nome[0]
            k = (cr, b)
            if k not in bracos:
                bracos[k] = [ini, fim]
            else:
                bracos[k][0] = min(bracos[k][0], ini)
                bracos[k][1] = max(bracos[k][1], fim)
    return bracos


BRACOS = carregar_bandas()


def perfil(caminho):
    """Mediana de log2 por janela, ajustada pela mediana autossomica da amostra."""
    jan = collections.defaultdict(list)
    autos = []
    with io.open(caminho, encoding='utf-8', errors='replace') as fh:
        next(fh)
        for linha in fh:
            c = linha.split('\t')
            if len(c) < 6:
                continue
            cr = c[0].replace('chr', '')
            try:
                ini, v = int(c[1]), float(c[5])
            except ValueError:
                continue
            if v < -15:                      # piso numerico do CNVkit
                v = -5.0
            jan[(cr, ini // JANELA)].append(v)
            if cr.isdigit():
                autos.append(v)
    base = st.median(autos) if autos else 0.0
    perf = collections.defaultdict(dict)
    for (cr, k), vs in jan.items():
        if len(vs) >= 1:
            perf[cr][k] = st.median(vs) - base
    return perf, base


def maior_bloco(perf, cr, braco, sentido):
    """Maior corrida contigua de janelas no sentido esperado. (ini, fim, n_com_dado).

    Uma janela SEM dado nao quebra a corrida: o exoma cobre os alvos do kit, e um
    vao sem sonda nao e evidencia de que o bloco terminou. So quebra a corrida uma
    janela que TEM dado e esta fora do sentido esperado. E o mesmo criterio de
    cnv_detect_focal_events.py, onde ignorar isso partia uma delecao unica em tres.
    """
    lim = BRACOS.get((cr, braco))
    if lim is None or cr not in perf:
        return None
    k0, k1 = lim[0] // JANELA, lim[1] // JANELA
    sinal = -1 if sentido == 'perda' else 1
    melhor = atual = None

    def fecha(a, m):
        if a and (m is None or (a[1] - a[0]) > (m[1] - m[0])):
            return a
        return m

    for k in range(k0, k1 + 1):
        v = perf[cr].get(k)
        if v is None:                      # vao sem sonda: nao decide nada
            continue
        if sinal * v >= LIMIAR:            # sustenta o bloco
            if atual is None:
                atual = [k, k, 1]
            else:
                atual[1], atual[2] = k, atual[2] + 1
        elif sinal * v <= -LIMIAR:         # contradiz: so isso encerra o bloco
            melhor = fecha(atual, melhor)
            atual = None
        # janela neutra, entre os dois limiares, e absorvida: um tremor isolado
        # nao e evidencia de que o bloco terminou, e quebrar nele fragmenta um
        # ganho de cromossomo inteiro em dezenas de pedacos.
    melhor = fecha(atual, melhor)
    if melhor is None:
        return None
    return melhor[0] * JANELA, (melhor[1] + 1) * JANELA, melhor[2]


def classificar(perf, cr, braco, sentido):
    lim = BRACOS.get((cr, braco))
    if lim is None:
        return 'sem_referencia', '', 0.0
    b = maior_bloco(perf, cr, braco, sentido)
    if b is None:
        return 'nao_reproduzido', '', 0.0
    ini, fim, n_dado = b
    ext = fim - ini
    arm = max(1, lim[1] - lim[0])
    frac = ext / arm
    if ext < BLOCO_MIN or n_dado < 8 or frac < FRACAO_MIN:
        return 'nao_reproduzido', '', frac
    # Passar do limiar janela a janela nao basta: regioes de cobertura pobre --
    # pericentromericas, sobretudo -- produzem corridas longas de desvio fraco.
    # Exige-se que a mediana do PROPRIO bloco corresponda a uma alteracao de
    # dosagem de verdade, e nao a uma deriva.
    janelas = list(range(ini // JANELA, fim // JANELA))
    vals = [perf[cr][k] for k in janelas if k in perf[cr]]
    if not vals or abs(st.median(vals)) < MAGNITUDE_MIN:
        return 'nao_reproduzido', '', frac
    # Um vao sem sonda nao quebra o bloco, mas um bloco feito sobretudo de vao
    # nao e medicao: e a reta que o grafico desenha entre dois pontos distantes.
    # O caso tipico e a heterocromatina de 9q12, o maior bloco heterocromatico
    # do genoma, praticamente sem alvo de exoma.
    if len(vals) / max(1, len(janelas)) < DENSIDADE_MIN:
        return 'nao_reproduzido', '', frac
    intervalo = f'{ini/1e6:.1f}–{fim/1e6:.1f} Mb'
    if frac >= FRACAO_BRACO:
        return 'concordante', intervalo, frac
    return 'concordante_com_refino', intervalo, frac


# ---- leitura do cariotipo ---------------------------------------------------
ALVO = re.compile(r'(?<![\w)])([+-])(\d{1,2}|X|Y)(?![\w(])')
DEL = re.compile(r'\b(del|dup|add)\((\d{1,2}|X|Y)\)\(([pq])')
ISO = re.compile(r'\bi\((\d{1,2}|X|Y)\)\(([pq])')
MAR = re.compile(r'\+(\d*)mar')


def eventos_do_cariotipo(texto):
    """Devolve lista de (notacao, cromossomo, braco ou None, sentido)."""
    ev = []
    corpo = texto.split('Molec')[0]
    for m in ALVO.finditer(corpo):
        sinal, cr = m.group(1), m.group(2)
        if cr == '' or 'mar' in corpo[m.end():m.end() + 4]:
            continue
        ev.append((m.group(0), cr, None, 'perda' if sinal == '-' else 'ganho'))
    for m in DEL.finditer(corpo):
        tipo, cr, br = m.groups()
        ev.append((m.group(0) + ')', cr, br, 'perda' if tipo == 'del' else 'ganho'))
    for m in ISO.finditer(corpo):
        cr, br = m.groups()
        ev.append((f'i({cr})({br}10)', cr, br, 'ganho'))
        ev.append((f'i({cr})({br}10)', cr, 'q' if br == 'p' else 'p', 'perda'))
    return ev
