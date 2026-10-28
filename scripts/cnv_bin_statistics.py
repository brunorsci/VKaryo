#!/usr/bin/env python3
"""Deriva estatísticas por bin a partir de uma coorte de referência, para que
a confiabilidade de uma chamada de CNV possa ser julgada DEPOIS, em qualquer
amostra isolada, sem precisar de coorte no momento da análise.

O problema que isto resolve
---------------------------
Com referência "flat" (sem amostras normais), o CNVkit não dispõe de variância
por bin e portanto não distingue "este bin nunca é bem capturado" de "esta
amostra perdeu esta região". Duas consequências, ambas medidas no painel
observada na pratica:

1. **Magnitude do log2 errada no lado das deleções.** Em uma deleção homozigótica, a
   profundidade medida no lócus de TP53 é 0,74x contra mediana de 29,3x na
   coorte -- razão que corresponde a log2 = -5,31. O CNVkit reporta -23,38:
   erro de 18 unidades de log2, um fator de 2^18. No lado das amplificações
   não há esse problema: a razão nativa concorda com a estimada pela
   razão de cobertura, e RI-1/BCL2 reporta +3,341 contra +3,09.

2. **Artefatos indistinguíveis de biologia.** 1.394 segmentos com log2 < -5,
   dos quais 72% recorrem nas MESMAS coordenadas em 5 ou mais linhagens (um
   local em 43 das 100), mais 500 segmentos em contigs não cromossômicos.

Por que anotar em vez de excluir
--------------------------------
- `--drop-low-coverage` removeria também as deleções homozigóticas
  verdadeiras, que por definição têm cobertura nula -- apagaria a deleção de
  um lócus supressor tumoral com cobertura residual contra dezenas de x, achado
  real e relevante em leucemia.
- Uma referência agrupada construída a partir da coorte não é aplicável a um
  paciente único, e o workflow precisa funcionar com uma amostra só (decisão
  de projeto).
- Excluir bins recorrentemente baixos apagaria CNV recorrente verdadeiro --
  deleção de CDKN2A, por exemplo, é recorrente em leucemia por biologia, não
  por artefato. Recorrência NÃO é prova de artefato.

Logo: nada é descartado. As estatísticas acompanham o workflow como arquivo de
referência e viram colunas na tabela de resultado, de modo que o analista (ou
a camada clínica) julgue cada segmento com a informação à vista.

Saída
-----
TSV com uma linha por bin:

  chrom, start, end, tipo_de_bin (alvo|antialvo),
  profundidade_mediana_normalizada  -- mediana, na coorte, da profundidade do
      bin dividida pela profundidade mediana da própria amostra. ~1,0 = bin
      típico; << 1 = bin mal capturado em geral.
  profundidade_mediana_absoluta     -- em x, para recalcular razões
  n_amostras_com_queda              -- em quantas amostras da coorte o bin
      ficou abaixo de 10% da mediana da amostra. Alto = o bin falha com
      frequência, logo a ausência de leituras ali é pouco informativa.
  n_amostras
"""
import argparse
import glob
import os
import statistics
import sys

NAO_CROMOSSOMICO = ('_random', '_alt', 'chrUn_', '_fix', '_decoy')
LIMIAR_QUEDA = 0.10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cnr-dir', required=True)
    ap.add_argument('--padrao', default='*_markduplicates_recal.cnr')
    ap.add_argument('--min-amostras', type=int, default=20)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    arquivos = sorted(glob.glob(os.path.join(args.cnr_dir, args.padrao)))
    if len(arquivos) < args.min_amostras:
        sys.exit(f'ERRO: {len(arquivos)} amostras, mínimo {args.min_amostras}')
    print(f'{len(arquivos)} amostras na coorte de referência', file=sys.stderr)

    coords = tipos = None
    norm_por_bin = None
    abs_por_bin = None
    quedas = None

    for n_arq, caminho in enumerate(arquivos):
        chrom_l, ini_l, fim_l, tipo_l, prof_l = [], [], [], [], []
        with open(caminho) as fh:
            next(fh)
            for linha in fh:
                c = linha.rstrip('\n').split('\t')
                chrom_l.append(c[0]); ini_l.append(c[1]); fim_l.append(c[2])
                tipo_l.append('antialvo' if c[3] in ('Antitarget', '-', '') else 'alvo')
                try:
                    prof_l.append(float(c[4]))
                except (ValueError, IndexError):
                    prof_l.append(0.0)
        atual = list(zip(chrom_l, ini_l, fim_l))
        if coords is None:
            coords, tipos = atual, tipo_l
            norm_por_bin = [[] for _ in coords]
            abs_por_bin = [[] for _ in coords]
            quedas = [0] * len(coords)
        elif atual != coords:
            sys.exit(f'ERRO: bins de {os.path.basename(caminho)} diferem -- '
                     'as estatísticas só valem se a referência for a mesma')
        mediana = statistics.median([p for p in prof_l if p > 0]) or 1.0
        for i, p in enumerate(prof_l):
            r = p / mediana
            norm_por_bin[i].append(r)
            abs_por_bin[i].append(p)
            if r < LIMIAR_QUEDA:
                quedas[i] += 1
        if (n_arq + 1) % 20 == 0:
            print(f'  {n_arq + 1}/{len(arquivos)}', file=sys.stderr)

    n = len(arquivos)
    with open(args.out, 'w') as fo:
        fo.write('chrom\tstart\tend\ttipo_de_bin\tprofundidade_mediana_normalizada\t'
                 'profundidade_mediana_absoluta\tn_amostras_com_queda\tn_amostras\n')
        for i, (chrom, inicio, fim) in enumerate(coords):
            tipo = ('contig_nao_cromossomico'
                    if any(m in chrom for m in NAO_CROMOSSOMICO) else tipos[i])
            fo.write(f'{chrom}\t{inicio}\t{fim}\t{tipo}\t'
                     f'{statistics.median(norm_por_bin[i]):.6f}\t'
                     f'{statistics.median(abs_por_bin[i]):.4f}\t'
                     f'{quedas[i]}\t{n}\n')

    print(f'\n{len(coords)} bins escritos em {args.out}', file=sys.stderr)
    from collections import Counter
    c = Counter('contig_nao_cromossomico'
                if any(m in ch for m in NAO_CROMOSSOMICO) else tipos[i]
                for i, (ch, _, _) in enumerate(coords))
    print(f'  composição: {dict(c)}', file=sys.stderr)
    baixos = sum(1 for i in range(len(coords))
                 if statistics.median(norm_por_bin[i]) < 0.05)
    freq = sum(1 for q in quedas if q >= 0.05 * n)
    print(f'  bins com mediana normalizada < 0,05: {baixos} '
          f'({100 * baixos / len(coords):.2f}%)', file=sys.stderr)
    print(f'  bins com queda em >= 5% das amostras: {freq} '
          f'({100 * freq / len(coords):.2f}%)', file=sys.stderr)


if __name__ == '__main__':
    main()
