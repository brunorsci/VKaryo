#!/usr/bin/env python3
"""Build the self-contained test dataset described in docs/example.md.

Writes a minimal reference (one synthetic chromosome), paired FASTQ for the
exome and the transcriptome branches, and the calls the run is expected to
produce. Two events are planted, so the answer is known by construction and
the example doubles as a regression check:

  * a homozygous deletion -- an interval left with zero reads;
  * a gene fusion -- reads whose two halves come from two distant genes.

Everything is generated from a fixed seed, so two runs give byte-identical
files. Nothing is downloaded.
"""

import argparse
import gzip
import os
import random

CHROM = 'chr1'
TAMANHO = 3_000_000          # cromossomo sintetico de 3 Mb
LEITURA = 100                # comprimento da leitura
INSERTO = 300                # tamanho do fragmento
COBERTURA = 30               # profundidade alvo nas regioes-alvo

# Genes do GTF minimo: (nome, inicio, fim)
GENES = [
    ('GENEA', 200_000, 230_000),
    ('GENEB', 700_000, 730_000),
    ('DELGENE', 1_500_000, 1_530_000),   # dentro da delecao plantada
    ('GENEC', 2_300_000, 2_330_000),
]
DELECAO = (1_480_000, 1_580_000)         # ~100 kb sem cobertura
FUSAO = ('GENEA', 'GENEB')               # par que as leituras quimericas unem


def genoma(rng):
    return ''.join(rng.choice('ACGT') for _ in range(TAMANHO))


def complemento_reverso(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def escrever_fasta(caminho, seq):
    with open(caminho, 'w') as fh:
        fh.write(f'>{CHROM}\n')
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + '\n')


def escrever_gtf(caminho):
    with open(caminho, 'w') as fh:
        for nome, ini, fim in GENES:
            attr = f'gene_id "{nome}"; gene_name "{nome}"; transcript_id "{nome}.1";'
            fh.write(f'{CHROM}\ttest\tgene\t{ini}\t{fim}\t.\t+\t.\t{attr}\n')
            fh.write(f'{CHROM}\ttest\ttranscript\t{ini}\t{fim}\t.\t+\t.\t{attr}\n')
            fh.write(f'{CHROM}\ttest\texon\t{ini}\t{fim}\t.\t+\t.\t{attr}\n')


def escrever_bed(caminho):
    """Alvos de captura: os genes, mais janelas espalhadas para o CNV ter
    pontos de medicao fora deles."""
    with open(caminho, 'w') as fh:
        for nome, ini, fim in GENES:
            fh.write(f'{CHROM}\t{ini}\t{fim}\t{nome}\n')
        for i in range(100_000, TAMANHO - 100_000, 50_000):
            fh.write(f'{CHROM}\t{i}\t{i + 2_000}\twindow_{i}\n')


def na_delecao(pos):
    return DELECAO[0] <= pos <= DELECAO[1]


def par(rng, seq, pos, nome, qual='I'):
    r1 = seq[pos:pos + LEITURA]
    inicio2 = pos + INSERTO - LEITURA
    r2 = complemento_reverso(seq[inicio2:inicio2 + LEITURA])
    q = qual * LEITURA
    return (f'@{nome}/1\n{r1}\n+\n{q}\n', f'@{nome}/2\n{r2}\n+\n{q}\n')


def gerar(seq, rng, alvos, com_fusao, com_delecao):
    """Devolve (linhas_r1, linhas_r2).

    So a amostra com a delecao plantada tem seus fragmentos suprimidos no
    intervalo: a deteccao focal compara a cobertura observada com a esperada
    derivada da coorte, de modo que e preciso haver amostra COBERTA naquele
    lócus para que a esperada nao seja zero. Com todas as amostras deletadas o
    evento desapareceria por construcao.
    """
    r1s, r2s = [], []
    n = 0
    for ini, fim in alvos:
        quantos = max(1, (fim - ini) * COBERTURA // (LEITURA * 2))
        for _ in range(quantos):
            pos = rng.randint(ini, max(ini, fim - INSERTO))
            if com_delecao and (na_delecao(pos) or na_delecao(pos + INSERTO)):
                continue
            a, b = par(rng, seq, pos, f'read{n}')
            r1s.append(a); r2s.append(b); n += 1
    if com_fusao:
        # leituras quimericas: metade em GENEA, metade em GENEB
        a_ini = dict((g[0], g[1]) for g in GENES)[FUSAO[0]]
        b_ini = dict((g[0], g[1]) for g in GENES)[FUSAO[1]]
        for i in range(60):
            esq = seq[a_ini + 10_000 + i: a_ini + 10_000 + i + LEITURA // 2]
            dir_ = seq[b_ini + 10_000 + i: b_ini + 10_000 + i + LEITURA // 2]
            q = 'I' * LEITURA
            r1s.append(f'@fusion{i}/1\n{esq + dir_}\n+\n{q}\n')
            r2s.append(f'@fusion{i}/2\n{complemento_reverso(esq + dir_)}\n+\n{q}\n')
    return r1s, r2s


def escrever_fastq(caminho, linhas):
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with gzip.open(caminho, 'wt') as fh:
        fh.writelines(linhas)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='test', help='diretorio de saida')
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    ref = os.path.join(args.outdir, 'reference')
    os.makedirs(ref, exist_ok=True)

    seq = genoma(rng)
    escrever_fasta(os.path.join(ref, 'mini.fa'), seq)
    escrever_gtf(os.path.join(ref, 'mini.gtf'))
    escrever_bed(os.path.join(ref, 'targets.bed'))

    alvos = [(g[1], g[2]) for g in GENES]
    alvos += [(i, i + 2_000) for i in range(100_000, TAMANHO - 100_000, 50_000)]

    # sample1 carrega os dois eventos; sample2 e o controle, sem nenhum.
    for amostra, fusao, delecao in (('sample1', True, True), ('sample2', False, False)):
        for ramo in ('wes', 'rna'):
            r1, r2 = gerar(seq, random.Random(args.seed), alvos,
                           fusao and ramo == 'rna', delecao)
            escrever_fastq(os.path.join(args.outdir, 'reads', ramo, f'{amostra}_1.fastq.gz'), r1)
            escrever_fastq(os.path.join(args.outdir, 'reads', ramo, f'{amostra}_2.fastq.gz'), r2)

    esp = os.path.join(args.outdir, 'expected')
    os.makedirs(esp, exist_ok=True)
    with open(os.path.join(esp, 'expected_calls.tsv'), 'w') as fh:
        fh.write('sample\ttype\tchrom\tstart\tend\tdetail\n')
        fh.write(f'sample1\tfusion\t{CHROM}\t-\t-\t{FUSAO[0]}--{FUSAO[1]}\n')
        fh.write(f'sample1\thomozygous_loss\t{CHROM}\t{DELECAO[0]}\t{DELECAO[1]}\tDELGENE\n')

    print(f'test dataset written to {args.outdir}/')
    print('next: build the mini indexes, as described in docs/example.md')


if __name__ == '__main__':
    main()
