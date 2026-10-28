#!/usr/bin/env python3
"""
Concatena CNVs (CNVkit .cns) e expressão (STAR ReadsPerGene.out.tab) de
TODAS as amostras DNA+RNA pareadas no formato esperado pelo pacote R
CNVExpression (github.com/brunorsci/CNVExpression_v1, instalado e testado
neste servidor -- ver inst/extdata/example_cnv.tsv e example_expression.tsv).
Lê da pasta de staging (output_vkaryo/cnv_expression/input/, ver
scripts/stage_cnv_expression_input.py), não direto de cnvkit/star_aligned --
mantém um snapshot auditável do que entrou na análise.

Por que concatenar amostras: o pacote (via CNVRanger::cnvEQTL) calcula log
fold-change de expressão por gene ENTRE grupos de amostras com estados de
número de cópias diferentes (CN0/CN1/CN3/CN4) -- precisa de várias amostras
por estado pra ter desvio/média (min_samples=3 por padrão no pacote). Uma
amostra sozinha não permite esse cálculo -- é uma limitação conhecida desta
primeira versão da análise (ver nota no rule file sobre extensão futura por
amostra usando expressão média de controles/referência).

Formato de saída CNV (sample, Chrom, Start, End, value):
    sample  Chrom  Start  End      value
    ERR...  chr1   100000 200000   -0.85

Formato de saída expressão (genes x amostras, contagens brutas do STAR,
coluna 'unstranded' -- mesma convenção já usada em validate_cnv_expression.py):
    Ensembl_ID  ERR...  ERR...
    DDX11L1     0       12

NOTA: apesar do nome "Ensembl_ID" na primeira coluna (é o nome de coluna
padrão esperado pela função load_expression_data() do pacote), os
identificadores aqui são SÍMBOLOS de gene (DDX11L1, TP53, ...), não IDs
ENSG -- o GTF usado pelo STAR neste pipeline foi indexado com gene_name
como identificador, não gene_id. Por isso o script R downstream chama
cnv_expression_pipeline(..., id_type = "symbol"), não "ensembl".
"""
import argparse
import csv
import glob
import os
import sys


def write_cnv_table(staging_dir, out_path):
    n_segments = 0
    samples = []
    with open(out_path, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['sample', 'Chrom', 'Start', 'End', 'value'])
        for cns_path in sorted(glob.glob(os.path.join(staging_dir, '*_cnv.cns'))):
            dna_id = os.path.basename(cns_path).replace('_cnv.cns', '')
            with open(cns_path) as f:
                header = f.readline().rstrip('\n').split('\t')
                idx = {h: i for i, h in enumerate(header)}
                ci, si, ei, li = idx.get('chromosome'), idx.get('start'), idx.get('end'), idx.get('log2')
                if None in (ci, si, ei, li):
                    print(f'  [aviso] colunas inesperadas em {cns_path}', file=sys.stderr)
                    continue
                had_rows = False
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) <= max(ci, si, ei, li):
                        continue
                    w.writerow([dna_id, fields[ci], fields[si], fields[ei], fields[li]])
                    n_segments += 1
                    had_rows = True
                if had_rows:
                    samples.append(dna_id)
    return n_segments, samples


def write_expression_table(staging_dir, out_path):
    matrix = {}  # gene -> {sample: count}
    samples_included = []
    for expr_path in sorted(glob.glob(os.path.join(staging_dir, '*_expression.tab'))):
        dna_id = os.path.basename(expr_path).replace('_expression.tab', '')
        samples_included.append(dna_id)
        with open(expr_path) as f:
            for i, line in enumerate(f):
                if i < 4:  # N_unmapped/N_multimapping/N_noFeature/N_ambiguous
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 2:
                    continue
                gene, unstranded = fields[0], fields[1]
                matrix.setdefault(gene, {})[dna_id] = unstranded

    with open(out_path, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(['Ensembl_ID'] + samples_included)
        for gene, counts in matrix.items():
            w.writerow([gene] + [counts.get(s, '0') for s in samples_included])

    return len(matrix), samples_included


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--staging-dir', required=True, help='output_vkaryo/cnv_expression/input/')
    ap.add_argument('--out-cnv', required=True)
    ap.add_argument('--out-expression', required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_cnv), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_expression), exist_ok=True)

    n_seg, cnv_samples = write_cnv_table(args.staging_dir, args.out_cnv)
    print(f'CNV: {n_seg} segmentos de {len(cnv_samples)} amostras -> {args.out_cnv}', file=sys.stderr)

    n_genes, expr_samples = write_expression_table(args.staging_dir, args.out_expression)
    print(f'Expressão: {n_genes} genes x {len(expr_samples)} amostras -> {args.out_expression}', file=sys.stderr)

    if len(cnv_samples) < 3 or len(expr_samples) < 3:
        print('AVISO: menos de 3 amostras em algum dos lados -- CNVRanger::cnvEQTL '
              'precisa de várias amostras por estado de número de cópias pra calcular '
              'desvio/média (ver docstring deste script).', file=sys.stderr)


if __name__ == '__main__':
    main()
