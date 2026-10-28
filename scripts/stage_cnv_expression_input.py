#!/usr/bin/env python3
"""
Copia os arquivos por amostra necessários pra análise de Expressão de CNVs
(CNVkit .cns + STAR ReadsPerGene.out.tab) pra uma pasta de input dedicada
dentro de output_vkaryo/cnv_expression/input/ -- mantém um snapshot auditável
de exatamente quais arquivos entraram na análise, separado das pastas gerais
do pipeline (cnvkit/, star_aligned/).

Só copia amostras com DNA+RNA pareados (a análise precisa de ambos).
"""
import argparse
import csv
import os
import shutil
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', required=True, help='output_vkaryo/ do projeto')
    ap.add_argument('--sample-mapping', required=True)
    ap.add_argument('--samples', default=None, help='IDs de DNA separados por vírgula (default: todos os pares)')
    ap.add_argument('--staging-dir', required=True, help='output_vkaryo/cnv_expression/input/')
    args = ap.parse_args()

    only = set(args.samples.split(',')) if args.samples else None
    os.makedirs(args.staging_dir, exist_ok=True)

    n_copied = 0
    with open(args.sample_mapping, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            dna, rna = row['DNA_sample'], row['RNA_sample']
            if only and dna not in only:
                continue

            cns_src = os.path.join(args.output_dir, 'cnvkit', f'{dna}_markduplicates_recal.cns')
            expr_src = os.path.join(args.output_dir, 'star_aligned', f'{rna}_ReadsPerGene.out.tab')
            if not (os.path.exists(cns_src) and os.path.exists(expr_src)):
                print(f'  [aviso] {dna}/{rna}: falta .cns ou ReadsPerGene.out.tab, pulando', file=sys.stderr)
                continue

            shutil.copy2(cns_src, os.path.join(args.staging_dir, f'{dna}_cnv.cns'))
            shutil.copy2(expr_src, os.path.join(args.staging_dir, f'{dna}_expression.tab'))
            n_copied += 1

    print(f'{n_copied} amostras copiadas para {args.staging_dir}', file=sys.stderr)
    if n_copied < 3:
        print('AVISO: menos de 3 amostras -- CNVRanger::cnvEQTL precisa de várias '
              'amostras por estado de número de cópias pra calcular desvio/média.', file=sys.stderr)


if __name__ == '__main__':
    main()
