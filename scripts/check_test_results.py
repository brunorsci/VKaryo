#!/usr/bin/env python3
"""Compare a run of the worked example against the calls it should produce.

Exits non-zero on any mismatch, so the example can be wired into CI. Reads
test/expected/expected_calls.tsv, written by make_test_dataset.py, and looks
for each call in the corresponding output table.
"""

import argparse
import csv
import glob
import io
import os
import sys

TOLERANCIA_PB = 50_000   # margem nas bordas do evento de numero de copias


def ler_esperado(caminho):
    with io.open(caminho, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def achou_fusao(outdir, amostra, detalhe):
    a, b = detalhe.split('--')
    for p in glob.glob(os.path.join(outdir, '**', f'{amostra}_fusions.tsv'), recursive=True):
        with io.open(p, encoding='utf-8', errors='replace') as fh:
            for linha in fh:
                if a in linha and b in linha:
                    return True
    return False


def achou_perda(outdir, amostra, chrom, inicio, fim):
    for p in glob.glob(os.path.join(outdir, '**', f'{amostra}_focal_events.tsv'), recursive=True):
        with io.open(p, encoding='utf-8', errors='replace') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r.get('chrom') != chrom:
                    continue
                if 'homozig' not in (r.get('classe_do_evento') or '').lower():
                    continue
                try:
                    s, e = int(r['start']), int(r['end'])
                except (KeyError, ValueError):
                    continue
                if s < fim + TOLERANCIA_PB and e > inicio - TOLERANCIA_PB:
                    return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', required=True, help='diretorio de saida da execucao')
    ap.add_argument('--expected', default='test/expected', help='diretorio com expected_calls.tsv')
    args = ap.parse_args()

    esperados = ler_esperado(os.path.join(args.expected, 'expected_calls.tsv'))
    falhas = []
    for e in esperados:
        amostra, tipo = e['sample'], e['type']
        if tipo == 'fusion':
            ok = achou_fusao(args.outdir, amostra, e['detail'])
        elif tipo == 'homozygous_loss':
            ok = achou_perda(args.outdir, amostra, e['chrom'], int(e['start']), int(e['end']))
        else:
            ok = False
        marca = 'OK  ' if ok else 'FAIL'
        print(f'{marca} {amostra:<9} {tipo:<18} {e["detail"]}')
        if not ok:
            falhas.append((amostra, tipo, e['detail']))

    print()
    if falhas:
        print(f'{len(falhas)} of {len(esperados)} expected calls were not found.')
        sys.exit(1)
    print(f'all {len(esperados)} expected calls found.')


if __name__ == '__main__':
    main()
