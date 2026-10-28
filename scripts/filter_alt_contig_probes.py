#!/usr/bin/env python3
"""
Filtra sondas do BED de captura (Agilent SureSelect V6) que referenciam
contigs alternativos ("_alt", "_random", "chrUn_...") ausentes do FASTA de
referência usado pelo pipeline (GRCh38_character_genome.fasta, 195 contigs,
vs. os 455 do hg38.fa completo usado para desenhar o kit de captura).

NOTA (Seção 6.6 da tese): essa incompatibilidade é a causa raiz do KeyError
de alt-contig que forçaria o CNVkit a rodar com `--no-gc --no-rmask`
(referência "flat", sem correção de conteúdo GC/RepeatMasker) na etapa
cnvkit_reference -- o que por sua vez produz um desvio sistemático de log2
GLOBAL (quase todo o genoma "ganho" ou "perdido" de uma vez, ver
karyotype_quality_flag em classify_leukemia.py). Apenas 151/244655 sondas
(0,06%) referenciam contigs ausentes; nenhum gene é perdido exclusivamente
por essas sondas (verificado manualmente: todas as sondas alt-contig são
cópias adicionais de regiões já cobertas por sondas no contig primário
correspondente).

Uso:
    python3 filter_alt_contig_probes.py \
        --bed reference/Agilent_SureSelect_V6_hg38.bed \
        --fai reference/GRCh38_character_genome.fasta.fai \
        --output reference/Agilent_SureSelect_V6_hg38.noalt.bed
"""

import argparse


def load_valid_contigs(fai_path):
    contigs = set()
    with open(fai_path) as f:
        for line in f:
            contigs.add(line.split('\t')[0])
    return contigs


def filter_bed(bed_path, valid_contigs, output_path):
    kept, dropped = 0, 0
    with open(bed_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#'):
                fout.write(line)
                continue
            chrom = line.split('\t', 1)[0]
            if chrom in valid_contigs:
                fout.write(line)
                kept += 1
            else:
                dropped += 1
    return kept, dropped


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bed', required=True)
    ap.add_argument('--fai', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    valid = load_valid_contigs(args.fai)
    kept, dropped = filter_bed(args.bed, valid, args.output)
    print(f"Contigs válidos no FASTA de referência: {len(valid)}")
    print(f"Sondas mantidas: {kept}")
    print(f"Sondas descartadas (contig ausente do FASTA): {dropped}")
