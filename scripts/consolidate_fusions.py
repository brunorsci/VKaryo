#!/usr/bin/env python3
"""
Consolida fusões gênicas de múltiplos callers (GeneFuse, Arriba, Manta)
Aumenta confiança através de evidência cruzada

Uso:
    python consolidate_fusions.py \
        --genefuse sample_fusions.json \
        --arriba sample_fusions.tsv \
        --manta sample_candidateSV.vcf.gz \
        --output sample_consolidated_fusions.tsv
"""

import argparse
import json
import pandas as pd
import pysam
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FusionConsolidator:
    """
    Consolida fusões de múltiplos callers e aumenta confiança
    """

    def __init__(self, min_support=1):
        self.min_support = min_support
        self.fusions = []

    def parse_genefuse(self, json_file):
        """
        Parse output do GeneFuse (formato JSON)
        """
        logger.info(f"Parsing GeneFuse: {json_file}")

        try:
            with open(json_file) as f:
                data = json.load(f)

            fusions = []
            for fusion in data.get('fusions', []):
                fusions.append({
                    'gene1': fusion['gene1'],
                    'gene2': fusion['gene2'],
                    'chr1': fusion['chr1'],
                    'chr2': fusion['chr2'],
                    'breakpoint1': fusion['break1'],
                    'breakpoint2': fusion['break2'],
                    'supporting_reads': fusion['unique_reads'],
                    'caller': 'GeneFuse',
                    'confidence': self._genefuse_confidence(fusion)
                })

            logger.info(f"  Found {len(fusions)} fusions")
            return fusions

        except Exception as e:
            logger.error(f"Error parsing GeneFuse: {e}")
            return []

    def _genefuse_confidence(self, fusion):
        """
        Calcula score de confiança do GeneFuse
        """
        reads = fusion.get('unique_reads', 0)

        if reads >= 10:
            return 'HIGH'
        elif reads >= 5:
            return 'MEDIUM'
        else:
            return 'LOW'

    def parse_arriba(self, tsv_file):
        """
        Parse output do Arriba (formato TSV)
        """
        logger.info(f"Parsing Arriba: {tsv_file}")

        try:
            df = pd.read_csv(tsv_file, sep='\t')

            fusions = []
            for _, row in df.iterrows():
                # Parse gene names (pode ter isoformas)
                gene1 = row['#gene1'].split('(')[0]
                gene2 = row['gene2'].split('(')[0]

                fusions.append({
                    'gene1': gene1,
                    'gene2': gene2,
                    'chr1': row['breakpoint1'].split(':')[0],
                    'chr2': row['breakpoint2'].split(':')[0],
                    'breakpoint1': int(row['breakpoint1'].split(':')[1]),
                    'breakpoint2': int(row['breakpoint2'].split(':')[1]),
                    'supporting_reads': row['split_reads1'] + row['split_reads2'],
                    'confidence': row['confidence'],
                    'caller': 'Arriba',
                    'direction1': row['direction1'],
                    'direction2': row['direction2'],
                    'reading_frame': row['reading_frame']
                })

            logger.info(f"  Found {len(fusions)} fusions")
            return fusions

        except Exception as e:
            logger.error(f"Error parsing Arriba: {e}")
            return []

    def parse_manta(self, vcf_file):
        """
        Parse SVs do Manta que podem ser fusões (BND type)
        """
        logger.info(f"Parsing Manta: {vcf_file}")

        try:
            vcf = pysam.VariantFile(vcf_file)
            fusions = []

            for variant in vcf:
                # Apenas BND (breakends) são potencialmente fusões
                if 'BND' not in variant.info.get('SVTYPE', ''):
                    continue

                # Parse ALT field para pegar parceiro
                alt = str(variant.alts[0]) if variant.alts else ''

                # Formato típico: ]chr:pos] or [chr:pos[
                if ']' in alt or '[' in alt:
                    # Extrair chr:pos do parceiro
                    mate_info = alt.replace(']', '').replace('[', '').split(':')

                    if len(mate_info) >= 2:
                        mate_chr = mate_info[0]
                        mate_pos = int(mate_info[1])

                        # Tentativa de anotar genes (precisaria GTF)
                        # Por ora, apenas reporta coordenadas
                        fusions.append({
                            'gene1': 'Unknown',  # Precisa anotação
                            'gene2': 'Unknown',
                            'chr1': variant.chrom,
                            'chr2': mate_chr,
                            'breakpoint1': variant.pos,
                            'breakpoint2': mate_pos,
                            'supporting_reads': variant.info.get('PR', [0])[0] + variant.info.get('SR', [0])[0],
                            'caller': 'Manta',
                            'confidence': 'MEDIUM'
                        })

            logger.info(f"  Found {len(fusions)} BND events")
            return fusions

        except Exception as e:
            logger.error(f"Error parsing Manta: {e}")
            return []

    def match_fusions(self, fusion1, fusion2, max_distance=10000):
        """
        Verifica se duas fusões são a mesma (com tolerância de distância)
        """
        # Genes devem ser os mesmos (ou reverso)
        genes_match = (
            (fusion1['gene1'] == fusion2['gene1'] and fusion1['gene2'] == fusion2['gene2']) or
            (fusion1['gene1'] == fusion2['gene2'] and fusion1['gene2'] == fusion2['gene1'])
        )

        if not genes_match:
            return False

        # Breakpoints devem estar próximos
        if fusion1['chr1'] == fusion2['chr1'] and fusion1['chr2'] == fusion2['chr2']:
            dist1 = abs(fusion1['breakpoint1'] - fusion2['breakpoint1'])
            dist2 = abs(fusion1['breakpoint2'] - fusion2['breakpoint2'])

            if dist1 <= max_distance and dist2 <= max_distance:
                return True

        return False

    def consolidate(self, genefuse_file=None, arriba_file=None, manta_file=None):
        """
        Consolida fusões de todos os callers
        """
        all_fusions = []

        # Parse cada caller
        if genefuse_file:
            all_fusions.extend(self.parse_genefuse(genefuse_file))

        if arriba_file:
            all_fusions.extend(self.parse_arriba(arriba_file))

        if manta_file:
            all_fusions.extend(self.parse_manta(manta_file))

        logger.info(f"Total de {len(all_fusions)} fusões de todos os callers")

        # Agrupar fusões similares
        consolidated = []
        used = set()

        for i, fusion in enumerate(all_fusions):
            if i in used:
                continue

            # Encontra todas as fusões que correspondem
            matches = [fusion]
            callers = [fusion['caller']]

            for j, other in enumerate(all_fusions[i+1:], start=i+1):
                if j in used:
                    continue

                if self.match_fusions(fusion, other):
                    matches.append(other)
                    callers.append(other['caller'])
                    used.add(j)

            # Cria fusão consolidada
            support_count = len(matches)
            total_reads = sum(f['supporting_reads'] for f in matches)

            # Usa informação do caller mais confiável (Arriba > GeneFuse > Manta)
            best_fusion = matches[0]
            if any(f['caller'] == 'Arriba' for f in matches):
                best_fusion = [f for f in matches if f['caller'] == 'Arriba'][0]

            consolidated_fusion = {
                **best_fusion,
                'support_count': support_count,
                'callers': ','.join(set(callers)),
                'total_supporting_reads': total_reads,
                'consolidated_confidence': self._calculate_confidence(support_count, total_reads)
            }

            # Adiciona informações clínicas se disponível
            consolidated_fusion['clinical_significance'] = self._get_clinical_significance(
                consolidated_fusion['gene1'],
                consolidated_fusion['gene2']
            )

            consolidated.append(consolidated_fusion)

        # Filtra por suporte mínimo
        consolidated = [f for f in consolidated if f['support_count'] >= self.min_support]

        # Ordena por confiança
        consolidated = sorted(
            consolidated,
            key=lambda x: (x['support_count'], x['total_supporting_reads']),
            reverse=True
        )

        logger.info(f"Consolidated to {len(consolidated)} high-confidence fusions")
        return consolidated

    def _calculate_confidence(self, support_count, total_reads):
        """
        Calcula confiança consolidada
        """
        if support_count >= 3 and total_reads >= 20:
            return 'VERY_HIGH'
        elif support_count >= 2 and total_reads >= 10:
            return 'HIGH'
        elif support_count >= 2 or total_reads >= 10:
            return 'MEDIUM'
        else:
            return 'LOW'

    def _get_clinical_significance(self, gene1, gene2):
        """
        Retorna significância clínica de fusões conhecidas
        """
        fusion_name = f"{gene1}-{gene2}"
        reverse_name = f"{gene2}-{gene1}"

        # Database de fusões recorrentes em leucemias
        known_fusions = {
            'BCR-ABL1': {
                'leukemia': 'CML, ALL Ph+',
                'significance': 'Diagnostic and therapeutic target',
                'therapy': 'TKI (Imatinib, Dasatinib, Nilotinib)',
                'prognosis': 'Good with TKI'
            },
            'PML-RARA': {
                'leukemia': 'APL (M3)',
                'significance': 'Diagnostic',
                'therapy': 'ATRA + Arsenic trioxide',
                'prognosis': 'Excellent with treatment'
            },
            'RUNX1-RUNX1T1': {
                'leukemia': 'AML-M2 with t(8;21)',
                'significance': 'Favorable prognostic',
                'therapy': 'Standard chemotherapy',
                'prognosis': 'Favorable'
            },
            'CBFB-MYH11': {
                'leukemia': 'AML with inv(16)',
                'significance': 'Favorable prognostic',
                'therapy': 'Standard chemotherapy + high-dose cytarabine',
                'prognosis': 'Favorable'
            },
            'KMT2A-MLLT3': {
                'leukemia': 'AML with t(9;11)',
                'significance': 'Intermediate prognosis',
                'therapy': 'Standard chemotherapy',
                'prognosis': 'Intermediate'
            },
            'DEK-NUP214': {
                'leukemia': 'AML with t(6;9)',
                'significance': 'Poor prognosis',
                'therapy': 'Intensive chemotherapy',
                'prognosis': 'Poor'
            },
            'ETV6-RUNX1': {
                'leukemia': 'B-ALL pediatric',
                'significance': 'Favorable prognosis',
                'therapy': 'Standard ALL protocol',
                'prognosis': 'Favorable'
            },
            'TCF3-PBX1': {
                'leukemia': 'B-ALL with t(1;19)',
                'significance': 'Intermediate prognosis',
                'therapy': 'Standard ALL protocol',
                'prognosis': 'Improved with modern therapy'
            },
        }

        if fusion_name in known_fusions:
            return known_fusions[fusion_name]
        elif reverse_name in known_fusions:
            return known_fusions[reverse_name]
        else:
            return {
                'leukemia': 'Unknown',
                'significance': 'Novel or rare fusion - requires investigation',
                'therapy': 'Standard protocol',
                'prognosis': 'Unknown'
            }

    def save_results(self, fusions, output_file):
        """
        Salva resultados consolidados
        """
        logger.info(f"Saving results to {output_file}")

        # Flatten clinical significance
        for fusion in fusions:
            clin = fusion.pop('clinical_significance', {})
            fusion['leukemia_type'] = clin.get('leukemia', 'Unknown')
            fusion['clinical_significance'] = clin.get('significance', 'Unknown')
            fusion['therapy_recommendation'] = clin.get('therapy', 'Standard')
            fusion['prognosis'] = clin.get('prognosis', 'Unknown')

        df = pd.DataFrame(fusions)

        # Reordena colunas
        column_order = [
            'gene1', 'gene2', 'chr1', 'chr2', 'breakpoint1', 'breakpoint2',
            'support_count', 'callers', 'total_supporting_reads',
            'consolidated_confidence', 'leukemia_type', 'clinical_significance',
            'therapy_recommendation', 'prognosis'
        ]

        # Adiciona colunas que existem
        column_order = [c for c in column_order if c in df.columns]
        other_columns = [c for c in df.columns if c not in column_order]
        df = df[column_order + other_columns]

        # Salva TSV
        df.to_csv(output_file, sep='\t', index=False)

        logger.info(f"✓ Saved {len(df)} consolidated fusions")

        # Estatísticas
        logger.info("\nSummary:")
        logger.info(f"  Very High confidence: {sum(df['consolidated_confidence'] == 'VERY_HIGH')}")
        logger.info(f"  High confidence: {sum(df['consolidated_confidence'] == 'HIGH')}")
        logger.info(f"  Medium confidence: {sum(df['consolidated_confidence'] == 'MEDIUM')}")
        logger.info(f"  Low confidence: {sum(df['consolidated_confidence'] == 'LOW')}")
        logger.info(f"\n  Known pathogenic fusions: {sum(df['clinical_significance'] != 'Novel or rare fusion - requires investigation')}")


def main():
    parser = argparse.ArgumentParser(
        description='Consolida fusões gênicas de múltiplos callers'
    )
    parser.add_argument('--genefuse', help='GeneFuse output JSON')
    parser.add_argument('--arriba', help='Arriba output TSV')
    parser.add_argument('--manta', help='Manta VCF')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--min-support', type=int, default=1,
                        help='Minimum number of callers supporting fusion')

    args = parser.parse_args()

    consolidator = FusionConsolidator(min_support=args.min_support)

    fusions = consolidator.consolidate(
        genefuse_file=args.genefuse,
        arriba_file=args.arriba,
        manta_file=args.manta
    )

    consolidator.save_results(fusions, args.output)


if __name__ == '__main__':
    main()
