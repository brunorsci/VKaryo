#!/usr/bin/env python3
"""
Valida CNVs usando expressão gênica do RNA-Seq
Ganhos devem ter expressão aumentada, perdas expressão diminuída

Uso:
    python validate_cnv_expression.py \
        --cnv output_vkaryo/cnvkit/sample.cns \
        --expression output_vkaryo/star_aligned/sample_ReadsPerGene.out.tab \
        --gtf reference/Homo_sapiens.GRCh38.113.gtf.gz \
        --output output_vkaryo/cnv_validated/sample_cnv_validated.tsv
"""

import pandas as pd
import numpy as np
import argparse
import logging
import gzip
from scipy import stats

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


class CNVExpressionValidator:
    """
    Valida CNVs comparando com expressão gênica
    """

    def __init__(self, gtf_file):
        self.gtf_file = gtf_file
        self.gene_coords = None

    def load_gene_coordinates(self):
        """
        Carrega coordenadas genômicas dos genes do GTF
        """
        logger.info("Carregando coordenadas dos genes do GTF...")

        genes = []

        open_func = gzip.open if self.gtf_file.endswith('.gz') else open

        with open_func(self.gtf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                if fields[2] != 'gene':  # Apenas features do tipo 'gene'
                    continue

                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]

                # Parse atributos
                attr_dict = {}
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    key_val = attr.split(' ', 1)
                    if len(key_val) == 2:
                        key, val = key_val
                        attr_dict[key] = val.strip('"')

                gene_id = attr_dict.get('gene_id', '')
                gene_name = attr_dict.get('gene_name', gene_id)

                genes.append({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                })

        self.gene_coords = pd.DataFrame(genes)
        logger.info(f"  Carregado {len(self.gene_coords)} genes")

        return self.gene_coords

    def load_cnvs(self, cnv_file):
        """
        Carrega CNVs do CNVkit (.cns file)
        """
        logger.info(f"Carregando CNVs: {cnv_file}")

        cnvs = pd.read_csv(cnv_file, sep='\t', comment='#')

        # Renomeia colunas se necessário
        if 'chromosome' in cnvs.columns:
            cnvs.rename(columns={'chromosome': 'chrom'}, inplace=True)

        # Classifica CNVs
        cnvs['cnv_type'] = 'neutral'
        cnvs.loc[cnvs['log2'] > 0.3, 'cnv_type'] = 'gain'
        cnvs.loc[cnvs['log2'] < -0.3, 'cnv_type'] = 'loss'

        # Filtra apenas alterações significativas
        cnvs_altered = cnvs[cnvs['cnv_type'] != 'neutral'].copy()

        logger.info(f"  CNVs totais: {len(cnvs)}")
        logger.info(f"  Ganhos: {(cnvs['cnv_type'] == 'gain').sum()}")
        logger.info(f"  Perdas: {(cnvs['cnv_type'] == 'loss').sum()}")

        return cnvs_altered

    def load_expression(self, expr_file):
        """
        Carrega expressão do STAR ReadsPerGene.out.tab
        """
        logger.info(f"Carregando expressão: {expr_file}")

        # STAR output: gene_id, unstranded, strand1, strand2
        expr = pd.read_csv(expr_file, sep='\t', header=None, skiprows=4,
                          names=['gene_id', 'unstranded', 'strand1', 'strand2'])

        # Usa coluna unstranded (ou strand1/2 se library for stranded)
        expr['counts'] = expr['unstranded']

        # Remove genes sem contagem
        expr = expr[expr['counts'] > 0].copy()

        # Normaliza (TPM-like simplificado: CPM)
        total_counts = expr['counts'].sum()
        expr['cpm'] = (expr['counts'] / total_counts) * 1e6

        # Log2 CPM (com pseudocount)
        expr['log2_cpm'] = np.log2(expr['cpm'] + 1)

        # Z-score (para detectar outliers)
        expr['zscore'] = stats.zscore(expr['log2_cpm'])

        logger.info(f"  Genes expressos: {len(expr)}")
        logger.info(f"  Mediana expressão (CPM): {expr['cpm'].median():.2f}")

        return expr

    def find_genes_in_cnv_region(self, cnv_region):
        """
        Encontra genes dentro de uma região CNV
        """
        chrom = cnv_region['chrom']
        start = cnv_region['start']
        end = cnv_region['end']

        # Busca genes que overlap com região CNV
        genes_in_region = self.gene_coords[
            (self.gene_coords['chrom'] == chrom) &
            (self.gene_coords['start'] <= end) &
            (self.gene_coords['end'] >= start)
        ].copy()

        return genes_in_region

    def validate_cnv_with_expression(self, cnvs, expression):
        """
        Valida cada CNV comparando com expressão dos genes na região
        """
        logger.info("Validando CNVs com expressão...")

        validated_cnvs = []

        for idx, cnv in cnvs.iterrows():
            # Genes na região
            genes_in_region = self.find_genes_in_cnv_region(cnv)

            if len(genes_in_region) == 0:
                # Sem genes na região (intergênica)
                validated_cnvs.append({
                    **cnv.to_dict(),
                    'n_genes': 0,
                    'expressed_genes': 0,
                    'median_expression_zscore': np.nan,
                    'expression_supports_cnv': 'NO_GENES',
                    'validation_status': 'INTERGENIC',
                    'confidence': 'LOW'
                })
                continue

            # Match expressão
            gene_ids = genes_in_region['gene_id'].tolist()
            expr_in_region = expression[expression['gene_id'].isin(gene_ids)]

            if len(expr_in_region) == 0:
                # Genes não expressos (pode ser tecido-específico)
                validated_cnvs.append({
                    **cnv.to_dict(),
                    'n_genes': len(genes_in_region),
                    'expressed_genes': 0,
                    'median_expression_zscore': np.nan,
                    'expression_supports_cnv': 'NOT_EXPRESSED',
                    'validation_status': 'NO_EXPRESSION',
                    'confidence': 'LOW'
                })
                continue

            # Calcula mediana de expressão (z-score)
            median_zscore = expr_in_region['zscore'].median()

            # Validação lógica:
            # GAIN deve ter expressão aumentada (zscore > 0)
            # LOSS deve ter expressão diminuída (zscore < 0)

            cnv_type = cnv['cnv_type']
            expression_supports = False
            validation_status = 'INVALID'
            confidence = 'LOW'

            if cnv_type == 'gain' and median_zscore > 0.5:
                expression_supports = True
                validation_status = 'VALIDATED'
                confidence = 'HIGH' if median_zscore > 1.0 else 'MEDIUM'

            elif cnv_type == 'loss' and median_zscore < -0.5:
                expression_supports = True
                validation_status = 'VALIDATED'
                confidence = 'HIGH' if median_zscore < -1.0 else 'MEDIUM'

            elif abs(median_zscore) < 0.5:
                # Expressão neutra - CNV pode ser falso positivo
                validation_status = 'CONFLICTING'
                confidence = 'LOW'

            else:
                # Expressão vai na direção oposta - likely FP
                validation_status = 'INVALID'
                confidence = 'VERY_LOW'

            # Lista genes
            gene_list = ';'.join(genes_in_region['gene_name'].tolist()[:10])  # Max 10
            if len(genes_in_region) > 10:
                gene_list += f";...+{len(genes_in_region)-10}_more"

            validated_cnvs.append({
                **cnv.to_dict(),
                'n_genes': len(genes_in_region),
                'expressed_genes': len(expr_in_region),
                'gene_list': gene_list,
                'median_expression_zscore': median_zscore,
                'expression_supports_cnv': 'YES' if expression_supports else 'NO',
                'validation_status': validation_status,
                'confidence': confidence
            })

        df_validated = pd.DataFrame(validated_cnvs)

        logger.info(f"\n=== RESULTADOS DA VALIDAÇÃO ===")
        logger.info(f"Total CNVs: {len(df_validated)}")
        logger.info(f"Validados (HIGH conf): {(df_validated['confidence'] == 'HIGH').sum()}")
        logger.info(f"Validados (MEDIUM conf): {(df_validated['confidence'] == 'MEDIUM').sum()}")
        logger.info(f"Conflitantes: {(df_validated['validation_status'] == 'CONFLICTING').sum()}")
        logger.info(f"Inválidos: {(df_validated['validation_status'] == 'INVALID').sum()}")

        return df_validated

    def summarize_by_chromosome(self, validated_cnvs):
        """
        Sumariza alterações por cromossomo
        """
        logger.info("\n=== ALTERAÇÕES POR CROMOSSOMO ===")

        # Filtra apenas validados
        validated = validated_cnvs[validated_cnvs['validation_status'] == 'VALIDATED']

        chr_summary = validated.groupby(['chrom', 'cnv_type']).agg({
            'log2': 'mean',
            'n_genes': 'sum',
            'confidence': lambda x: (x.isin(['HIGH', 'MEDIUM'])).sum()
        }).reset_index()

        for chrom in sorted(validated['chrom'].unique()):
            chr_data = chr_summary[chr_summary['chrom'] == chrom]

            gains = chr_data[chr_data['cnv_type'] == 'gain']
            losses = chr_data[chr_data['cnv_type'] == 'loss']

            line = f"{chrom}: "
            if len(gains) > 0:
                line += f"+{gains['n_genes'].sum():.0f} genes gained"
            if len(losses) > 0:
                if len(gains) > 0:
                    line += " | "
                line += f"-{losses['n_genes'].sum():.0f} genes lost"

            logger.info(line)


def main():
    parser = argparse.ArgumentParser(
        description='Valida CNVs usando expressão gênica do RNA-Seq'
    )
    parser.add_argument('--cnv', required=True, help='Arquivo CNVkit .cns')
    parser.add_argument('--expression', required=True, help='Arquivo STAR ReadsPerGene.out.tab')
    parser.add_argument('--gtf', required=True, help='Arquivo GTF com anotação de genes')
    parser.add_argument('--output', required=True, help='Arquivo de saída TSV')
    parser.add_argument('--output-high-conf', help='Arquivo apenas CNVs alta confiança (opcional)')

    args = parser.parse_args()

    # Inicializa validador
    validator = CNVExpressionValidator(args.gtf)

    # Carrega dados
    validator.load_gene_coordinates()
    cnvs = validator.load_cnvs(args.cnv)
    expression = validator.load_expression(args.expression)

    # Valida
    validated = validator.validate_cnv_with_expression(cnvs, expression)

    # Sumariza
    validator.summarize_by_chromosome(validated)

    # Salva resultados completos
    validated.to_csv(args.output, sep='\t', index=False)
    logger.info(f"\n✓ Resultados salvos: {args.output}")

    # Salva apenas alta confiança
    if args.output_high_conf:
        high_conf = validated[validated['confidence'].isin(['HIGH', 'MEDIUM'])]
        high_conf.to_csv(args.output_high_conf, sep='\t', index=False)
        logger.info(f"✓ Alta confiança salva: {args.output_high_conf} ({len(high_conf)} CNVs)")

    logger.info("\n" + "="*60)
    logger.info("VALIDAÇÃO COMPLETA!")
    logger.info("="*60)


if __name__ == '__main__':
    main()
