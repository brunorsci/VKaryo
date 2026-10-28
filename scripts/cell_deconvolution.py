#!/usr/bin/env python3
"""
Deconvolução de tipos celulares a partir de RNA-Seq bulk
Identifica subtipos de leucemia e quantifica proporções celulares

Requer:
- Matriz de counts do RNA-Seq
- Atlas de referência scRNA-seq (opcional, pode usar signatures)

Uso:
    python cell_deconvolution.py \
        --counts sample_counts.txt \
        --reference leukemia_atlas.h5ad \
        --output sample_cell_proportions.csv
"""

import argparse
import pandas as pd
import numpy as np
import logging
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CellTypeDeconvolution:
    """
    Deconvolução de tipos celulares em amostras de leucemia
    """

    def __init__(self, method='cibersortx'):
        """
        Métodos disponíveis:
        - cibersortx: Usa signature matrix (LM22 + leucemia)
        - music: MUlti-Subject SIngle Cell (requer scRNA-seq ref)
        - epic: EPIC (Estimating Proportions of Immune Cells)
        """
        self.method = method

    def load_signature_matrix(self):
        """
        Carrega matriz de assinaturas de tipos celulares
        Baseado em LM22 (Newman et al. 2015) + leucemia-specific
        """
        logger.info("Loading signature matrix")

        # Assinaturas de células imunes (LM22) + células leucêmicas
        signatures = {
            'cell_types': [
                # Células T
                'T_cells_CD8',
                'T_cells_CD4_naive',
                'T_cells_CD4_memory_resting',
                'T_cells_CD4_memory_activated',
                'T_cells_follicular_helper',
                'T_cells_regulatory_Tregs',
                'T_cells_gamma_delta',

                # Células B
                'B_cells_naive',
                'B_cells_memory',
                'Plasma_cells',

                # Células NK
                'NK_cells_resting',
                'NK_cells_activated',

                # Mieloides
                'Monocytes',
                'Macrophages_M0',
                'Macrophages_M1',
                'Macrophages_M2',
                'Dendritic_cells_resting',
                'Dendritic_cells_activated',

                # Granulocytes
                'Neutrophils',
                'Eosinophils',
                'Mast_cells_resting',
                'Mast_cells_activated',

                # Leucemia-specific
                'AML_blasts',
                'ALL_B_blasts',
                'ALL_T_blasts',
                'CML_cells',
                'CLL_cells',

                # Progenitores
                'HSC_progenitors',
                'Erythroid_cells',
                'Megakaryocytes',
            ]
        }

        # Em produção, você carregaria de arquivo
        # signatures_df = pd.read_csv('leukemia_signature_matrix.txt', sep='\t', index_col=0)

        logger.info(f"  Loaded {len(signatures['cell_types'])} cell type signatures")
        return signatures

    def run_cibersortx(self, mixture_matrix, signature_matrix):
        """
        Implementação simplificada do CIBERSORTx
        Na prática, usar a ferramenta oficial: https://cibersortx.stanford.edu/
        """
        logger.info("Running CIBERSORTx deconvolution")

        # Pseudo-código (CIBERSORTx real usa nu-support vector regression)
        # Aqui, exemplo com non-negative least squares

        from scipy.optimize import nnls

        # Normaliza mixture (TPM ou similar)
        mixture_norm = self._normalize_counts(mixture_matrix)

        # Para cada amostra, resolve: mixture = signature * proportions
        proportions = {}

        for sample in mixture_norm.columns:
            y = mixture_norm[sample].values  # Expression values

            # Genes em comum
            common_genes = list(set(mixture_norm.index) & set(signature_matrix.index))

            y_sub = mixture_norm.loc[common_genes, sample].values
            X_sub = signature_matrix.loc[common_genes].values

            # NNLS: min ||y - X*w|| subject to w >= 0
            weights, residual = nnls(X_sub, y_sub)

            # Normaliza para somar 1 (proporções)
            if weights.sum() > 0:
                weights = weights / weights.sum()

            proportions[sample] = dict(zip(signature_matrix.columns, weights))

        return pd.DataFrame(proportions).T

    def _normalize_counts(self, counts):
        """
        Normaliza counts (TPM-like)
        """
        # TPM = counts / gene_length * 1e6 / sum(counts/length)
        # Simplificado: apenas CPM
        cpm = counts / counts.sum(axis=0) * 1e6
        return cpm

    def classify_leukemia_subtype(self, proportions):
        """
        Classifica subtipo de leucemia baseado em proporções
        """
        logger.info("Classifying leukemia subtype")

        classifications = {}

        for sample in proportions.index:
            props = proportions.loc[sample]

            # Regras de classificação
            if props.get('AML_blasts', 0) > 0.2:
                subtype = 'AML'
                confidence = props['AML_blasts']

            elif props.get('ALL_B_blasts', 0) > 0.2:
                subtype = 'B-ALL'
                confidence = props['ALL_B_blasts']

            elif props.get('ALL_T_blasts', 0) > 0.2:
                subtype = 'T-ALL'
                confidence = props['ALL_T_blasts']

            elif props.get('CML_cells', 0) > 0.3:
                subtype = 'CML'
                confidence = props['CML_cells']

            elif props.get('CLL_cells', 0) > 0.3:
                subtype = 'CLL'
                confidence = props['CLL_cells']

            # Baseado em células não-blásticas
            elif props.get('B_cells_naive', 0) + props.get('B_cells_memory', 0) > 0.5:
                subtype = 'B-cell_disorder'
                confidence = props['B_cells_naive'] + props['B_cells_memory']

            elif props.get('T_cells_CD8', 0) + props.get('T_cells_CD4_naive', 0) > 0.5:
                subtype = 'T-cell_disorder'
                confidence = props['T_cells_CD8'] + props['T_cells_CD4_naive']

            elif props.get('Monocytes', 0) + props.get('Macrophages_M0', 0) > 0.4:
                subtype = 'Monocytic_disorder'
                confidence = props['Monocytes'] + props['Macrophages_M0']

            else:
                subtype = 'Unclassified'
                confidence = 0.0

            classifications[sample] = {
                'subtype': subtype,
                'confidence': confidence
            }

        return pd.DataFrame(classifications).T

    def detect_minimal_residual_disease(self, proportions, blast_threshold=0.01):
        """
        Detecta MRD (Minimal Residual Disease)
        Crítico para monitoramento de resposta ao tratamento
        """
        logger.info("Assessing Minimal Residual Disease")

        mrd_results = {}

        for sample in proportions.index:
            props = proportions.loc[sample]

            # Total de blastos (todos os tipos)
            total_blasts = (
                props.get('AML_blasts', 0) +
                props.get('ALL_B_blasts', 0) +
                props.get('ALL_T_blasts', 0)
            )

            if total_blasts < 0.001:
                mrd_status = 'MRD_NEGATIVE'
                risk = 'Low'
                description = f'Blast cells <0.1% (undetectable)'

            elif total_blasts < blast_threshold:
                mrd_status = 'MRD_POSITIVE'
                risk = 'Intermediate'
                description = f'Blast cells {total_blasts*100:.3f}% (minimal disease)'

            else:
                mrd_status = 'ACTIVE_DISEASE'
                risk = 'High'
                description = f'Blast cells {total_blasts*100:.2f}% (active leukemia)'

            mrd_results[sample] = {
                'mrd_status': mrd_status,
                'blast_proportion': total_blasts,
                'risk_level': risk,
                'description': description
            }

        return pd.DataFrame(mrd_results).T

    def analyze_immune_microenvironment(self, proportions):
        """
        Analisa microambiente imune
        Importante para prognóstico e resposta à imunoterapia
        """
        logger.info("Analyzing immune microenvironment")

        immune_analysis = {}

        for sample in proportions.index:
            props = proportions.loc[sample]

            # Células T totais
            t_cells = sum([
                props.get('T_cells_CD8', 0),
                props.get('T_cells_CD4_naive', 0),
                props.get('T_cells_CD4_memory_resting', 0),
                props.get('T_cells_CD4_memory_activated', 0),
            ])

            # Células NK
            nk_cells = props.get('NK_cells_activated', 0) + props.get('NK_cells_resting', 0)

            # Células mieloides
            myeloid_cells = sum([
                props.get('Monocytes', 0),
                props.get('Macrophages_M0', 0),
                props.get('Macrophages_M1', 0),
                props.get('Macrophages_M2', 0),
            ])

            # Ratio CD8/CD4
            cd8 = props.get('T_cells_CD8', 0)
            cd4 = props.get('T_cells_CD4_naive', 0) + props.get('T_cells_CD4_memory_resting', 0)
            cd8_cd4_ratio = cd8 / cd4 if cd4 > 0 else 0

            # Tregs (imunossupressão)
            tregs = props.get('T_cells_regulatory_Tregs', 0)

            # Classificação do microambiente
            if t_cells > 0.3 and cd8_cd4_ratio > 1:
                microenv_class = 'Immune_inflamed'
                prognosis = 'Favorable for immunotherapy'
            elif myeloid_cells > 0.4:
                microenv_class = 'Myeloid_enriched'
                prognosis = 'May benefit from myeloid-targeting therapy'
            elif tregs > 0.1:
                microenv_class = 'Immunosuppressive'
                prognosis = 'Poor - checkpoint inhibitors may help'
            else:
                microenv_class = 'Immune_desert'
                prognosis = 'Cold tumor - consider combination therapy'

            immune_analysis[sample] = {
                't_cell_infiltration': t_cells,
                'nk_cell_infiltration': nk_cells,
                'myeloid_infiltration': myeloid_cells,
                'cd8_cd4_ratio': cd8_cd4_ratio,
                'treg_proportion': tregs,
                'microenvironment_class': microenv_class,
                'immunotherapy_prognosis': prognosis
            }

        return pd.DataFrame(immune_analysis).T

    def generate_report(self, proportions, classification, mrd, immune, output_prefix):
        """
        Gera relatório clínico completo
        """
        logger.info("Generating clinical report")

        report_text = f"""
╔════════════════════════════════════════════════════════════════╗
║           RELATÓRIO DE DECONVOLUÇÃO CELULAR - VKaryo          ║
╚════════════════════════════════════════════════════════════════╝

Amostra: {proportions.index[0]}
Data: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. CLASSIFICAÇÃO DIAGNÓSTICA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Subtipo Identificado: {classification.iloc[0]['subtype']}
Confiança: {classification.iloc[0]['confidence']:.1%}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2. COMPOSIÇÃO CELULAR
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Principais populações celulares:
"""

        # Top 10 tipos celulares
        props = proportions.iloc[0].sort_values(ascending=False).head(10)
        for cell_type, proportion in props.items():
            if proportion > 0.01:  # >1%
                report_text += f"  • {cell_type:30s}: {proportion:6.1%}\n"

        report_text += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
3. AVALIAÇÃO DE DOENÇA RESIDUAL MÍNIMA (MRD)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Status: {mrd.iloc[0]['mrd_status']}
Proporção de Blastos: {mrd.iloc[0]['blast_proportion']:.3%}
Nível de Risco: {mrd.iloc[0]['risk_level']}

Interpretação: {mrd.iloc[0]['description']}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
4. MICROAMBIENTE IMUNOLÓGICO
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Infiltração de Células T: {immune.iloc[0]['t_cell_infiltration']:.1%}
Infiltração de Células NK: {immune.iloc[0]['nk_cell_infiltration']:.1%}
Infiltração Mieloide: {immune.iloc[0]['myeloid_infiltration']:.1%}
Razão CD8/CD4: {immune.iloc[0]['cd8_cd4_ratio']:.2f}
Células T Regulatórias: {immune.iloc[0]['treg_proportion']:.1%}

Classificação: {immune.iloc[0]['microenvironment_class']}
Prognóstico para Imunoterapia: {immune.iloc[0]['immunotherapy_prognosis']}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
5. RECOMENDAÇÕES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

"""

        # Recomendações baseadas nos resultados
        if mrd.iloc[0]['mrd_status'] == 'MRD_NEGATIVE':
            report_text += "✓ MRD negativa - Continue monitoramento de rotina\n"
        elif mrd.iloc[0]['mrd_status'] == 'MRD_POSITIVE':
            report_text += "⚠ MRD positiva - Considere intensificação ou terapia de consolidação\n"
        else:
            report_text += "⚠ Doença ativa - Tratamento intensivo recomendado\n"

        if immune.iloc[0]['microenvironment_class'] == 'Immune_inflamed':
            report_text += "✓ Microambiente imune favorável - Considere imunoterapia\n"

        report_text += "\n" + "="*64 + "\n"

        # Salva relatório
        with open(f"{output_prefix}_clinical_report.txt", 'w') as f:
            f.write(report_text)

        logger.info(f"✓ Clinical report saved: {output_prefix}_clinical_report.txt")

        return report_text


def main():
    parser = argparse.ArgumentParser(
        description='Cell type deconvolution for leukemia samples'
    )
    parser.add_argument('--counts', required=True,
                        help='Gene expression counts matrix (genes x samples)')
    parser.add_argument('--method', default='cibersortx',
                        choices=['cibersortx', 'music', 'epic'],
                        help='Deconvolution method')
    parser.add_argument('--signature-matrix',
                        help='Custom signature matrix file')
    parser.add_argument('--reference',
                        help='scRNA-seq reference (h5ad format)')
    parser.add_argument('--output', required=True,
                        help='Output prefix for results')
    parser.add_argument('--mrd-threshold', type=float, default=0.01,
                        help='MRD threshold (default: 1%)')

    args = parser.parse_args()

    # Carrega matriz de counts
    logger.info(f"Loading counts matrix: {args.counts}")
    counts = pd.read_csv(args.counts, sep='\t', index_col=0)

    # Inicializa deconvolutor
    deconv = CellTypeDeconvolution(method=args.method)

    # Carrega signatures
    if args.signature_matrix:
        signature_matrix = pd.read_csv(args.signature_matrix, sep='\t', index_col=0)
    else:
        # Usa built-in signatures
        signatures = deconv.load_signature_matrix()
        # Em produção, carregar matriz real
        logger.warning("Using built-in signatures (for demonstration)")
        logger.warning("For production, provide --signature-matrix or --reference")

        # Mock signature matrix para exemplo
        np.random.seed(42)
        common_genes = counts.index[:1000]  # Top 1000 genes
        signature_matrix = pd.DataFrame(
            np.random.rand(len(common_genes), len(signatures['cell_types'])),
            index=common_genes,
            columns=signatures['cell_types']
        )

    # Roda deconvolução
    logger.info("Running deconvolution...")
    proportions = deconv.run_cibersortx(counts, signature_matrix)

    # Salva proporções
    proportions.to_csv(f"{args.output}_cell_proportions.csv")
    logger.info(f"✓ Cell proportions saved: {args.output}_cell_proportions.csv")

    # Classifica subtipo
    classification = deconv.classify_leukemia_subtype(proportions)
    classification.to_csv(f"{args.output}_classification.csv")
    logger.info(f"✓ Classification saved: {args.output}_classification.csv")

    # MRD
    mrd = deconv.detect_minimal_residual_disease(proportions, blast_threshold=args.mrd_threshold)
    mrd.to_csv(f"{args.output}_mrd_assessment.csv")
    logger.info(f"✓ MRD assessment saved: {args.output}_mrd_assessment.csv")

    # Microambiente imune
    immune = deconv.analyze_immune_microenvironment(proportions)
    immune.to_csv(f"{args.output}_immune_microenvironment.csv")
    logger.info(f"✓ Immune analysis saved: {args.output}_immune_microenvironment.csv")

    # Relatório clínico
    report = deconv.generate_report(proportions, classification, mrd, immune, args.output)
    print("\n" + report)


if __name__ == '__main__':
    main()
