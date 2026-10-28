#!/usr/bin/env python3
"""
Gera relatório PDF completo de cariotipagem virtual
Integra TODOS os resultados: fusões, CNVs, SNVs, aneuploidias, classificação

Uso:
    python generate_pdf_report.py \
        --sample SAMPLE_ID \
        --classification output_vkaryo/classification/sample_classification.json \
        --fusions output_vkaryo/arriba_annotated/sample_fusions_annotated.tsv \
        --cnvs output_vkaryo/cnv_validated/sample_cnv_validated.tsv \
        --ideogram output_vkaryo/plots/sample_aneuploidy_ideogram.pdf \
        --output output_vkaryo/reports/sample_report.pdf
"""

import argparse
import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from datetime import datetime
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


class PDFReportGenerator:
    """
    Gera relatório PDF completo
    """

    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.report_data = {}

    def load_classification(self, classification_file):
        """Carrega classificação JSON"""
        logger.info(f"Carregando classificação: {classification_file}")

        with open(classification_file) as f:
            self.report_data['classification'] = json.load(f)

    def load_fusions(self, fusion_file):
        """Carrega fusões anotadas"""
        logger.info(f"Carregando fusões: {fusion_file}")

        fusions = pd.read_csv(fusion_file, sep='\t')
        self.report_data['fusions'] = fusions

    def load_cnvs(self, cnv_file):
        """Carrega CNVs validados"""
        logger.info(f"Carregando CNVs: {cnv_file}")

        cnvs = pd.read_csv(cnv_file, sep='\t')
        self.report_data['cnvs'] = cnvs

    def create_cover_page(self, pdf):
        """
        Página de capa
        """
        fig = plt.figure(figsize=(8.5, 11))
        ax = fig.add_subplot(111)
        ax.axis('off')

        # Título principal
        ax.text(0.5, 0.85, 'RELATÓRIO DE',
                ha='center', va='center', fontsize=24, fontweight='bold',
                transform=ax.transAxes)

        ax.text(0.5, 0.80, 'CARIOTIPAGEM VIRTUAL',
                ha='center', va='center', fontsize=28, fontweight='bold',
                transform=ax.transAxes, color='#2c3e50')

        # Info da amostra
        ax.text(0.5, 0.65, f'Amostra: {self.sample_id}',
                ha='center', va='center', fontsize=18,
                transform=ax.transAxes)

        # Classificação (se disponível)
        if 'classification' in self.report_data:
            classification = self.report_data['classification']
            subtype = classification['classification']['subtype']
            confidence = classification['classification']['confidence']

            ax.text(0.5, 0.55, 'Classificação:',
                    ha='center', va='center', fontsize=14, fontweight='bold',
                    transform=ax.transAxes)

            ax.text(0.5, 0.50, subtype,
                    ha='center', va='center', fontsize=16, color='#e74c3c',
                    transform=ax.transAxes, fontweight='bold')

            ax.text(0.5, 0.45, f'Confiança: {confidence:.0%}',
                    ha='center', va='center', fontsize=12,
                    transform=ax.transAxes)

            # Risco
            risk_group = classification['risk_stratification']['group']
            risk_color = {
                'Favorable': '#27ae60',
                'Intermediate': '#f39c12',
                'Adverse': '#e74c3c'
            }.get(risk_group, '#95a5a6')

            ax.text(0.5, 0.35, f'Grupo de Risco: {risk_group}',
                    ha='center', va='center', fontsize=14,
                    transform=ax.transAxes, color=risk_color, fontweight='bold')

        # Rodapé
        ax.text(0.5, 0.15, 'VKaryo - Virtual Karyotyping Workflow',
                ha='center', va='center', fontsize=12,
                transform=ax.transAxes, style='italic')

        ax.text(0.5, 0.10, f'Data: {datetime.now().strftime("%d/%m/%Y %H:%M")}',
                ha='center', va='center', fontsize=10,
                transform=ax.transAxes, color='gray')

        ax.text(0.5, 0.05, 'Este relatório foi gerado automaticamente',
                ha='center', va='center', fontsize=8,
                transform=ax.transAxes, color='gray', style='italic')

        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

    def create_summary_page(self, pdf):
        """
        Página de sumário executivo
        """
        fig = plt.figure(figsize=(8.5, 11))
        ax = fig.add_subplot(111)
        ax.axis('off')

        y_pos = 0.95

        # Título
        ax.text(0.5, y_pos, 'SUMÁRIO EXECUTIVO',
                ha='center', va='top', fontsize=18, fontweight='bold',
                transform=ax.transAxes)

        y_pos -= 0.08

        # Classificação
        if 'classification' in self.report_data:
            classification = self.report_data['classification']

            ax.text(0.05, y_pos, '1. CLASSIFICAÇÃO DIAGNÓSTICA',
                    ha='left', va='top', fontsize=14, fontweight='bold',
                    transform=ax.transAxes)
            y_pos -= 0.04

            ax.text(0.05, y_pos,
                    f"Subtipo: {classification['classification']['subtype']}",
                    ha='left', va='top', fontsize=11,
                    transform=ax.transAxes)
            y_pos -= 0.03

            ax.text(0.05, y_pos,
                    f"Classificação WHO: {classification['classification']['classification_who']}",
                    ha='left', va='top', fontsize=10,
                    transform=ax.transAxes)
            y_pos -= 0.03

            ax.text(0.05, y_pos,
                    f"Confiança: {classification['classification']['confidence']:.0%}",
                    ha='left', va='top', fontsize=10,
                    transform=ax.transAxes)

            y_pos -= 0.05

            # Evidências
            ax.text(0.05, y_pos, 'Evidências principais:',
                    ha='left', va='top', fontsize=10, fontweight='bold',
                    transform=ax.transAxes)
            y_pos -= 0.03

            for evidence in classification['classification']['evidence_summary'][:5]:
                ax.text(0.08, y_pos, f"• {evidence}",
                        ha='left', va='top', fontsize=9,
                        transform=ax.transAxes)
                y_pos -= 0.025

            y_pos -= 0.04

            # Risco
            ax.text(0.05, y_pos, '2. ESTRATIFICAÇÃO DE RISCO',
                    ha='left', va='top', fontsize=14, fontweight='bold',
                    transform=ax.transAxes)
            y_pos -= 0.04

            risk = classification['risk_stratification']

            ax.text(0.05, y_pos,
                    f"Grupo: {risk['group']} (Score: {risk['score']})",
                    ha='left', va='top', fontsize=11,
                    transform=ax.transAxes)
            y_pos -= 0.03

            ax.text(0.05, y_pos, 'Fatores prognósticos:',
                    ha='left', va='top', fontsize=10, fontweight='bold',
                    transform=ax.transAxes)
            y_pos -= 0.03

            for factor in risk['factors'][:6]:
                ax.text(0.08, y_pos, f"• {factor}",
                        ha='left', va='top', fontsize=9,
                        transform=ax.transAxes)
                y_pos -= 0.025

            y_pos -= 0.04

            # Recomendações
            ax.text(0.05, y_pos, '3. RECOMENDAÇÕES TERAPÊUTICAS',
                    ha='left', va='top', fontsize=14, fontweight='bold',
                    transform=ax.transAxes)
            y_pos -= 0.04

            for rec in classification['recommendations'][:5]:
                ax.text(0.05, y_pos,
                        f"[{rec['category']}] {rec['recommendation']}",
                        ha='left', va='top', fontsize=9,
                        transform=ax.transAxes)
                y_pos -= 0.02

                ax.text(0.08, y_pos,
                        f"Indicação: {rec['indication']}",
                        ha='left', va='top', fontsize=8, style='italic',
                        transform=ax.transAxes, color='gray')
                y_pos -= 0.03

        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

    def create_fusions_page(self, pdf):
        """
        Página de fusões gênicas
        """
        if 'fusions' not in self.report_data:
            return

        fusions = self.report_data['fusions']

        fig = plt.figure(figsize=(8.5, 11))
        ax = fig.add_subplot(111)
        ax.axis('off')

        y_pos = 0.95

        # Título
        ax.text(0.5, y_pos, 'FUSÕES GÊNICAS DETECTADAS',
                ha='center', va='top', fontsize=16, fontweight='bold',
                transform=ax.transAxes)

        y_pos -= 0.06

        ax.text(0.05, y_pos, f'Total de fusões: {len(fusions)}',
                ha='left', va='top', fontsize=11,
                transform=ax.transAxes)

        y_pos -= 0.05

        # Lista fusões (top 15)
        for idx, (_, fusion) in enumerate(fusions.head(15).iterrows(), 1):
            if y_pos < 0.15:
                break

            fusion_name = fusion.get('fusion', f"{fusion['gene1']}--{fusion['gene2']}")
            confidence = fusion.get('confidence', 'N/A')

            # Nome da fusão
            ax.text(0.05, y_pos, f"{idx}. {fusion_name}",
                    ha='left', va='top', fontsize=10, fontweight='bold',
                    transform=ax.transAxes)
            y_pos -= 0.025

            # Tipo e confiança
            leuk_type = fusion.get('annot_leukemia_type', 'N/A')
            ax.text(0.08, y_pos,
                    f"Tipo: {leuk_type} | Confiança: {confidence}",
                    ha='left', va='top', fontsize=8,
                    transform=ax.transAxes)
            y_pos -= 0.02

            # Significância clínica
            clinical_sig = fusion.get('annot_clinical_significance', '')
            if clinical_sig and clinical_sig != 'N/A':
                ax.text(0.08, y_pos, f"Significância: {clinical_sig}",
                        ha='left', va='top', fontsize=8, style='italic',
                        transform=ax.transAxes, color='#e74c3c')
                y_pos -= 0.02

            # Prognóstico
            prognosis = fusion.get('annot_prognosis', '')
            if prognosis and prognosis != 'N/A' and prognosis != 'Desconhecido':
                ax.text(0.08, y_pos, f"Prognóstico: {prognosis}",
                        ha='left', va='top', fontsize=8,
                        transform=ax.transAxes)
                y_pos -= 0.02

            # Terapia
            therapy = fusion.get('annot_therapy', '')
            if therapy and therapy != 'N/A' and therapy != 'Protocolo padrão baseado em tipo celular':
                ax.text(0.08, y_pos, f"Terapia: {therapy}",
                        ha='left', va='top', fontsize=7,
                        transform=ax.transAxes, color='#27ae60')
                y_pos -= 0.02

            y_pos -= 0.015  # Espaço entre fusões

        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

    def create_cnv_summary_page(self, pdf):
        """
        Página de sumário de CNVs
        """
        if 'cnvs' not in self.report_data:
            return

        cnvs = self.report_data['cnvs']

        fig = plt.figure(figsize=(8.5, 11))
        ax = fig.add_subplot(111)
        ax.axis('off')

        y_pos = 0.95

        # Título
        ax.text(0.5, y_pos, 'ALTERAÇÕES NO NÚMERO DE CÓPIAS (CNVs)',
                ha='center', va='top', fontsize=16, fontweight='bold',
                transform=ax.transAxes)

        y_pos -= 0.06

        # Estatísticas
        total = len(cnvs)
        validated = cnvs[cnvs['validation_status'] == 'VALIDATED']
        gains = cnvs[cnvs['cnv_type'] == 'gain']
        losses = cnvs[cnvs['cnv_type'] == 'loss']

        ax.text(0.05, y_pos, f'Total de CNVs: {total}',
                ha='left', va='top', fontsize=11,
                transform=ax.transAxes)
        y_pos -= 0.03

        ax.text(0.05, y_pos, f'Validados por expressão: {len(validated)} ({len(validated)/total*100:.1f}%)',
                ha='left', va='top', fontsize=10,
                transform=ax.transAxes)
        y_pos -= 0.03

        ax.text(0.05, y_pos, f'Ganhos: {len(gains)} | Perdas: {len(losses)}',
                ha='left', va='top', fontsize=10,
                transform=ax.transAxes)

        y_pos -= 0.05

        # CNVs validados por cromossomo
        ax.text(0.05, y_pos, 'CNVs validados por cromossomo:',
                ha='left', va='top', fontsize=11, fontweight='bold',
                transform=ax.transAxes)
        y_pos -= 0.03

        chr_summary = validated.groupby(['chrom', 'cnv_type']).agg({
            'n_genes': 'sum',
            'log2': 'mean'
        }).reset_index()

        for chrom in sorted(validated['chrom'].unique()):
            chr_data = chr_summary[chr_summary['chrom'] == chrom]

            gains_chr = chr_data[chr_data['cnv_type'] == 'gain']
            losses_chr = chr_data[chr_data['cnv_type'] == 'loss']

            line = f"{chrom}: "
            parts = []

            if len(gains_chr) > 0:
                n_genes = gains_chr['n_genes'].sum()
                parts.append(f"+{int(n_genes)} genes")

            if len(losses_chr) > 0:
                n_genes = losses_chr['n_genes'].sum()
                parts.append(f"-{int(n_genes)} genes")

            line += " | ".join(parts)

            ax.text(0.08, y_pos, line,
                    ha='left', va='top', fontsize=9,
                    transform=ax.transAxes)
            y_pos -= 0.025

        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

    def embed_ideogram(self, pdf, ideogram_file):
        """
        Inclui ideograma de aneuploidias
        """
        if not ideogram_file:
            return

        logger.info(f"Incluindo ideograma: {ideogram_file}")

        # Carrega e inclui PDF do ideograma
        try:
            from PyPDF2 import PdfReader, PdfWriter

            # Adiciona páginas do ideograma ao relatório
            reader = PdfReader(ideogram_file)
            for page in reader.pages:
                pdf.savePage(page)

        except ImportError:
            logger.warning("PyPDF2 não instalado, incluindo como figura")

            # Alternativa: carrega como imagem
            from PIL import Image
            import numpy as np

            try:
                # Converte PDF para imagem (requer pdf2image)
                from pdf2image import convert_from_path

                images = convert_from_path(ideogram_file)

                for img in images:
                    fig = plt.figure(figsize=(8.5, 11))
                    ax = fig.add_subplot(111)
                    ax.imshow(np.array(img))
                    ax.axis('off')

                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close()

            except ImportError:
                logger.warning("pdf2image não instalado, pulando ideograma")

    def generate_report(self, output_file, ideogram_file=None):
        """
        Gera relatório PDF completo
        """
        logger.info("Gerando relatório PDF...")

        with PdfPages(output_file) as pdf:
            # Página de capa
            self.create_cover_page(pdf)

            # Sumário executivo
            self.create_summary_page(pdf)

            # Fusões
            self.create_fusions_page(pdf)

            # CNVs
            self.create_cnv_summary_page(pdf)

            # Ideograma (se disponível)
            # Note: incluir PDFs em PDFs requer PyPDF2
            # Por simplicidade, vamos criar página de referência
            if ideogram_file:
                fig = plt.figure(figsize=(8.5, 11))
                ax = fig.add_subplot(111)
                ax.axis('off')

                ax.text(0.5, 0.5,
                        'IDEOGRAMA DE ANEUPLOIDIAS\n\nVer arquivo separado:\n' + ideogram_file,
                        ha='center', va='center', fontsize=14,
                        transform=ax.transAxes)

                pdf.savefig(fig, bbox_inches='tight')
                plt.close()

            # Metadados
            d = pdf.infodict()
            d['Title'] = f'Relatório de Cariotipagem Virtual - {self.sample_id}'
            d['Author'] = 'VKaryo Workflow'
            d['Subject'] = 'Virtual Karyotyping Report'
            d['Keywords'] = 'Leukemia, Karyotype, NGS, CNV, Fusion'
            d['CreationDate'] = datetime.now()

        logger.info(f"✓ Relatório PDF salvo: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Gera relatório PDF completo de cariotipagem virtual'
    )
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--classification', required=True,
                       help='Arquivo JSON de classificação')
    parser.add_argument('--fusions', help='Fusões anotadas TSV')
    parser.add_argument('--cnvs', help='CNVs validados TSV')
    parser.add_argument('--ideogram', help='Ideograma PDF')
    parser.add_argument('--output', required=True, help='Arquivo PDF de saída')

    args = parser.parse_args()

    # Inicializa gerador
    generator = PDFReportGenerator(args.sample)

    # Carrega dados
    generator.load_classification(args.classification)

    if args.fusions:
        generator.load_fusions(args.fusions)

    if args.cnvs:
        generator.load_cnvs(args.cnvs)

    # Gera relatório
    generator.generate_report(args.output, args.ideogram)

    logger.info("\n" + "="*60)
    logger.info("RELATÓRIO PDF GERADO COM SUCESSO!")
    logger.info("="*60)
    logger.info(f"\nOutput: {args.output}")


if __name__ == '__main__':
    main()
