#!/usr/bin/env python3
"""
Gera relatório clínico consolidado integrando todos os resultados VKaryo

Integra:
- CNVs (CNVkit)
- SVs (Manta + Delly)
- SNVs (Strelka + ANNOVAR)
- Fusões gênicas (consolidadas)
- Deconvolução celular (se disponível)

Uso:
    python generate_clinical_report.py \
        --sample SAMPLE_ID \
        --output-dir output_vkaryo/ \
        --output report_SAMPLE_ID.pdf
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ClinicalReportGenerator:
    """
    Gera relatório clínico integrado para cariotipagem virtual
    """

    def __init__(self, sample_id, output_dir):
        self.sample_id = sample_id
        self.output_dir = Path(output_dir)
        self.report_data = {}

    def load_cnvs(self):
        """
        Carrega e processa CNVs do CNVkit
        """
        logger.info("Loading CNV data")

        try:
            cns_file = self.output_dir / f"cnvkit/{self.sample_id}_markduplicates_recal.cns"
            cnvs = pd.read_csv(cns_file, sep='\t', comment='#')

            # Filtra alterações significativas
            cnvs['status'] = 'neutral'
            cnvs.loc[cnvs['log2'] > 0.3, 'status'] = 'gain'
            cnvs.loc[cnvs['log2'] < -0.3, 'status'] = 'loss'

            altered = cnvs[cnvs['status'] != 'neutral']

            # Identifica alterações cromossômicas completas
            chr_alterations = self._detect_chromosome_alterations(cnvs)

            # Calcula aneuploidia
            total_bases = (cnvs['end'] - cnvs['start']).sum()
            altered_bases = (altered['end'] - altered['start']).sum()
            aneuploidy_score = altered_bases / total_bases

            self.report_data['cnvs'] = {
                'total_cnvs': len(altered),
                'gains': (cnvs['status'] == 'gain').sum(),
                'losses': (cnvs['status'] == 'loss').sum(),
                'chromosome_alterations': chr_alterations,
                'aneuploidy_score': aneuploidy_score,
                'data': altered
            }

            logger.info(f"  Found {len(altered)} significant CNVs")

        except Exception as e:
            logger.error(f"Error loading CNVs: {e}")
            self.report_data['cnvs'] = None

    def _detect_chromosome_alterations(self, cnvs):
        """
        Detecta ganhos/perdas cromossômicas completas
        """
        chr_status = {}

        for chrom in cnvs['chromosome'].unique():
            chr_cnvs = cnvs[cnvs['chromosome'] == chrom]

            # Calcula proporção do cromossomo alterado
            chr_length = (chr_cnvs['end'] - chr_cnvs['start']).sum()
            gained = chr_cnvs[chr_cnvs['log2'] > 0.3]
            lost = chr_cnvs[chr_cnvs['log2'] < -0.3]

            gained_prop = (gained['end'] - gained['start']).sum() / chr_length
            lost_prop = (lost['end'] - lost['start']).sum() / chr_length

            # Se >70% do cromossomo está alterado
            if gained_prop > 0.7:
                chr_status[chrom] = 'trisomy' if gained_prop > 0.9 else 'partial_gain'
            elif lost_prop > 0.7:
                chr_status[chrom] = 'monosomy' if lost_prop > 0.9 else 'partial_loss'

        return chr_status

    def load_svs(self):
        """
        Carrega variantes estruturais
        """
        logger.info("Loading structural variants")

        try:
            # Prioriza resultado anotado (SANSA ou AnnotSV)
            sv_file = self.output_dir / f"sansa/manta/{self.sample_id}_manta_sansa.tsv.gz"

            if not sv_file.exists():
                sv_file = self.output_dir / f"annotsv/manta/{self.sample_id}_annotsv_manta.tsv"

            if sv_file.exists():
                svs = pd.read_csv(sv_file, sep='\t')

                # Filtra por confiança
                if 'FILTER' in svs.columns:
                    svs_pass = svs[svs['FILTER'] == 'PASS']
                else:
                    svs_pass = svs

                self.report_data['svs'] = {
                    'total': len(svs_pass),
                    'data': svs_pass
                }

                logger.info(f"  Found {len(svs_pass)} structural variants")
            else:
                logger.warning("  No SV file found")
                self.report_data['svs'] = None

        except Exception as e:
            logger.error(f"Error loading SVs: {e}")
            self.report_data['svs'] = None

    def load_snvs(self):
        """
        Carrega SNVs e indels anotados
        """
        logger.info("Loading SNVs/indels")

        try:
            snv_file = self.output_dir / f"strelka/annovar/{self.sample_id}_annotation_updated.hg38_multianno.csv"

            if snv_file.exists():
                snvs = pd.read_csv(snv_file)

                # Filtra variantes clinicamente relevantes
                relevant_genes = [
                    'FLT3', 'NPM1', 'CEBPA', 'IDH1', 'IDH2', 'TET2', 'DNMT3A',
                    'TP53', 'ASXL1', 'RUNX1', 'KIT', 'JAK2', 'CALR', 'MPL',
                    'SF3B1', 'SRSF2', 'U2AF1', 'ZRSR2', 'EZH2', 'NRAS', 'KRAS'
                ]

                if 'Gene.refGene' in snvs.columns:
                    clinically_relevant = snvs[snvs['Gene.refGene'].isin(relevant_genes)]
                else:
                    clinically_relevant = snvs

                # Calcula TMB (tumor mutational burden)
                tmb = len(snvs) / 30  # Assumindo ~30Mb de exoma

                self.report_data['snvs'] = {
                    'total': len(snvs),
                    'clinically_relevant': len(clinically_relevant),
                    'tmb': tmb,
                    'data': clinically_relevant
                }

                logger.info(f"  Found {len(clinically_relevant)} clinically relevant variants")
            else:
                logger.warning("  No SNV file found")
                self.report_data['snvs'] = None

        except Exception as e:
            logger.error(f"Error loading SNVs: {e}")
            self.report_data['snvs'] = None

    def load_fusions(self):
        """
        Carrega fusões gênicas consolidadas
        """
        logger.info("Loading gene fusions")

        try:
            fusion_file = self.output_dir / f"consolidated_fusions/{self.sample_id}_consolidated_fusions.tsv"

            if fusion_file.exists():
                fusions = pd.read_csv(fusion_file, sep='\t')

                # Filtra alta confiança
                high_conf = fusions[fusions['consolidated_confidence'].isin(['HIGH', 'VERY_HIGH'])]

                self.report_data['fusions'] = {
                    'total': len(fusions),
                    'high_confidence': len(high_conf),
                    'data': high_conf if len(high_conf) > 0 else fusions
                }

                logger.info(f"  Found {len(high_conf)} high-confidence fusions")
            else:
                logger.warning("  No fusion file found")
                self.report_data['fusions'] = None

        except Exception as e:
            logger.error(f"Error loading fusions: {e}")
            self.report_data['fusions'] = None

    def load_cell_deconvolution(self):
        """
        Carrega resultados de deconvolução celular
        """
        logger.info("Loading cell deconvolution data")

        try:
            deconv_dir = self.output_dir / "deconvolution"

            # Cell proportions
            props_file = deconv_dir / f"{self.sample_id}_cell_proportions.csv"
            if props_file.exists():
                proportions = pd.read_csv(props_file, index_col=0)

                # Classification
                class_file = deconv_dir / f"{self.sample_id}_classification.csv"
                classification = pd.read_csv(class_file, index_col=0) if class_file.exists() else None

                # MRD
                mrd_file = deconv_dir / f"{self.sample_id}_mrd_assessment.csv"
                mrd = pd.read_csv(mrd_file, index_col=0) if mrd_file.exists() else None

                self.report_data['deconvolution'] = {
                    'proportions': proportions,
                    'classification': classification,
                    'mrd': mrd
                }

                logger.info("  Deconvolution data loaded")
            else:
                logger.warning("  No deconvolution data found")
                self.report_data['deconvolution'] = None

        except Exception as e:
            logger.error(f"Error loading deconvolution: {e}")
            self.report_data['deconvolution'] = None

    def classify_leukemia(self):
        """
        Classifica tipo de leucemia baseado em todos os dados
        """
        logger.info("Classifying leukemia subtype")

        classification = {
            'subtype': 'Unclassified',
            'confidence': 0.0,
            'evidence': []
        }

        # Prioridade 1: Fusões patognomônicas
        if self.report_data.get('fusions'):
            for _, fusion in self.report_data['fusions']['data'].iterrows():
                fusion_name = f"{fusion['gene1']}-{fusion['gene2']}"

                if 'BCR-ABL1' in fusion_name or 'ABL1-BCR' in fusion_name:
                    classification['subtype'] = 'CML / Ph+ ALL'
                    classification['confidence'] = 0.95
                    classification['evidence'].append('BCR-ABL1 fusion detected')
                    return classification

                if 'PML-RARA' in fusion_name:
                    classification['subtype'] = 'APL (M3)'
                    classification['confidence'] = 0.95
                    classification['evidence'].append('PML-RARA fusion detected')
                    return classification

                if 'RUNX1-RUNX1T1' in fusion_name:
                    classification['subtype'] = 'AML with t(8;21)'
                    classification['confidence'] = 0.90
                    classification['evidence'].append('RUNX1-RUNX1T1 fusion detected')
                    return classification

        # Prioridade 2: Deconvolução celular
        if self.report_data.get('deconvolution'):
            deconv_class = self.report_data['deconvolution'].get('classification')
            if deconv_class is not None:
                classification['subtype'] = deconv_class.iloc[0]['subtype']
                classification['confidence'] = deconv_class.iloc[0]['confidence']
                classification['evidence'].append('Cell type deconvolution')

        # Prioridade 3: Mutações específicas
        if self.report_data.get('snvs'):
            snvs = self.report_data['snvs']['data']

            if 'Gene.refGene' in snvs.columns:
                genes_mutated = snvs['Gene.refGene'].unique()

                if 'NPM1' in genes_mutated and 'FLT3' not in genes_mutated:
                    classification['subtype'] = 'AML with mutated NPM1'
                    classification['confidence'] = 0.70
                    classification['evidence'].append('NPM1 mutation without FLT3-ITD')

                if 'JAK2' in genes_mutated:
                    classification['subtype'] = 'Myeloproliferative neoplasm'
                    classification['confidence'] = 0.65
                    classification['evidence'].append('JAK2 mutation')

        # Prioridade 4: Alterações cromossômicas
        if self.report_data.get('cnvs'):
            chr_alt = self.report_data['cnvs']['chromosome_alterations']

            if 'chr8' in chr_alt and chr_alt['chr8'] == 'trisomy':
                classification['subtype'] = 'AML (trisomy 8 associated)'
                classification['confidence'] = 0.50
                classification['evidence'].append('Trisomy 8 detected')

        return classification

    def assess_risk(self, classification):
        """
        Avalia grupo de risco (ELN 2022 para AML, ou equivalente)
        """
        logger.info("Assessing risk stratification")

        risk = {
            'group': 'Intermediate',
            'score': 0,
            'factors': []
        }

        # Fatores favoráveis
        if classification['subtype'] == 'AML with t(8;21)':
            risk['score'] -= 2
            risk['factors'].append('Favorable: t(8;21)')

        if classification['subtype'] == 'APL (M3)':
            risk['score'] -= 3
            risk['factors'].append('Favorable: APL (excellent prognosis with ATRA)')

        if self.report_data.get('snvs'):
            snvs = self.report_data['snvs']['data']
            if 'Gene.refGene' in snvs.columns:
                genes = snvs['Gene.refGene'].unique()

                if 'NPM1' in genes and 'FLT3' not in genes:
                    risk['score'] -= 1
                    risk['factors'].append('Favorable: NPM1 mutation without FLT3-ITD')

                if 'CEBPA' in genes:
                    risk['score'] -= 1
                    risk['factors'].append('Favorable: CEBPA mutation')

        # Fatores adversos
        if self.report_data.get('snvs'):
            snvs = self.report_data['snvs']['data']
            if 'Gene.refGene' in snvs.columns:
                genes = snvs['Gene.refGene'].unique()

                if 'TP53' in genes:
                    risk['score'] += 3
                    risk['factors'].append('Adverse: TP53 mutation')

                if 'ASXL1' in genes or 'RUNX1' in genes:
                    risk['score'] += 1
                    risk['factors'].append('Adverse: ASXL1/RUNX1 mutation')

        if self.report_data.get('cnvs'):
            chr_alt = self.report_data['cnvs']['chromosome_alterations']

            if 'chr7' in chr_alt and chr_alt['chr7'] in ['monosomy', 'partial_loss']:
                risk['score'] += 2
                risk['factors'].append('Adverse: -7/del(7q)')

            if 'chr5' in chr_alt and chr_alt['chr5'] in ['monosomy', 'partial_loss']:
                risk['score'] += 2
                risk['factors'].append('Adverse: -5/del(5q)')

        # Aneuplodia complexa
        if self.report_data.get('cnvs'):
            if self.report_data['cnvs']['aneuploidy_score'] > 0.3:
                risk['score'] += 1
                risk['factors'].append('Adverse: Complex karyotype')

        # Classificação final
        if risk['score'] <= -2:
            risk['group'] = 'Favorable'
        elif risk['score'] >= 2:
            risk['group'] = 'Adverse'
        else:
            risk['group'] = 'Intermediate'

        return risk

    def generate_text_report(self, output_file):
        """
        Gera relatório em texto
        """
        logger.info("Generating text report")

        # Carrega todos os dados
        self.load_cnvs()
        self.load_svs()
        self.load_snvs()
        self.load_fusions()
        self.load_cell_deconvolution()

        # Classificação
        classification = self.classify_leukemia()
        risk = self.assess_risk(classification)

        # Gera relatório
        report = f"""
╔══════════════════════════════════════════════════════════════════════════╗
║                  RELATÓRIO DE CARIOTIPAGEM VIRTUAL                       ║
║                            VKaryo v2.0                                   ║
╚══════════════════════════════════════════════════════════════════════════╝

Amostra: {self.sample_id}
Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. CLASSIFICAÇÃO DIAGNÓSTICA
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Subtipo: {classification['subtype']}
Confiança: {classification['confidence']:.1%}

Evidências:
"""

        for evidence in classification['evidence']:
            report += f"  • {evidence}\n"

        report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2. ESTRATIFICAÇÃO DE RISCO
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Grupo de Risco: {risk['group']}
Score: {risk['score']}

Fatores Prognósticos:
"""

        for factor in risk['factors']:
            report += f"  • {factor}\n"

        # CNVs
        if self.report_data.get('cnvs'):
            cnv_data = self.report_data['cnvs']
            report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
3. ALTERAÇÕES NO NÚMERO DE CÓPIAS (CNVs)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Total de CNVs: {cnv_data['total_cnvs']}
  Ganhos: {cnv_data['gains']}
  Perdas: {cnv_data['losses']}

Score de Aneuploidia: {cnv_data['aneuploidy_score']:.2%}

Alterações Cromossômicas:
"""
            if cnv_data['chromosome_alterations']:
                for chrom, status in cnv_data['chromosome_alterations'].items():
                    report += f"  • {chrom}: {status}\n"
            else:
                report += "  Nenhuma alteração cromossômica completa detectada\n"

        # SVs
        if self.report_data.get('svs'):
            sv_data = self.report_data['svs']
            report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
4. VARIANTES ESTRUTURAIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Total de SVs: {sv_data['total']}
"""

        # Fusões
        if self.report_data.get('fusions'):
            fusion_data = self.report_data['fusions']
            report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
5. FUSÕES GÊNICAS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Total de fusões: {fusion_data['total']}
Alta confiança: {fusion_data['high_confidence']}

Fusões detectadas:
"""
            for _, fusion in fusion_data['data'].head(10).iterrows():
                report += f"  • {fusion['gene1']}-{fusion['gene2']}"
                report += f" (suporte: {fusion.get('support_count', 'N/A')} callers"
                report += f", confiança: {fusion.get('consolidated_confidence', 'N/A')})\n"

                if 'clinical_significance' in fusion and fusion['clinical_significance']:
                    report += f"    Significância: {fusion.get('leukemia_type', 'Unknown')}\n"

        # SNVs
        if self.report_data.get('snvs'):
            snv_data = self.report_data['snvs']
            report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
6. VARIANTES DE NUCLEOTÍDEO ÚNICO (SNVs/Indels)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Total de variantes: {snv_data['total']}
Clinicamente relevantes: {snv_data['clinically_relevant']}
TMB (Tumor Mutational Burden): {snv_data['tmb']:.2f} mut/Mb

Mutações em genes driver:
"""
            for _, var in snv_data['data'].head(20).iterrows():
                gene = var.get('Gene.refGene', 'Unknown')
                mutation = f"{var.get('AAChange.refGene', 'N/A')}"
                report += f"  • {gene}: {mutation}\n"

        # Deconvolução
        if self.report_data.get('deconvolution'):
            deconv_data = self.report_data['deconvolution']
            report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
7. COMPOSIÇÃO CELULAR (Deconvolução)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Principais tipos celulares:
"""
            props = deconv_data['proportions'].iloc[0].sort_values(ascending=False).head(10)
            for cell_type, proportion in props.items():
                if proportion > 0.01:
                    report += f"  • {cell_type}: {proportion:.1%}\n"

            if deconv_data.get('mrd'):
                mrd = deconv_data['mrd'].iloc[0]
                report += f"\nDoença Residual Mínima (MRD):\n"
                report += f"  Status: {mrd['mrd_status']}\n"
                report += f"  Proporção de blastos: {mrd['blast_proportion']:.3%}\n"
                report += f"  Nível de risco: {mrd['risk_level']}\n"

        # Recomendações
        report += f"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
8. RECOMENDAÇÕES TERAPÊUTICAS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

        recommendations = self._generate_recommendations(classification, risk)
        for rec in recommendations:
            report += f"  • {rec}\n"

        report += "\n" + "="*78 + "\n"
        report += "Este relatório foi gerado automaticamente pelo VKaryo.\n"
        report += "Interpretação clínica deve ser feita por profissional qualificado.\n"
        report += "="*78 + "\n"

        # Salva
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(report)

        logger.info(f"✓ Report saved: {output_file}")

        return report

    def _generate_recommendations(self, classification, risk):
        """
        Gera recomendações terapêuticas
        """
        recommendations = []

        # Baseado no subtipo
        if 'BCR-ABL1' in classification['subtype']:
            recommendations.append("Inibidores de tirosina quinase (TKIs): Imatinib, Dasatinib, ou Nilotinib")
            recommendations.append("Monitoramento molecular BCR-ABL1 a cada 3 meses")

        if 'APL' in classification['subtype']:
            recommendations.append("ATRA (All-trans retinoic acid) + Arsenic trioxide")
            recommendations.append("Prognóstico excelente com tratamento adequado (>90% cura)")

        # Baseado em risco
        if risk['group'] == 'Favorable':
            recommendations.append("Quimioterapia de indução padrão")
            recommendations.append("Consolidação com altas doses de citarabina")

        elif risk['group'] == 'Adverse':
            recommendations.append("Quimioterapia intensiva")
            recommendations.append("Considerar transplante alogênico de células tronco em primeira remissão")

        # Baseado em mutações específicas
        if self.report_data.get('snvs'):
            snvs = self.report_data['snvs']['data']
            if 'Gene.refGene' in snvs.columns:
                genes = snvs['Gene.refGene'].unique()

                if 'FLT3' in genes:
                    recommendations.append("Considerar inibidores de FLT3 (Midostaurin, Gilteritinib)")

                if 'IDH1' in genes or 'IDH2' in genes:
                    recommendations.append("Considerar inibidores de IDH (Ivosidenib, Enasidenib)")

        if not recommendations:
            recommendations.append("Seguir protocolo padrão baseado no subtipo e risco")

        return recommendations


def main():
    parser = argparse.ArgumentParser(
        description='Generate consolidated clinical report for virtual karyotyping'
    )
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output-dir', required=True, help='VKaryo output directory')
    parser.add_argument('--output', required=True, help='Output report file')

    args = parser.parse_args()

    generator = ClinicalReportGenerator(args.sample, args.output_dir)
    report = generator.generate_text_report(args.output)

    print("\n" + "="*78)
    print("RELATÓRIO GERADO COM SUCESSO!")
    print("="*78)
    print(f"\nArquivo: {args.output}")
    print("\nResumo:")
    print(report[:500] + "...")


if __name__ == '__main__':
    main()
