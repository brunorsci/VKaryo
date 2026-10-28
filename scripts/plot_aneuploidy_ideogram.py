#!/usr/bin/env python3
"""
Gera ideograma de aneuploidias (genoma completo)
Visualiza CNVs, alterações cromossômicas e aneuploidias

Uso:
    python plot_aneuploidy_ideogram.py \
        --cnv output_vkaryo/cnvkit/sample.cns \
        --validated output_vkaryo/cnv_validated/sample_cnv_validated.tsv \
        --sample SAMPLE_ID \
        --output output_vkaryo/plots/sample_aneuploidy_ideogram.pdf
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


# Comprimentos cromossômicos (GRCh38 - bp)
CHR_LENGTHS = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415,
}

# Cores
COLOR_GAIN = '#e74c3c'      # Vermelho
COLOR_LOSS = '#3498db'      # Azul
COLOR_NEUTRAL = '#ecf0f1'   # Cinza claro
COLOR_HIGH_CONF = '#27ae60' # Verde (validado)


class AneuploidyPlotter:
    """
    Plota ideograma de aneuploidias genoma-completo
    """

    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.chr_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

    def load_cnvs(self, cnv_file):
        """Carrega CNVs"""
        logger.info(f"Carregando CNVs: {cnv_file}")

        cnvs = pd.read_csv(cnv_file, sep='\t', comment='#')

        if 'chromosome' in cnvs.columns:
            cnvs.rename(columns={'chromosome': 'chrom'}, inplace=True)

        # Classifica
        cnvs['cnv_type'] = 'neutral'
        cnvs.loc[cnvs['log2'] > 0.3, 'cnv_type'] = 'gain'
        cnvs.loc[cnvs['log2'] < -0.3, 'cnv_type'] = 'loss'

        logger.info(f"  Total: {len(cnvs)}, Alterados: {(cnvs['cnv_type'] != 'neutral').sum()}")

        return cnvs

    def load_validated(self, validated_file):
        """Carrega CNVs validados (opcional)"""
        if not validated_file:
            return None

        logger.info(f"Carregando CNVs validados: {validated_file}")

        validated = pd.read_csv(validated_file, sep='\t')

        logger.info(f"  Validados: {len(validated)}")

        return validated

    def calculate_chromosome_aneuploidy(self, cnvs):
        """
        Calcula se cromossomo tem aneuplodia (>50% alterado)
        """
        chr_status = {}

        for chrom in self.chr_order:
            if chrom not in CHR_LENGTHS:
                continue

            chr_cnvs = cnvs[cnvs['chrom'] == chrom]

            if len(chr_cnvs) == 0:
                chr_status[chrom] = 'neutral'
                continue

            chr_length = CHR_LENGTHS[chrom]

            # Calcula proporção alterada
            gained = chr_cnvs[chr_cnvs['cnv_type'] == 'gain']
            lost = chr_cnvs[chr_cnvs['cnv_type'] == 'loss']

            gained_bases = (gained['end'] - gained['start']).sum()
            lost_bases = (lost['end'] - lost['start']).sum()

            gained_prop = gained_bases / chr_length
            lost_prop = lost_bases / chr_length

            # Classifica
            if gained_prop > 0.7:
                chr_status[chrom] = 'trisomy' if gained_prop > 0.9 else 'partial_gain'
            elif lost_prop > 0.7:
                chr_status[chrom] = 'monosomy' if lost_prop > 0.9 else 'partial_loss'
            else:
                chr_status[chrom] = 'neutral'

        return chr_status

    def plot_genome_ideogram(self, cnvs, validated=None, output_file=None):
        """
        Plota ideograma genoma completo
        """
        logger.info("Gerando ideograma...")

        fig, axes = plt.subplots(
            nrows=24, ncols=1,
            figsize=(16, 20),
            gridspec_kw={'hspace': 0.3}
        )

        # Calcula aneuploidias
        chr_status = self.calculate_chromosome_aneuploidy(cnvs)

        # Plota cada cromossomo
        for idx, chrom in enumerate(self.chr_order):
            ax = axes[idx]

            if chrom not in CHR_LENGTHS:
                ax.axis('off')
                continue

            chr_length = CHR_LENGTHS[chrom]
            chr_cnvs = cnvs[cnvs['chrom'] == chrom]

            # Background (cromossomo neutro)
            ax.add_patch(Rectangle(
                (0, 0), chr_length / 1e6, 1,
                facecolor=COLOR_NEUTRAL, edgecolor='black', linewidth=0.5
            ))

            # Plota CNVs
            for _, cnv in chr_cnvs.iterrows():
                start_mb = cnv['start'] / 1e6
                width_mb = (cnv['end'] - cnv['start']) / 1e6

                if cnv['cnv_type'] == 'gain':
                    color = COLOR_GAIN
                    y_pos = 0
                    height = 1
                elif cnv['cnv_type'] == 'loss':
                    color = COLOR_LOSS
                    y_pos = 0
                    height = 1
                else:
                    continue

                # Marca CNVs validados com borda verde
                if validated is not None:
                    cnv_id = f"{cnv['chrom']}:{cnv['start']}-{cnv['end']}"
                    is_validated = any(
                        (validated['chrom'] == cnv['chrom']) &
                        (validated['start'] == cnv['start']) &
                        (validated['validation_status'] == 'VALIDATED')
                    )
                    edgecolor = COLOR_HIGH_CONF if is_validated else 'none'
                    linewidth = 2 if is_validated else 0
                else:
                    edgecolor = 'none'
                    linewidth = 0

                ax.add_patch(Rectangle(
                    (start_mb, y_pos), width_mb, height,
                    facecolor=color, alpha=0.8,
                    edgecolor=edgecolor, linewidth=linewidth
                ))

            # Configurações do eixo
            ax.set_xlim(0, chr_length / 1e6)
            ax.set_ylim(0, 1)
            ax.set_yticks([])

            # Label do cromossomo + status de aneuplodia
            status = chr_status.get(chrom, 'neutral')
            status_label = ''

            if status == 'trisomy':
                status_label = ' [TRISSOMIA]'
                ax.set_facecolor('#ffe6e6')  # Fundo vermelho claro
            elif status == 'monosomy':
                status_label = ' [MONOSSOMIA]'
                ax.set_facecolor('#e6f2ff')  # Fundo azul claro
            elif status == 'partial_gain':
                status_label = ' [Ganho parcial]'
            elif status == 'partial_loss':
                status_label = ' [Perda parcial]'

            chr_label = chrom.replace('chr', '')
            ax.set_ylabel(f'{chr_label}{status_label}', fontsize=10, fontweight='bold')

            # Eixo X apenas no último
            if idx == len(self.chr_order) - 1:
                ax.set_xlabel('Posição (Mb)', fontsize=12)
            else:
                ax.set_xticklabels([])

        # Título
        fig.suptitle(
            f'Ideograma de Aneuploidias - {self.sample_id}',
            fontsize=16, fontweight='bold', y=0.995
        )

        # Legenda
        legend_elements = [
            mpatches.Patch(facecolor=COLOR_GAIN, label='Ganho (log2 > 0.3)'),
            mpatches.Patch(facecolor=COLOR_LOSS, label='Perda (log2 < -0.3)'),
            mpatches.Patch(facecolor=COLOR_HIGH_CONF, edgecolor=COLOR_HIGH_CONF,
                          linewidth=2, fill=False, label='Validado por expressão'),
        ]

        fig.legend(
            handles=legend_elements,
            loc='upper right',
            bbox_to_anchor=(0.98, 0.99),
            fontsize=10
        )

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"✓ Ideograma salvo: {output_file}")

        return fig

    def plot_cnv_summary_circular(self, cnvs, output_file=None):
        """
        Plota sumário circular (tipo Circos)
        """
        logger.info("Gerando plot circular...")

        fig = plt.figure(figsize=(14, 14))
        ax = fig.add_subplot(111, projection='polar')

        # Calcula ângulos para cada cromossomo
        total_length = sum(CHR_LENGTHS.values())
        chr_angles = {}
        cumulative = 0

        for chrom in self.chr_order:
            if chrom not in CHR_LENGTHS:
                continue

            start_angle = (cumulative / total_length) * 2 * np.pi
            chr_length = CHR_LENGTHS[chrom]
            end_angle = ((cumulative + chr_length) / total_length) * 2 * np.pi

            chr_angles[chrom] = {
                'start': start_angle,
                'end': end_angle,
                'mid': (start_angle + end_angle) / 2
            }

            cumulative += chr_length

        # Plota cromossomos (fundo)
        for chrom, angles in chr_angles.items():
            theta = np.linspace(angles['start'], angles['end'], 100)
            r = np.ones_like(theta) * 1.0

            ax.plot(theta, r, color='black', linewidth=2)
            ax.fill_between(theta, 0.95, 1.05, alpha=0.1, color='gray')

            # Label
            chr_label = chrom.replace('chr', '')
            ax.text(
                angles['mid'], 1.15, chr_label,
                ha='center', va='center', fontsize=10, fontweight='bold'
            )

        # Plota CNVs
        for _, cnv in cnvs.iterrows():
            chrom = cnv['chrom']

            if chrom not in chr_angles:
                continue

            chr_info = chr_angles[chrom]
            chr_length = CHR_LENGTHS[chrom]

            # Posição angular do CNV
            start_frac = cnv['start'] / chr_length
            end_frac = cnv['end'] / chr_length

            theta_start = chr_info['start'] + start_frac * (chr_info['end'] - chr_info['start'])
            theta_end = chr_info['start'] + end_frac * (chr_info['end'] - chr_info['start'])

            theta_cnv = np.linspace(theta_start, theta_end, 50)

            if cnv['cnv_type'] == 'gain':
                r_cnv = np.ones_like(theta_cnv) * (1.0 + cnv['log2'] * 0.3)
                color = COLOR_GAIN
            elif cnv['cnv_type'] == 'loss':
                r_cnv = np.ones_like(theta_cnv) * (1.0 + cnv['log2'] * 0.3)
                color = COLOR_LOSS
            else:
                continue

            ax.plot(theta_cnv, r_cnv, color=color, linewidth=2, alpha=0.7)

        # Configurações
        ax.set_ylim(0, 1.5)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.spines['polar'].set_visible(False)
        ax.grid(False)

        plt.title(
            f'Perfil Circular de CNVs - {self.sample_id}',
            fontsize=16, fontweight='bold', pad=20
        )

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"✓ Plot circular salvo: {output_file}")

        return fig

    def generate_aneuploidy_summary(self, cnvs):
        """
        Gera sumário textual de aneuploidias
        """
        chr_status = self.calculate_chromosome_aneuploidy(cnvs)

        summary = {
            'trisomies': [],
            'monosomies': [],
            'partial_gains': [],
            'partial_losses': [],
            'aneuploidy_score': 0
        }

        for chrom, status in chr_status.items():
            if status == 'trisomy':
                summary['trisomies'].append(chrom)
            elif status == 'monosomy':
                summary['monosomies'].append(chrom)
            elif status == 'partial_gain':
                summary['partial_gains'].append(chrom)
            elif status == 'partial_loss':
                summary['partial_losses'].append(chrom)

        # Score de aneuploidia (proporção do genoma alterado)
        altered_bases = 0
        total_bases = 0

        for chrom in self.chr_order:
            if chrom not in CHR_LENGTHS:
                continue

            chr_cnvs = cnvs[cnvs['chrom'] == chrom]
            chr_length = CHR_LENGTHS[chrom]

            altered = chr_cnvs[chr_cnvs['cnv_type'] != 'neutral']
            altered_bases += (altered['end'] - altered['start']).sum()
            total_bases += chr_length

        summary['aneuploidy_score'] = altered_bases / total_bases if total_bases > 0 else 0

        return summary


def main():
    parser = argparse.ArgumentParser(
        description='Gera ideograma de aneuploidias genoma completo'
    )
    parser.add_argument('--cnv', required=True, help='Arquivo CNVkit .cns')
    parser.add_argument('--validated', help='Arquivo CNVs validados (opcional)')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Arquivo PDF de saída')
    parser.add_argument('--output-circular', help='Plot circular (opcional)')
    parser.add_argument('--output-summary', help='Sumário textual (opcional)')

    args = parser.parse_args()

    # Inicializa plotter
    plotter = AneuploidyPlotter(args.sample)

    # Carrega dados
    cnvs = plotter.load_cnvs(args.cnv)
    validated = plotter.load_validated(args.validated) if args.validated else None

    # Gera ideograma
    plotter.plot_genome_ideogram(cnvs, validated, args.output)

    # Plot circular (opcional)
    if args.output_circular:
        plotter.plot_cnv_summary_circular(cnvs, args.output_circular)

    # Sumário textual
    summary = plotter.generate_aneuploidy_summary(cnvs)

    logger.info("\n" + "="*60)
    logger.info("SUMÁRIO DE ANEUPLOIDIAS")
    logger.info("="*60)
    logger.info(f"Score de aneuploidia: {summary['aneuploidy_score']:.2%}")

    if summary['trisomies']:
        logger.info(f"\nTrissomias: {', '.join(summary['trisomies'])}")
    if summary['monosomies']:
        logger.info(f"Monossomias: {', '.join(summary['monosomies'])}")
    if summary['partial_gains']:
        logger.info(f"Ganhos parciais: {', '.join(summary['partial_gains'])}")
    if summary['partial_losses']:
        logger.info(f"Perdas parciais: {', '.join(summary['partial_losses'])}")

    if args.output_summary:
        with open(args.output_summary, 'w') as f:
            f.write(f"Sample: {args.sample}\n")
            f.write(f"Aneuploidy Score: {summary['aneuploidy_score']:.4f}\n\n")
            f.write(f"Trisomies: {', '.join(summary['trisomies']) if summary['trisomies'] else 'None'}\n")
            f.write(f"Monosomies: {', '.join(summary['monosomies']) if summary['monosomies'] else 'None'}\n")
            f.write(f"Partial Gains: {', '.join(summary['partial_gains']) if summary['partial_gains'] else 'None'}\n")
            f.write(f"Partial Losses: {', '.join(summary['partial_losses']) if summary['partial_losses'] else 'None'}\n")

        logger.info(f"\n✓ Sumário salvo: {args.output_summary}")

    logger.info("\n✓ PLOTS GERADOS COM SUCESSO!")


if __name__ == '__main__':
    main()
