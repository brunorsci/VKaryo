#!/usr/bin/env python3
"""
Rastreia variantes de MRD e gera visualizações
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Parâmetros do Snakemake
mrd_report = snakemake.input.mrd_report
annotated = snakemake.input.annotated
output_tracking = snakemake.output.tracking
output_plot = snakemake.output.plot

def parse_mrd_report(report_path):
    """Extrai variantes do relatório MRD"""
    variants = []
    current_gene = None

    with open(report_path, 'r') as f:
        in_variants_section = False
        for line in f:
            line = line.strip()

            if "MRD VARIANTS DETECTED" in line:
                in_variants_section = True
                continue

            if "INTERPRETAÇÃO CLÍNICA" in line:
                break

            if in_variants_section and line and not line.startswith("="):
                if not line.startswith(" "):
                    # Nova seção de gene
                    current_gene = line.rstrip(":")
                elif ":" in line and ">" in line:
                    # Linha de variante: chr:pos ref>alt
                    var_info = line.strip().split()
                    position = var_info[0]
                    alleles = var_info[1]
                elif "VAF:" in line:
                    vaf_str = line.split("VAF:")[1].strip().rstrip("%")
                    vaf = float(vaf_str)

                    if current_gene:
                        variants.append({
                            'gene': current_gene,
                            'position': position,
                            'alleles': alleles,
                            'vaf': vaf
                        })

    return variants

def load_annovar(annovar_path):
    """Carrega anotações ANNOVAR"""
    try:
        df = pd.read_csv(annovar_path, sep='\t', low_memory=False)
        return df
    except:
        return None

# Parse variantes MRD
mrd_variants = parse_mrd_report(mrd_report)

# Carrega anotações ANNOVAR se disponível
annovar_df = load_annovar(annotated)

# Cria tabela de tracking
tracking_data = []
for var in mrd_variants:
    gene = var['gene']
    vaf = var['vaf']

    # Busca anotação adicional se disponível
    cosmic = "."
    clinvar = "."
    func = "."

    if annovar_df is not None:
        # Tenta encontrar variante no ANNOVAR
        matching = annovar_df[annovar_df['Gene.refGene'] == gene]
        if not matching.empty:
            row = matching.iloc[0]
            cosmic = row.get('cosmic70', '.')
            clinvar = row.get('CLINSIG', '.')
            func = row.get('Func.refGene', '.')

    tracking_data.append({
        'Gene': gene,
        'Position': var['position'],
        'Alleles': var['alleles'],
        'VAF (%)': vaf,
        'Function': func,
        'COSMIC': cosmic,
        'ClinVar': clinvar,
        'MRD_Level': classify_mrd_level(vaf)
    })

def classify_mrd_level(vaf):
    """Classifica nível de MRD"""
    if vaf < 0.01:
        return "Negative"
    elif vaf < 0.1:
        return "Low"
    elif vaf < 1.0:
        return "Moderate"
    else:
        return "High"

# Salva tracking table
tracking_df = pd.DataFrame(tracking_data)
tracking_df = tracking_df.sort_values('VAF (%)', ascending=False)
tracking_df.to_csv(output_tracking, sep='\t', index=False)

# Gera plot de MRD
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: VAF por gene
if not tracking_df.empty:
    genes = tracking_df['Gene'].values
    vafs = tracking_df['VAF (%)'].values

    # Cores por nível de MRD
    colors = []
    for vaf in vafs:
        if vaf >= 1.0:
            colors.append('#d62728')  # Vermelho - Alto
        elif vaf >= 0.1:
            colors.append('#ff7f0e')  # Laranja - Moderado
        elif vaf >= 0.01:
            colors.append('#ffdd57')  # Amarelo - Baixo
        else:
            colors.append('#2ca02c')  # Verde - Negativo

    bars = ax1.barh(genes, vafs, color=colors, edgecolor='black', linewidth=1.5)
    ax1.set_xlabel('VAF (%)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Gene', fontsize=12, fontweight='bold')
    ax1.set_title('MRD Variants - Allele Frequency', fontsize=14, fontweight='bold')
    ax1.axvline(x=0.01, color='green', linestyle='--', alpha=0.5, label='MRD threshold (0.01%)')
    ax1.axvline(x=0.1, color='orange', linestyle='--', alpha=0.5, label='Moderate (0.1%)')
    ax1.axvline(x=1.0, color='red', linestyle='--', alpha=0.5, label='High (1.0%)')
    ax1.legend(loc='lower right', fontsize=9)
    ax1.set_xscale('log')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Distribuição de MRD levels
    mrd_levels = tracking_df['MRD_Level'].value_counts()
    level_order = ['High', 'Moderate', 'Low', 'Negative']
    level_colors = {'High': '#d62728', 'Moderate': '#ff7f0e', 'Low': '#ffdd57', 'Negative': '#2ca02c'}

    levels_plot = []
    counts_plot = []
    colors_plot = []
    for level in level_order:
        if level in mrd_levels.index:
            levels_plot.append(level)
            counts_plot.append(mrd_levels[level])
            colors_plot.append(level_colors[level])

    ax2.bar(levels_plot, counts_plot, color=colors_plot, edgecolor='black', linewidth=1.5)
    ax2.set_xlabel('MRD Level', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Variants', fontsize=12, fontweight='bold')
    ax2.set_title('MRD Level Distribution', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

    # Adiciona valores nas barras
    for i, (level, count) in enumerate(zip(levels_plot, counts_plot)):
        ax2.text(i, count + 0.1, str(count), ha='center', va='bottom', fontweight='bold')

else:
    ax1.text(0.5, 0.5, 'No MRD variants detected\n(MRD Negative)',
             ha='center', va='center', fontsize=14, transform=ax1.transAxes,
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    ax2.text(0.5, 0.5, 'Excellent molecular response',
             ha='center', va='center', fontsize=14, transform=ax2.transAxes,
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

plt.tight_layout()
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ MRD tracking table saved: {output_tracking}")
print(f"✓ MRD plot saved: {output_plot}")
print(f"✓ Total MRD variants tracked: {len(tracking_data)}")
