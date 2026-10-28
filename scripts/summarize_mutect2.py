#!/usr/bin/env python3
"""
Sumariza variantes somáticas do Mutect2
Foco em driver mutations para leucemias
"""

import pandas as pd
import gzip
import os
import argparse
from collections import defaultdict

# Parâmetros via argparse ou Snakemake
try:
    annotated = snakemake.input.annotated
    vcf_file = snakemake.input.vcf
    output_summary = snakemake.output.summary
    output_drivers = snakemake.output.drivers
except NameError:
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True)
    parser.add_argument('--annotated', required=True)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--drivers', required=True)
    args = parser.parse_args()
    vcf_file = args.vcf
    annotated = args.annotated
    output_summary = args.summary
    output_drivers = args.drivers

# Genes drivers em leucemias
LEUKEMIA_DRIVERS = {
    # AML
    'FLT3', 'NPM1', 'DNMT3A', 'TET2', 'IDH1', 'IDH2', 'TP53', 'RUNX1',
    'CEBPA', 'KIT', 'NRAS', 'KRAS', 'PTPN11', 'WT1', 'PHF6', 'ASXL1',
    'EZH2', 'CBL', 'GATA2', 'STAG2', 'SMC1A', 'SMC3', 'RAD21',
    # ALL
    'NOTCH1', 'FBXW7', 'PTEN', 'CDKN2A', 'PAX5', 'IKZF1', 'JAK2', 'JAK1',
    'IL7R', 'CRLF2', 'RAS', 'NT5C2',
    # CML/MPN
    'ABL1', 'BCR', 'JAK2', 'CALR', 'MPL', 'CSF3R',
    # CLL
    'SF3B1', 'ATM', 'BIRC3', 'NOTCH1', 'MYD88', 'TP53'
}

# Categorias de impacto
HIGH_IMPACT = {'frameshift_deletion', 'frameshift_insertion', 'stopgain',
               'stoploss', 'splicing', 'nonsynonymous_SNV'}
MODERATE_IMPACT = {'nonframeshift_deletion', 'nonframeshift_insertion'}

def parse_vcf_stats(vcf_path):
    """Extrai estatísticas básicas do VCF"""
    total = 0
    by_type = defaultdict(int)

    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            total += 1
            fields = line.strip().split('\t')
            ref = fields[3]
            alt = fields[4]

            # Classifica tipo
            if len(ref) == 1 and len(alt) == 1:
                by_type['SNV'] += 1
            elif len(ref) < len(alt):
                by_type['INS'] += 1
            elif len(ref) > len(alt):
                by_type['DEL'] += 1
            else:
                by_type['MNV'] += 1

    return total, by_type

def classify_variant_impact(func, exonic_func):
    """Classifica impacto funcional"""
    if func in ['exonic', 'exonic;splicing']:
        if exonic_func in HIGH_IMPACT:
            return 'HIGH'
        elif exonic_func in MODERATE_IMPACT:
            return 'MODERATE'
        else:
            return 'LOW'
    elif 'splicing' in func:
        return 'HIGH'
    elif func in ['UTR3', 'UTR5']:
        return 'LOW'
    elif func == 'intronic':
        return 'MODIFIER'
    else:
        return 'MODIFIER'

# Carrega anotações ANNOVAR (se disponível)
df = pd.DataFrame()
if annotated and os.path.exists(annotated):
    try:
        df = pd.read_csv(annotated, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Warning: Could not load ANNOVAR annotations: {e}")
        df = pd.DataFrame()
else:
    print("Note: ANNOVAR annotations not available, skipping driver mutation analysis")

# Estatísticas do VCF
total_variants, variant_types = parse_vcf_stats(vcf_file)

# Classifica origem da variante (somática vs germinativa)
def classify_origin(gnomad_af):
    """Classifica origem baseado na frequência alélica no gnomAD"""
    try:
        af = float(gnomad_af)
        if af >= 0.01:
            return 'Germline'
        elif af >= 0.001:
            return 'Likely_Germline'
        else:
            return 'Somatic'
    except (ValueError, TypeError):
        return 'Somatic'  # Se não há dado de gnomAD, assume somática

# Identifica driver mutations
drivers = []
if not df.empty and 'Gene.refGene' in df.columns:
    for idx, row in df.iterrows():
        gene = row.get('Gene.refGene', '.')

        # Checa se é driver
        if gene in LEUKEMIA_DRIVERS:
            func = row.get('Func.refGene', '.')
            exonic = row.get('ExonicFunc.refGene', '.')
            cosmic = row.get('cosmic101', row.get('cosmic70', '.'))
            clinvar = row.get('CLNSIG', row.get('clinvar_20240917', '.'))
            gnomad_af = row.get('gnomad41_genome_AF', '.')

            impact = classify_variant_impact(func, exonic)
            origin = classify_origin(gnomad_af)

            # Apenas HIGH/MODERATE impact
            if impact in ['HIGH', 'MODERATE']:
                drivers.append({
                    'Gene': gene,
                    'Chr': row.get('Chr', '.'),
                    'Start': row.get('Start', '.'),
                    'Ref': row.get('Ref', '.'),
                    'Alt': row.get('Alt', '.'),
                    'Function': func,
                    'ExonicFunc': exonic,
                    'Impact': impact,
                    'Origin': origin,
                    'gnomAD_AF': gnomad_af,
                    'COSMIC': cosmic,
                    'ClinVar': clinvar,
                    'AAChange': row.get('AAChange.refGene', '.')
                })

# Salva drivers
drivers_df = pd.DataFrame(drivers)
if not drivers_df.empty:
    drivers_df = drivers_df.sort_values(['Origin', 'Impact', 'Gene'])
drivers_df.to_csv(output_drivers, sep='\t', index=False)

# Separa somáticas e germinativas
somatic_drivers = drivers_df[drivers_df['Origin'] == 'Somatic'] if not drivers_df.empty else pd.DataFrame()
germline_drivers = drivers_df[drivers_df['Origin'].isin(['Germline', 'Likely_Germline'])] if not drivers_df.empty else pd.DataFrame()

# Gera relatório sumário
with open(output_summary, 'w') as report:
    report.write("=" * 60 + "\n")
    report.write("MUTECT2 VARIANTS SUMMARY\n")
    report.write("=" * 60 + "\n\n")

    report.write(f"Total PASS variants: {total_variants}\n\n")

    report.write("Variant types:\n")
    for vtype, count in sorted(variant_types.items()):
        report.write(f"  {vtype}: {count} ({count/total_variants*100:.1f}%)\n")
    report.write("\n")

    report.write("=" * 60 + "\n")
    report.write("SOMATIC DRIVER MUTATIONS\n")
    report.write("=" * 60 + "\n\n")

    if not somatic_drivers.empty:
        report.write(f"Total somatic drivers: {len(somatic_drivers)}\n")
        report.write(f"  HIGH impact: {sum(somatic_drivers['Impact'] == 'HIGH')}\n")
        report.write(f"  MODERATE impact: {sum(somatic_drivers['Impact'] == 'MODERATE')}\n\n")

        by_gene = somatic_drivers.groupby('Gene').size().sort_values(ascending=False)

        for gene, count in by_gene.items():
            report.write(f"\n{gene} ({count} variant(s)):\n")
            gene_vars = somatic_drivers[somatic_drivers['Gene'] == gene]

            for _, var in gene_vars.iterrows():
                report.write(f"  {var['Chr']}:{var['Start']} {var['Ref']}>{var['Alt']}\n")
                report.write(f"    Function: {var['ExonicFunc']}\n")
                report.write(f"    Impact: {var['Impact']}\n")
                report.write(f"    gnomAD AF: {var['gnomAD_AF']}\n")
                if var['AAChange'] != '.':
                    report.write(f"    AA Change: {var['AAChange']}\n")
                if str(var['COSMIC']) != '.':
                    report.write(f"    COSMIC: {var['COSMIC']}\n")
                if str(var['ClinVar']) != '.':
                    report.write(f"    ClinVar: {var['ClinVar']}\n")
    else:
        report.write("No somatic driver mutations detected.\n")

    report.write("\n" + "=" * 60 + "\n")
    report.write("GERMLINE VARIANTS IN DRIVER GENES\n")
    report.write("=" * 60 + "\n\n")

    if not germline_drivers.empty:
        report.write(f"Total germline variants in driver genes: {len(germline_drivers)}\n\n")
        by_gene = germline_drivers.groupby('Gene').size().sort_values(ascending=False)

        for gene, count in by_gene.items():
            report.write(f"\n{gene} ({count} variant(s)):\n")
            gene_vars = germline_drivers[germline_drivers['Gene'] == gene]

            for _, var in gene_vars.iterrows():
                report.write(f"  {var['Chr']}:{var['Start']} {var['Ref']}>{var['Alt']}\n")
                report.write(f"    Function: {var['ExonicFunc']}\n")
                report.write(f"    gnomAD AF: {var['gnomAD_AF']}\n")
                if str(var['ClinVar']) != '.':
                    report.write(f"    ClinVar: {var['ClinVar']}\n")
    else:
        report.write("No germline variants in driver genes.\n")

    if not somatic_drivers.empty:
        report.write("\n" + "=" * 60 + "\n")
        report.write("CLINICAL INTERPRETATION\n")
        report.write("=" * 60 + "\n\n")

        by_gene = somatic_drivers.groupby('Gene').size().sort_values(ascending=False)

        if 'FLT3' in by_gene:
            report.write("• FLT3: Targetable with FLT3 inhibitors (midostaurin, gilteritinib)\n")
        if 'IDH1' in by_gene or 'IDH2' in by_gene:
            report.write("• IDH1/2: Targetable with IDH inhibitors (ivosidenib, enasidenib)\n")
        if 'NPM1' in by_gene:
            report.write("• NPM1: Favorable prognosis marker in AML\n")
        if 'TP53' in by_gene:
            report.write("• TP53: Poor prognosis, consider alternative therapies\n")
        if 'JAK2' in by_gene:
            report.write("• JAK2: Consider JAK inhibitors (ruxolitinib)\n")
        if 'BCR' in by_gene or 'ABL1' in by_gene:
            report.write("• BCR-ABL1: Use tyrosine kinase inhibitors (imatinib, etc.)\n")

        report.write("\nRecommendations:\n")
        report.write("- Correlate with clinical presentation and cytogenetics\n")
        report.write("- Consider molecular profiling panel for therapy selection\n")
        report.write("- Assess for targeted therapy eligibility\n")
        report.write("- Use for MRD monitoring post-treatment\n")

    else:
        report.write("No driver mutations detected in known leukemia genes.\n\n")
        report.write("Note:\n")
        report.write("- May indicate favorable cytogenetics-only disease\n")
        report.write("- Review complete variant list for novel alterations\n")
        report.write("- Consider RNA-seq for fusion detection\n")

print(f"Summary saved to {output_summary}")
print(f"Driver mutations saved to {output_drivers}")
print(f"Total variants: {total_variants}")
print(f"Total drivers: {len(drivers_df)} (Somatic: {len(somatic_drivers)}, Germline: {len(germline_drivers)})")
