#!/usr/bin/env python3
"""
Filtra variantes VarDict para MRD clínico em leucemias
Foco em genes recorrentes e baixa VAF
"""

import sys
import gzip
import subprocess
import os
from pathlib import Path

# Parâmetros do Snakemake
input_vcf = snakemake.input.vcf
output_filtered = snakemake.output.filtered
output_report = snakemake.output.mrd_report

# Parâmetros
mrd_genes = set(snakemake.params.mrd_genes.split(","))
min_vaf_mrd = float(snakemake.params.min_vaf_mrd)
max_vaf_germline = float(snakemake.params.max_vaf_germline)

def parse_info(info_str):
    """Parse VCF INFO field"""
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict

def get_vaf(info_dict):
    """Extrai VAF do INFO field VarDict"""
    # VarDict usa AF para allele frequency
    if "AF" in info_dict:
        return float(info_dict["AF"])
    return 0.0

def get_gene(info_dict):
    """Extrai nome do gene"""
    # VarDict pode usar GENE tag
    if "GENE" in info_dict:
        return info_dict["GENE"]
    return None

def is_actionable_variant(chrom, pos, ref, alt, info_dict, gene):
    """
    Define se variante é acionável para MRD
    Critérios:
    - Gene recorrente em leucemia
    - VAF entre min_mrd e max_germline
    - Qualidade adequada
    """
    vaf = get_vaf(info_dict)

    # Filtro VAF: evita germline e ruído
    if vaf < min_vaf_mrd or vaf > max_vaf_germline:
        return False, vaf

    # Filtro por gene
    if gene and gene in mrd_genes:
        return True, vaf

    return False, vaf


# Processa VCF
mrd_variants = []
total_variants = 0
filtered_count = 0

opener = gzip.open if input_vcf.endswith('.gz') else open

# Cria VCF temporário (não comprimido)
temp_vcf = str(output_filtered).replace('.vcf.gz', '.tmp.vcf')

with opener(input_vcf, 'rt') as vcf_in:
    # Escreve VCF filtrado temporário (não comprimido)
    with open(temp_vcf, 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith('#'):
                vcf_out.write(line)
                continue

            total_variants += 1
            fields = line.strip().split('\t')
            chrom, pos, id_col, ref, alt = fields[0:5]
            qual, filter_col, info = fields[5:8]

            # Parse INFO
            info_dict = parse_info(info)
            gene = get_gene(info_dict)

            # Filtra para MRD
            is_mrd, vaf = is_actionable_variant(chrom, pos, ref, alt, info_dict, gene)

            if is_mrd:
                vcf_out.write(line)
                filtered_count += 1
                mrd_variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gene': gene,
                    'vaf': vaf,
                    'qual': qual
                })

# Comprimir com bgzip (BGZF format, compatível com bcftools/tabix)
subprocess.run(['bgzip', '-f', temp_vcf], check=True)
# bgzip adiciona .gz automaticamente, então renomear se necessário
import os
if os.path.exists(temp_vcf + '.gz'):
    os.rename(temp_vcf + '.gz', output_filtered)

# Gera relatório MRD
with open(output_report, 'w') as report:
    report.write("=" * 60 + "\n")
    report.write("MINIMAL RESIDUAL DISEASE (MRD) - VarDict\n")
    report.write("=" * 60 + "\n\n")

    report.write(f"Total variants called: {total_variants}\n")
    report.write(f"MRD-relevant variants: {filtered_count}\n")
    report.write(f"Filter rate: {(1 - filtered_count/total_variants)*100:.1f}%\n\n")

    if mrd_variants:
        report.write("=" * 60 + "\n")
        report.write("MRD VARIANTS DETECTED\n")
        report.write("=" * 60 + "\n\n")

        # Agrupa por gene
        by_gene = {}
        for var in mrd_variants:
            gene = var['gene']
            if gene not in by_gene:
                by_gene[gene] = []
            by_gene[gene].append(var)

        for gene, variants in sorted(by_gene.items()):
            report.write(f"\n{gene}:\n")
            for var in variants:
                report.write(f"  {var['chrom']}:{var['pos']} {var['ref']}>{var['alt']}\n")
                report.write(f"    VAF: {var['vaf']*100:.3f}%\n")
                report.write(f"    Quality: {var['qual']}\n")

        report.write("\n" + "=" * 60 + "\n")
        report.write("INTERPRETAÇÃO CLÍNICA\n")
        report.write("=" * 60 + "\n")
        report.write("\nMRD Levels:\n")
        report.write("  < 0.01% (10^-4): MRD negativo molecular\n")
        report.write("  0.01-0.1%: MRD baixo\n")
        report.write("  0.1-1%: MRD moderado\n")
        report.write("  > 1%: MRD alto (considerar doença residual)\n\n")

        # Classifica MRD level
        max_vaf = max(v['vaf'] for v in mrd_variants) * 100
        if max_vaf < 0.01:
            mrd_status = "NEGATIVO"
        elif max_vaf < 0.1:
            mrd_status = "BAIXO"
        elif max_vaf < 1.0:
            mrd_status = "MODERADO"
        else:
            mrd_status = "ALTO"

        report.write(f"MRD Status: {mrd_status}\n")
        report.write(f"Highest VAF detected: {max_vaf:.3f}%\n\n")

        report.write("Recomendação:\n")
        if mrd_status == "NEGATIVO":
            report.write("- Continuar monitoramento\n")
            report.write("- Considerar espaçamento de análises\n")
        elif mrd_status == "BAIXO":
            report.write("- Monitoramento frequente (30-60 dias)\n")
            report.write("- Avaliar resposta ao tratamento\n")
        elif mrd_status == "MODERADO":
            report.write("- Monitoramento intensivo (15-30 dias)\n")
            report.write("- Considerar intensificação terapêutica\n")
        else:
            report.write("- URGENTE: Avaliar falha terapêutica\n")
            report.write("- Considerar mudança de estratégia\n")
            report.write("- Descartar recaída molecular\n")
    else:
        report.write("\nNENHUMA variante MRD-relevante detectada.\n")
        report.write("Status: MRD NEGATIVO\n\n")
        report.write("Interpretação:\n")
        report.write("- Excelente resposta molecular\n")
        report.write("- Continuar monitoramento conforme protocolo\n")

print(f"✓ Filtered {filtered_count}/{total_variants} MRD-relevant variants")
print(f"✓ Report saved to {output_report}")
