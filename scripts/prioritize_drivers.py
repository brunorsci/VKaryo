#!/usr/bin/env python3
"""
Prioritize driver gene variants from merged Mutect2 + VarDict calls

NOTA: o nome do gene não pode ser extraído de campos INFO (GENE /
Gene.refGene / ANN) -- esses campos não existem nos VCFs brutos de
Mutect2/VarDict (essa anotação só aparece depois do ANNOVAR, em outro
arquivo). Usa sobreposição de coordenada genômica contra
driver_genes_leukemia.txt (mesma estratégia da consolidação de resultados).
"""

import pysam
import pandas as pd
import subprocess
import sys


def load_driver_genes(driver_file):
    """Load driver gene coordinates: [(gene, chrom, start, end), ...]"""
    genes = []
    with open(driver_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            gene, chrom, start, end = parts[0], parts[1], int(parts[2]), int(parts[3])
            genes.append((gene, chrom, start, end))
    return genes


def gene_at(chrom, pos, driver_coords):
    """Return gene name if (chrom, pos) falls inside a driver gene's body, else None."""
    for gene, g_chrom, g_start, g_end in driver_coords:
        if chrom == g_chrom and g_start <= pos <= g_end:
            return gene
    return None


def get_vaf(record, sample_name):
    """Extract VAF from VCF record (Mutect2 AF, VarDict AF/FREQ, or computed from AD)"""
    try:
        sample = record.samples[sample_name]
        if 'AF' in sample and sample['AF'] is not None:
            af = sample['AF']
            return float(af[0] if isinstance(af, (tuple, list)) else af)
        if 'FREQ' in sample and sample['FREQ'] is not None:
            freq = sample['FREQ']
            if isinstance(freq, str):
                return float(freq.replace('%', '')) / 100
            return float(freq)
        if 'AD' in sample and sample['AD'] is not None:
            ad = sample['AD']
            if isinstance(ad, (tuple, list)) and len(ad) == 2:
                ref_count, alt_count = ad
                total = (ref_count or 0) + (alt_count or 0)
                if total > 0:
                    return alt_count / total
    except Exception:
        pass
    return None


def prioritize_variants(merged_vcf, driver_coords, min_vaf, output_vcf, output_tsv, log):
    vcf_in = pysam.VariantFile(merged_vcf)
    sample_name = list(vcf_in.header.samples)[0]

    # pysam auto-bgzips when the output path ends in .gz, but never writes a
    # .tbi -- that has to be done explicitly afterwards.
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)

    variants_list = []
    driver_count = 0
    non_driver_count = 0

    for record in vcf_in:
        gene = gene_at(record.chrom, record.pos, driver_coords)
        vaf = get_vaf(record, sample_name)
        is_driver = gene is not None

        priority = 0
        if is_driver:
            priority += 10
            driver_count += 1
        else:
            non_driver_count += 1
        if vaf and vaf >= min_vaf:
            priority += int(vaf * 10)

        depth = None
        try:
            if 'DP' in record.samples[sample_name]:
                depth = record.samples[sample_name]['DP']
        except Exception:
            pass

        # Only the driver-gene subset goes into the prioritized VCF; the TSV
        # keeps every variant (ranked) for full traceability.
        if is_driver:
            vcf_out.write(record)

        variants_list.append({
            'Chromosome': record.chrom,
            'Position': record.pos,
            'Ref': record.ref,
            'Alt': ','.join(str(a) for a in record.alts) if record.alts else '',
            'Gene': gene if gene else 'Unknown',
            'VAF': f"{vaf:.4f}" if vaf is not None else 'NA',
            'Depth': depth if depth is not None else 'NA',
            'Is_Driver': 'Yes' if is_driver else 'No',
            'Priority': priority,
        })

    vcf_in.close()
    vcf_out.close()

    # Index the (bgzipped) driver-only VCF
    subprocess.run(['tabix', '-f', '-p', 'vcf', output_vcf], check=True)

    df = pd.DataFrame(variants_list)
    if len(df) > 0:
        df = df.sort_values('Priority', ascending=False)
    df.to_csv(output_tsv, sep='\t', index=False)

    with open(log, 'w') as f:
        f.write(f"Driver gene variants: {driver_count}\n")
        f.write(f"Non-driver variants: {non_driver_count}\n")
        f.write(f"Total variants: {len(variants_list)}\n")

    print("Prioritization complete:")
    print(f"  Driver genes: {driver_count}")
    print(f"  Non-driver: {non_driver_count}")
    print(f"  Total: {len(variants_list)}")


if __name__ == "__main__":
    merged_vcf = snakemake.input.merged
    driver_genes_file = snakemake.input.drivers
    output_vcf = snakemake.output.drivers_vcf
    output_tsv = snakemake.output.drivers_tsv
    min_vaf = snakemake.params.min_vaf
    log = snakemake.log[0]

    driver_coords = load_driver_genes(driver_genes_file)
    prioritize_variants(merged_vcf, driver_coords, min_vaf, output_vcf, output_tsv, log)
