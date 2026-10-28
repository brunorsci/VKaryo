#!/usr/bin/env python3
"""
Create consensus variant calls from Mutect2 + VarDict
Classify variants by confidence based on caller agreement and driver gene status
"""

import pysam
import pandas as pd
from collections import defaultdict


def load_driver_genes(driver_file):
    """Load list of driver genes"""
    with open(driver_file, 'r') as f:
        drivers = set(line.strip() for line in f if line.strip())
    return drivers


def load_vcf_variants(vcf_file):
    """Load all variants from VCF into a set of (chrom, pos, ref, alt) tuples"""
    variants = set()
    vcf = pysam.VariantFile(vcf_file)

    for record in vcf:
        for alt in record.alts:
            variant_key = (record.chrom, record.pos, record.ref, str(alt))
            variants.add(variant_key)

    vcf.close()
    return variants


def extract_gene_from_record(record):
    """Extract gene name from VCF record"""
    info = record.info

    if 'GENE' in info:
        return info['GENE']
    elif 'Gene.refGene' in info:
        return info['Gene.refGene']
    elif 'ANN' in info:
        ann = info['ANN'][0] if isinstance(info['ANN'], list) else info['ANN']
        fields = ann.split('|')
        if len(fields) > 3:
            return fields[3]

    return None


def consensus_calling(merged_vcf, mutect2_vcf, vardict_vcf, driver_genes_file,
                     high_conf_vcf, medium_conf_vcf, summary_file, log_file):
    """
    Create consensus calls:
    - HIGH: Detected by both callers
    - MEDIUM: Detected by one caller in driver gene
    - LOW: Detected by one caller in non-driver gene
    """

    # Load driver genes
    driver_genes = load_driver_genes(driver_genes_file)

    # Load variants from each caller
    print("Loading Mutect2 variants...")
    mutect2_variants = load_vcf_variants(mutect2_vcf)

    print("Loading VarDict variants...")
    vardict_variants = load_vcf_variants(vardict_vcf)

    # Open merged VCF
    merged = pysam.VariantFile(merged_vcf)
    header = merged.header

    # Create output VCFs
    high_conf = pysam.VariantFile(high_conf_vcf, 'w', header=header)
    medium_conf = pysam.VariantFile(medium_conf_vcf, 'w', header=header)

    # Statistics
    stats = {
        'high_confidence': 0,
        'medium_confidence': 0,
        'low_confidence': 0,
        'total': 0
    }

    variant_details = []

    # Process each variant in merged VCF
    for record in merged:
        stats['total'] += 1

        # Create variant key
        for alt in record.alts:
            variant_key = (record.chrom, record.pos, record.ref, str(alt))

            # Check which callers detected this variant
            in_mutect2 = variant_key in mutect2_variants
            in_vardict = variant_key in vardict_variants

            # Extract gene
            gene = extract_gene_from_record(record)
            is_driver = gene in driver_genes if gene else False

            # Classify confidence
            if in_mutect2 and in_vardict:
                # HIGH: Both callers agree
                confidence = 'HIGH'
                high_conf.write(record)
                stats['high_confidence'] += 1
            elif is_driver:
                # MEDIUM: One caller, but in driver gene
                confidence = 'MEDIUM'
                medium_conf.write(record)
                stats['medium_confidence'] += 1
            else:
                # LOW: One caller, non-driver gene
                confidence = 'LOW'
                stats['low_confidence'] += 1

            # Collect details
            variant_details.append({
                'Chromosome': record.chrom,
                'Position': record.pos,
                'Ref': record.ref,
                'Alt': str(alt),
                'Gene': gene if gene else 'Unknown',
                'Is_Driver': 'Yes' if is_driver else 'No',
                'In_Mutect2': 'Yes' if in_mutect2 else 'No',
                'In_VarDict': 'Yes' if in_vardict else 'No',
                'Confidence': confidence
            })

    # Close VCFs
    merged.close()
    high_conf.close()
    medium_conf.close()

    # Write summary
    with open(summary_file, 'w') as f:
        f.write("Consensus Variant Calling Summary\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"Total variants in merged callset: {stats['total']}\n\n")

        f.write("Confidence Breakdown:\n")
        f.write(f"  HIGH confidence (both callers): {stats['high_confidence']} "
                f"({stats['high_confidence']/stats['total']*100:.1f}%)\n")
        f.write(f"  MEDIUM confidence (driver gene): {stats['medium_confidence']} "
                f"({stats['medium_confidence']/stats['total']*100:.1f}%)\n")
        f.write(f"  LOW confidence (single caller): {stats['low_confidence']} "
                f"({stats['low_confidence']/stats['total']*100:.1f}%)\n\n")

        f.write("Caller Statistics:\n")
        f.write(f"  Mutect2 unique variants: {len(mutect2_variants)}\n")
        f.write(f"  VarDict unique variants: {len(vardict_variants)}\n")
        f.write(f"  Overlapping variants: {stats['high_confidence']}\n\n")

    # Write detailed TSV
    df = pd.DataFrame(variant_details)
    tsv_file = summary_file.replace('.txt', '_detailed.tsv')
    df.to_csv(tsv_file, sep='\t', index=False)

    # Log
    with open(log_file, 'w') as f:
        f.write(f"Consensus calling completed\n")
        f.write(f"HIGH: {stats['high_confidence']}\n")
        f.write(f"MEDIUM: {stats['medium_confidence']}\n")
        f.write(f"LOW: {stats['low_confidence']}\n")

    print(f"\nConsensus calling complete:")
    print(f"  HIGH confidence: {stats['high_confidence']}")
    print(f"  MEDIUM confidence: {stats['medium_confidence']}")
    print(f"  LOW confidence: {stats['low_confidence']}")


if __name__ == "__main__":
    # Snakemake inputs
    merged_vcf = snakemake.input.merged
    mutect2_vcf = snakemake.input.mutect2
    vardict_vcf = snakemake.input.vardict
    driver_genes_file = snakemake.input.drivers

    # Snakemake outputs
    high_conf_vcf = snakemake.output.high_conf
    medium_conf_vcf = snakemake.output.medium_conf
    summary_file = snakemake.output.summary

    # Log
    log_file = snakemake.log[0]

    # Run consensus calling
    consensus_calling(merged_vcf, mutect2_vcf, vardict_vcf, driver_genes_file,
                     high_conf_vcf, medium_conf_vcf, summary_file, log_file)
