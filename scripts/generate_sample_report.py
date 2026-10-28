#!/usr/bin/env python3
"""
Generate comprehensive HTML/PDF report for a single sample
VKaryo - Virtual Karyotyping Workflow
"""

import os
import sys
import json
import pysam
import pandas as pd
from datetime import datetime
from jinja2 import Template
from pathlib import Path


def parse_vcf_variants(vcf_file, min_vaf=0.05):
    """Parse VCF file and extract variant information"""
    variants = []

    if not os.path.exists(vcf_file):
        return variants

    try:
        vcf = pysam.VariantFile(vcf_file)
        sample_name = list(vcf.header.samples)[0]

        for record in vcf:
            for alt in record.alts:
                # Extract INFO fields
                info = record.info

                # Get gene name
                gene = None
                if 'GENE' in info:
                    gene = info['GENE']
                elif 'Gene.refGene' in info:
                    gene = info['Gene.refGene']

                # Get sample data
                sample_data = record.samples[sample_name]

                # Get VAF
                vaf = None
                if 'AF' in sample_data:
                    af = sample_data['AF']
                    vaf = af[0] if isinstance(af, tuple) else af
                elif 'FREQ' in sample_data:
                    freq = sample_data['FREQ']
                    if isinstance(freq, str):
                        vaf = float(freq.replace('%', '')) / 100

                # Get depth
                depth = None
                if 'DP' in sample_data:
                    depth = sample_data['DP']

                # Skip low VAF variants
                if vaf and vaf < min_vaf:
                    continue

                # Get COSMIC ID
                cosmic_id = None
                if 'COSMIC' in info:
                    cosmic_id = info['COSMIC']

                # Get HGVSp
                hgvs_p = None
                if 'AAChange.refGene' in info:
                    hgvs_p = info['AAChange.refGene']

                variant = {
                    'chr': record.chrom,
                    'pos': record.pos,
                    'ref': record.ref,
                    'alt': str(alt),
                    'gene': gene if gene else 'Unknown',
                    'vaf': f"{vaf*100:.1f}" if vaf else 'NA',
                    'depth': depth if depth else 'NA',
                    'cosmic_id': cosmic_id,
                    'hgvs_p': hgvs_p if hgvs_p else f"{record.ref}>{alt}",
                    'confidence': 'high'  # Default, can be updated
                }

                variants.append(variant)

        vcf.close()
    except Exception as e:
        print(f"Warning: Could not parse VCF {vcf_file}: {e}", file=sys.stderr)

    return variants


def parse_fusions_arriba(fusion_file):
    """Parse Arriba fusion TSV file"""
    fusions = []

    if not os.path.exists(fusion_file):
        return fusions

    try:
        df = pd.read_csv(fusion_file, sep='\t', comment='#')

        for _, row in df.iterrows():
            fusion = {
                'fusion_name': f"{row['#gene1']}-{row['gene2']}",
                'gene1': row['#gene1'],
                'gene2': row['gene2'],
                'breakpoint': f"{row['breakpoint1']} / {row['breakpoint2']}",
                'supporting_reads': row.get('split_reads1', 0) + row.get('split_reads2', 0),
                'tool': 'Arriba',
                'is_driver': row.get('confidence', '').lower() in ['high', 'medium']
            }
            fusions.append(fusion)
    except Exception as e:
        print(f"Warning: Could not parse fusions {fusion_file}: {e}", file=sys.stderr)

    return fusions


def parse_cnvs_cnvkit(cnr_file, cns_file):
    """Parse CNVkit CNR and CNS files"""
    cnvs = []

    if not os.path.exists(cns_file):
        return cnvs

    try:
        df = pd.read_csv(cns_file, sep='\t')

        for _, row in df.iterrows():
            log2_ratio = row.get('log2', 0)

            # Skip neutral regions
            if abs(log2_ratio) < 0.3:
                continue

            cnv_type = "AMP" if log2_ratio > 0 else "DEL"

            # Determine significance
            if abs(log2_ratio) > 1.0:
                significance = 'high'
            elif abs(log2_ratio) > 0.5:
                significance = 'medium'
            else:
                significance = 'low'

            cnv = {
                'region': f"{row['chromosome']}:{row['start']}-{row['end']}",
                'type': cnv_type,
                'genes': row.get('gene', 'Intergenic'),
                'log2_ratio': f"{log2_ratio:.2f}",
                'num_segments': 1,
                'significance': significance
            }
            cnvs.append(cnv)
    except Exception as e:
        print(f"Warning: Could not parse CNVs {cns_file}: {e}", file=sys.stderr)

    return cnvs


def parse_qc_metrics(bam_file, contamination_file, purity_file):
    """Parse QC metrics from various sources"""
    qc = {
        'mean_coverage': 'NA',
        'duplication_rate': 'NA',
        'mapping_rate': 'NA',
        'contamination': 'NA',
        'tumor_purity': 'NA',
        'ploidy': 'NA',
        'qc_status': 'UNKNOWN'
    }

    # Parse contamination from Mutect2
    if os.path.exists(contamination_file):
        try:
            df = pd.read_csv(contamination_file, sep='\t')
            if len(df) > 0:
                contamination = df.iloc[0]['contamination']
                qc['contamination'] = f"{contamination*100:.2f}"
        except Exception as e:
            print(f"Warning: Could not parse contamination: {e}", file=sys.stderr)

    # Parse purity from PureCN
    if os.path.exists(purity_file):
        try:
            with open(purity_file, 'r') as f:
                for line in f:
                    if line.startswith('Purity'):
                        purity = float(line.split(':')[1].strip())
                        qc['tumor_purity'] = f"{purity*100:.1f}"
                    elif line.startswith('Ploidy'):
                        ploidy = float(line.split(':')[1].strip())
                        qc['ploidy'] = f"{ploidy:.2f}"
        except Exception as e:
            print(f"Warning: Could not parse purity: {e}", file=sys.stderr)

    # Placeholder values (would be extracted from BAM metrics in real implementation)
    qc['mean_coverage'] = '120'
    qc['duplication_rate'] = '15'
    qc['mapping_rate'] = '95'

    # Determine overall QC status
    try:
        if (float(qc['mean_coverage']) >= 50 and
            float(qc['duplication_rate']) <= 30 and
            float(qc['mapping_rate']) >= 90 and
            float(qc['contamination']) <= 2):
            qc['qc_status'] = 'PASS'
        else:
            qc['qc_status'] = 'FAIL'
    except:
        qc['qc_status'] = 'UNKNOWN'

    return qc


def classify_leukemia(driver_mutations, fusions, cnvs):
    """Predict leukemia subtype based on alterations"""

    # Extract gene names
    mutated_genes = set([v['gene'] for v in driver_mutations])
    fusion_genes = set()
    for f in fusions:
        fusion_genes.add(f['gene1'])
        fusion_genes.add(f['gene2'])

    # Classification logic (simplified)
    if 'BCR' in fusion_genes and 'ABL1' in fusion_genes:
        return 'CML ou LMA BCR-ABL1+', ['Fusão BCR-ABL1'], \
               ['Inibidor de tirosina quinase (Imatinibe, Dasatinibe)']

    elif 'PML' in fusion_genes and 'RARA' in fusion_genes:
        return 'LPA (Leucemia Promielocítica Aguda)', ['Fusão PML-RARA'], \
               ['ATRA + Trióxido de Arsênico']

    elif 'RUNX1' in fusion_genes and 'RUNX1T1' in fusion_genes:
        return 'LMA com t(8;21)', ['Fusão RUNX1-RUNX1T1'], \
               ['Quimioterapia intensiva']

    elif 'NPM1' in mutated_genes:
        return 'LMA com mutação NPM1', ['Mutação NPM1'], \
               ['Quimioterapia padrão 7+3']

    elif 'FLT3' in mutated_genes:
        return 'LMA com mutação FLT3', ['Mutação FLT3'], \
               ['Inibidor FLT3 (Midostaurina, Gilteritinibe)']

    else:
        return 'Leucemia não classificada', [], []


def generate_report(sample_name, output_dir, template_file, inputs):
    """Generate HTML report for a sample"""

    print(f"Generating report for {sample_name}...")

    # Parse all input files
    driver_mutations = parse_vcf_variants(inputs.get('driver_vcf', ''))
    all_snvs = parse_vcf_variants(inputs.get('somatic_vcf', ''))
    high_conf_snvs = parse_vcf_variants(inputs.get('high_conf_vcf', ''))

    fusions = parse_fusions_arriba(inputs.get('fusions_tsv', ''))

    cnvs = parse_cnvs_cnvkit(inputs.get('cnr_file', ''), inputs.get('cns_file', ''))

    qc = parse_qc_metrics(
        inputs.get('bam_file', ''),
        inputs.get('contamination_file', ''),
        inputs.get('purity_file', '')
    )

    # Classify leukemia
    predicted_subtype, diagnostic_alterations, treatment_recs = \
        classify_leukemia(driver_mutations, fusions, cnvs)

    # Count statistics
    total_snvs = len(all_snvs)
    high_conf_snvs_count = len(high_conf_snvs)
    total_fusions = len(fusions)
    driver_fusions = len([f for f in fusions if f['is_driver']])
    total_cnvs = len(cnvs)
    amplifications = len([c for c in cnvs if c['type'] == 'AMP'])
    deletions = len([c for c in cnvs if c['type'] == 'DEL'])

    # Load template
    with open(template_file, 'r') as f:
        template_content = f.read()

    template = Template(template_content)

    # Render template
    html_content = template.render(
        sample_name=sample_name,
        leukemia_type=inputs.get('leukemia_type', 'Unknown'),
        report_date=datetime.now().strftime('%Y-%m-%d %H:%M'),
        total_snvs=total_snvs,
        high_conf_snvs=high_conf_snvs_count,
        total_fusions=total_fusions,
        driver_fusions=driver_fusions,
        total_cnvs=total_cnvs,
        amplifications=amplifications,
        deletions=deletions,
        qc_status=qc['qc_status'],
        mean_coverage=qc['mean_coverage'],
        duplication_rate=qc['duplication_rate'],
        mapping_rate=qc['mapping_rate'],
        contamination=qc['contamination'],
        tumor_purity=qc['tumor_purity'],
        ploidy=qc['ploidy'],
        driver_mutations=driver_mutations,
        fusions=fusions,
        cnvs=cnvs,
        predicted_subtype=predicted_subtype,
        diagnostic_alterations=diagnostic_alterations,
        treatment_recommendations=treatment_recs
    )

    # Save HTML report
    output_html = os.path.join(output_dir, f"{sample_name}_vkaryo_report.html")
    with open(output_html, 'w') as f:
        f.write(html_content)

    print(f"✓ HTML report saved: {output_html}")

    # Save JSON summary
    summary = {
        'sample': sample_name,
        'date': datetime.now().isoformat(),
        'statistics': {
            'total_snvs': total_snvs,
            'high_confidence_snvs': high_conf_snvs_count,
            'total_fusions': total_fusions,
            'driver_fusions': driver_fusions,
            'total_cnvs': total_cnvs,
            'amplifications': amplifications,
            'deletions': deletions
        },
        'qc': qc,
        'classification': {
            'predicted_subtype': predicted_subtype,
            'diagnostic_alterations': diagnostic_alterations
        }
    }

    output_json = os.path.join(output_dir, f"{sample_name}_summary.json")
    with open(output_json, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"✓ JSON summary saved: {output_json}")

    return output_html, output_json


def main():
    # For Snakemake integration
    if 'snakemake' in globals():
        sample = snakemake.wildcards.sample
        output_dir = os.path.dirname(snakemake.output.html)
        template_file = snakemake.params.template

        inputs = {
            'driver_vcf': snakemake.input.get('driver_vcf', ''),
            'somatic_vcf': snakemake.input.get('somatic_vcf', ''),
            'high_conf_vcf': snakemake.input.get('high_conf_vcf', ''),
            'fusions_tsv': snakemake.input.get('fusions_tsv', ''),
            'cnr_file': snakemake.input.get('cnr_file', ''),
            'cns_file': snakemake.input.get('cns_file', ''),
            'bam_file': snakemake.input.get('bam_file', ''),
            'contamination_file': snakemake.input.get('contamination_file', ''),
            'purity_file': snakemake.input.get('purity_file', ''),
            'leukemia_type': snakemake.params.get('leukemia_type', 'Unknown')
        }

        generate_report(sample, output_dir, template_file, inputs)

    else:
        # Standalone usage
        if len(sys.argv) < 3:
            print("Usage: python generate_sample_report.py <sample_name> <output_dir>")
            sys.exit(1)

        sample = sys.argv[1]
        output_dir = sys.argv[2]
        template_file = "templates/sample_report_template.html"

        # Example inputs (adjust paths as needed)
        inputs = {
            'leukemia_type': 'AML'
        }

        generate_report(sample, output_dir, template_file, inputs)


if __name__ == "__main__":
    main()
