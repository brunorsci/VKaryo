#!/usr/bin/env python3
"""
Filtra variantes germinativas usando databases públicos
Compensa falta de amostra normal pareada

Uso:
    python filter_germline_variants.py \
        --vcf sample.vcf.gz \
        --gnomad gnomad.genomes.v4.0.sites.vcf.gz \
        --dbsnp dbsnp_156.vcf.gz \
        --output-somatic sample_somatic.vcf.gz \
        --output-germline sample_germline.vcf.gz \
        --af-threshold 0.001
"""

import argparse
import pysam
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class GermlineFilter:
    """
    Filtra variantes germinativas sem amostra normal pareada
    """

    def __init__(self, af_threshold=0.001, gnomad_vcf=None, dbsnp_vcf=None):
        self.af_threshold = af_threshold
        self.gnomad_vcf = gnomad_vcf
        self.dbsnp_vcf = dbsnp_vcf

        self.stats = defaultdict(int)

    def classify_variant(self, variant, gnomad_db=None, dbsnp_db=None):
        """
        Classifica variante como somática ou germinativa
        """
        self.stats['total'] += 1

        # Critério 1: Frequência populacional (gnomAD)
        if gnomad_db:
            gnomad_af = self._get_gnomad_af(variant, gnomad_db)
            if gnomad_af and gnomad_af > self.af_threshold:
                self.stats['filtered_gnomad'] += 1
                return 'GERMLINE', f'gnomAD_AF={gnomad_af:.4f}'

        # Critério 2: Presente em dbSNP (comum)
        if dbsnp_db:
            in_dbsnp, rs_id = self._check_dbsnp(variant, dbsnp_db)
            if in_dbsnp:
                # Se dbSNP e não clinicamente relevante, provavelmente germinativa
                if not self._is_clinically_relevant(variant):
                    self.stats['filtered_dbsnp'] += 1
                    return 'GERMLINE', f'dbSNP={rs_id}'

        # Critério 3: VAF (Variant Allele Frequency)
        vaf = self._get_vaf(variant)
        if vaf:
            # VAF ~50% ou ~100% sugere germinativa (heterozigoto/homozigoto)
            # VAF baixo (<30%) mais provável somática
            if 0.4 <= vaf <= 0.6:
                self.stats['suspected_het_germline'] += 1
                return 'LIKELY_GERMLINE', f'VAF={vaf:.2f}_het'
            elif vaf > 0.9:
                self.stats['suspected_hom_germline'] += 1
                return 'LIKELY_GERMLINE', f'VAF={vaf:.2f}_hom'

        # Critério 4: Genes drivers conhecidos de leucemia
        if self._is_known_driver(variant):
            # Mesmo que em dbSNP, mutações driver são relevantes
            self.stats['somatic_driver'] += 1
            return 'SOMATIC', 'known_driver'

        # Default: sem evidência forte, considera somática (conservador)
        self.stats['somatic_no_evidence'] += 1
        return 'SOMATIC', 'no_population_evidence'

    def _get_gnomad_af(self, variant, gnomad_db):
        """
        Busca frequência alélica no gnomAD
        """
        try:
            # Busca variante exata
            for gnomad_var in gnomad_db.fetch(variant.contig, variant.pos-1, variant.pos+1):
                if (gnomad_var.pos == variant.pos and
                    gnomad_var.ref == variant.ref and
                    gnomad_var.alts and variant.alts and
                    gnomad_var.alts[0] == variant.alts[0]):

                    # Pega maior AF entre populações
                    af_fields = ['AF', 'AF_popmax', 'AF_nfe', 'AF_afr', 'AF_eas', 'AF_amr']
                    afs = []

                    for field in af_fields:
                        if field in gnomad_var.info:
                            val = gnomad_var.info[field]
                            if isinstance(val, (list, tuple)):
                                afs.extend(val)
                            else:
                                afs.append(val)

                    if afs:
                        return max(afs)

        except Exception as e:
            logger.debug(f"Error fetching gnomAD: {e}")

        return None

    def _check_dbsnp(self, variant, dbsnp_db):
        """
        Verifica se variante está no dbSNP
        """
        try:
            for dbsnp_var in dbsnp_db.fetch(variant.contig, variant.pos-1, variant.pos+1):
                if (dbsnp_var.pos == variant.pos and
                    dbsnp_var.ref == variant.ref and
                    dbsnp_var.alts and variant.alts and
                    dbsnp_var.alts[0] == variant.alts[0]):

                    rs_id = dbsnp_var.id or 'rs_unknown'
                    return True, rs_id

        except Exception as e:
            logger.debug(f"Error fetching dbSNP: {e}")

        return False, None

    def _get_vaf(self, variant):
        """
        Extrai VAF (Variant Allele Frequency) do VCF
        """
        try:
            # Tenta pegar do campo AF (se existe)
            if 'AF' in variant.info:
                af = variant.info['AF']
                if isinstance(af, (list, tuple)):
                    return af[0]
                return af

            # Se não, calcula de AD (Allele Depth)
            sample = list(variant.samples.values())[0]

            if 'AD' in sample:
                ad = sample['AD']
                if isinstance(ad, (list, tuple)) and len(ad) >= 2:
                    ref_depth = ad[0]
                    alt_depth = ad[1]
                    total = ref_depth + alt_depth

                    if total > 0:
                        return alt_depth / total

        except Exception as e:
            logger.debug(f"Error getting VAF: {e}")

        return None

    def _is_known_driver(self, variant):
        """
        Verifica se variante está em gene driver conhecido de leucemia
        """
        # Lista de genes driver em leucemias
        leukemia_drivers = {
            'FLT3', 'NPM1', 'CEBPA', 'IDH1', 'IDH2', 'TET2', 'DNMT3A',
            'TP53', 'ASXL1', 'RUNX1', 'KIT', 'JAK2', 'CALR', 'MPL',
            'SF3B1', 'SRSF2', 'U2AF1', 'ZRSR2',  # Splicing
            'EZH2', 'SUZ12', 'SETD2',  # Chromatin
            'NRAS', 'KRAS', 'PTPN11',  # RAS pathway
            'PHF6', 'BCOR', 'BCORL1', 'STAG2'  # Others
        }

        # Pega nome do gene (pode estar em diferentes campos)
        gene = None

        # INFO field pode ter gene annotation
        if 'GENEINFO' in variant.info:
            gene = variant.info['GENEINFO'].split(':')[0]
        elif 'Gene.refGene' in variant.info:
            gene = variant.info['Gene.refGene']

        # Ou no ID (algumas annotations)
        elif variant.id and '|' in variant.id:
            parts = variant.id.split('|')
            if len(parts) > 1:
                gene = parts[1]

        if gene and gene in leukemia_drivers:
            return True

        return False

    def _is_clinically_relevant(self, variant):
        """
        Verifica se variante tem relevância clínica conhecida (ClinVar)
        """
        # Verifica CLNSIG (ClinVar significance)
        if 'CLNSIG' in variant.info:
            clnsig = variant.info['CLNSIG']
            if isinstance(clnsig, str):
                clnsig = [clnsig]

            # Pathogenic ou Likely_pathogenic
            for sig in clnsig:
                if 'pathogenic' in sig.lower() and 'benign' not in sig.lower():
                    return True

        return False

    def filter_vcf(self, input_vcf, output_somatic, output_germline, output_report):
        """
        Filtra VCF e separa somáticas de germinativas
        """
        logger.info(f"Filtering VCF: {input_vcf}")

        # Abre databases
        gnomad_db = pysam.VariantFile(self.gnomad_vcf) if self.gnomad_vcf else None
        dbsnp_db = pysam.VariantFile(self.dbsnp_vcf) if self.dbsnp_vcf else None

        # Abre input/output VCFs
        vcf_in = pysam.VariantFile(input_vcf)

        # Adiciona header info
        vcf_in.header.add_line('##INFO=<ID=CLASSIFICATION,Number=1,Type=String,Description="Variant classification: SOMATIC, GERMLINE, LIKELY_GERMLINE">')
        vcf_in.header.add_line('##INFO=<ID=FILTER_REASON,Number=1,Type=String,Description="Reason for classification">')

        vcf_somatic = pysam.VariantFile(output_somatic, 'w', header=vcf_in.header)
        vcf_germline = pysam.VariantFile(output_germline, 'w', header=vcf_in.header)

        # Processa cada variante
        for variant in vcf_in:
            classification, reason = self.classify_variant(variant, gnomad_db, dbsnp_db)

            # Adiciona info ao variant
            variant.info['CLASSIFICATION'] = classification
            variant.info['FILTER_REASON'] = reason

            # Escreve no arquivo apropriado
            if classification == 'SOMATIC':
                vcf_somatic.write(variant)
            else:  # GERMLINE ou LIKELY_GERMLINE
                vcf_germline.write(variant)

        # Fecha arquivos
        vcf_in.close()
        vcf_somatic.close()
        vcf_germline.close()

        if gnomad_db:
            gnomad_db.close()
        if dbsnp_db:
            dbsnp_db.close()

        # Relatório
        self._write_report(output_report)

        logger.info(f"✓ Filtering complete")
        logger.info(f"  Somatic variants: {output_somatic}")
        logger.info(f"  Germline variants: {output_germline}")
        logger.info(f"  Report: {output_report}")

    def _write_report(self, report_file):
        """
        Escreve relatório de filtragem
        """
        with open(report_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("GERMLINE FILTERING REPORT\n")
            f.write("=" * 60 + "\n\n")

            f.write(f"Total variants: {self.stats['total']}\n\n")

            f.write("Germline filtering:\n")
            f.write(f"  Filtered by gnomAD (AF > {self.af_threshold}): {self.stats['filtered_gnomad']}\n")
            f.write(f"  Filtered by dbSNP: {self.stats['filtered_dbsnp']}\n")
            f.write(f"  Suspected heterozygous germline (VAF ~50%): {self.stats['suspected_het_germline']}\n")
            f.write(f"  Suspected homozygous germline (VAF >90%): {self.stats['suspected_hom_germline']}\n\n")

            total_germline = (
                self.stats['filtered_gnomad'] +
                self.stats['filtered_dbsnp'] +
                self.stats['suspected_het_germline'] +
                self.stats['suspected_hom_germline']
            )

            f.write(f"Total germline/likely germline: {total_germline}\n\n")

            f.write("Somatic variants:\n")
            f.write(f"  Known driver genes: {self.stats['somatic_driver']}\n")
            f.write(f"  No population evidence: {self.stats['somatic_no_evidence']}\n\n")

            total_somatic = self.stats['somatic_driver'] + self.stats['somatic_no_evidence']

            f.write(f"Total somatic: {total_somatic}\n\n")

            f.write("=" * 60 + "\n")

        logger.info(f"Report statistics:")
        logger.info(f"  Total: {self.stats['total']}")
        logger.info(f"  Germline: {total_germline} ({total_germline/self.stats['total']*100:.1f}%)")
        logger.info(f"  Somatic: {total_somatic} ({total_somatic/self.stats['total']*100:.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description='Filtra variantes germinativas sem amostra normal'
    )
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--gnomad', help='gnomAD VCF database')
    parser.add_argument('--dbsnp', help='dbSNP VCF database')
    parser.add_argument('--output-somatic', required=True, help='Output somatic VCF')
    parser.add_argument('--output-germline', required=True, help='Output germline VCF')
    parser.add_argument('--output-report', help='Output report file')
    parser.add_argument('--af-threshold', type=float, default=0.001,
                        help='Population AF threshold for germline (default: 0.001)')

    args = parser.parse_args()

    if not args.output_report:
        args.output_report = args.output_somatic.replace('.vcf', '_report.txt')

    filter_tool = GermlineFilter(
        af_threshold=args.af_threshold,
        gnomad_vcf=args.gnomad,
        dbsnp_vcf=args.dbsnp
    )

    filter_tool.filter_vcf(
        args.vcf,
        args.output_somatic,
        args.output_germline,
        args.output_report
    )


if __name__ == '__main__':
    main()
