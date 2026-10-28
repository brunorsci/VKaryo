#!/usr/bin/env python3
"""
Script para criar Panel of Normals (PoN) para CNVkit
Baixa e processa amostras normais públicas para compensar falta de controles pareados

Uso:
    python create_pon_cnvkit.py --output reference_pon.cnn --targets targets.bed
"""

import subprocess
import os
import sys
import argparse
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PanelOfNormalsBuilder:
    """
    Constrói Panel of Normals para CNVkit usando amostras públicas
    """

    def __init__(self, output_dir="pon_data", ref_genome=None, targets_bed=None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.ref_genome = ref_genome
        self.targets_bed = targets_bed

        # Amostras normais públicas recomendadas
        self.normal_samples = {
            # 1000 Genomes - Trios CEPH (alta qualidade)
            '1000G': [
                'NA12878',  # Filha do trio CEPH/Utah
                'NA12891',  # Pai
                'NA12892',  # Mãe
                'HG00096',  # GBR
                'HG00097',  # GBR
                'HG00099',  # GBR
            ],
            # GIAB (Genome in a Bottle) - gold standard
            'GIAB': [
                'HG001',  # NA12878
                'HG002',  # Filho Ashkenazi
                'HG003',  # Pai Ashkenazi
                'HG004',  # Mãe Ashkenazi
            ]
        }

    def download_1000g_bam(self, sample_id):
        """
        Baixa BAM do 1000 Genomes Project
        """
        logger.info(f"Downloading 1000G sample: {sample_id}")

        base_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"

        bam_file = self.output_dir / f"{sample_id}.bam"

        if bam_file.exists():
            logger.info(f"  {sample_id} já existe, pulando download")
            return str(bam_file)

        # URL específica pode variar, esse é um template
        # Na prática, você precisa consultar o catálogo do 1000G
        url = f"{base_url}{sample_id}/{sample_id}.alt_bwamem_GRCh38DH.20150826.*.bam"

        cmd = f"wget -P {self.output_dir} {url}"

        try:
            subprocess.run(cmd, shell=True, check=True)
            logger.info(f"  ✓ Download completo: {sample_id}")
            return str(bam_file)
        except subprocess.CalledProcessError:
            logger.warning(f"  ✗ Falha no download: {sample_id}")
            return None

    def download_giab_bam(self, sample_id):
        """
        Baixa BAM do Genome in a Bottle
        """
        logger.info(f"Downloading GIAB sample: {sample_id}")

        base_url = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"

        bam_file = self.output_dir / f"{sample_id}.bam"

        if bam_file.exists():
            logger.info(f"  {sample_id} já existe, pulando download")
            return str(bam_file)

        # Mapear sample_id para path específico
        giab_paths = {
            'HG001': 'NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_AJtrio_novoalign_bams/',
            'HG002': 'AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/',
            'HG003': 'AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/',
            'HG004': 'AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/'
        }

        if sample_id not in giab_paths:
            logger.error(f"Sample {sample_id} não encontrado no GIAB")
            return None

        # Note: URLs exatas devem ser verificadas no site do GIAB
        url = f"{base_url}{giab_paths[sample_id]}"

        logger.info(f"  URL: {url}")
        logger.info(f"  ⚠️  ATENÇÃO: Verifique a URL exata no site do GIAB")
        logger.info(f"     https://www.nist.gov/programs-projects/genome-bottle")

        # Para demonstração, não vamos baixar automaticamente
        logger.info(f"  Para baixar manualmente:")
        logger.info(f"     wget -P {self.output_dir} [URL_COMPLETA_DO_BAM]")

        return None

    def use_local_normals(self, bam_list_file):
        """
        Usa BAMs normais locais (se já tiver amostras normais disponíveis)
        """
        logger.info("Usando amostras normais locais")

        with open(bam_list_file) as f:
            bams = [line.strip() for line in f if line.strip()]

        return bams

    def create_cnvkit_reference(self, normal_bams, output_cnn, n_processes=8):
        """
        Cria referência CNVkit a partir dos BAMs normais
        """
        logger.info(f"Criando referência CNVkit com {len(normal_bams)} amostras normais")

        normal_bams_str = ' '.join(normal_bams)

        # Comando CNVkit batch para criar referência
        cmd = f"""
        cnvkit.py batch \
            --normal {normal_bams_str} \
            --targets {self.targets_bed} \
            --fasta {self.ref_genome} \
            --output-reference {output_cnn} \
            --output-dir {self.output_dir} \
            --processes {n_processes} \
            --method wgs
        """

        logger.info("Executando CNVkit...")
        logger.info(f"Comando: {cmd}")

        try:
            subprocess.run(cmd, shell=True, check=True)
            logger.info(f"✓ Referência criada com sucesso: {output_cnn}")
            return output_cnn
        except subprocess.CalledProcessError as e:
            logger.error(f"✗ Erro ao criar referência: {e}")
            return None

    def create_flat_reference(self, output_cnn):
        """
        Cria referência flat (fallback quando não há normais)
        Menos ideal, mas funciona
        """
        logger.warning("Criando referência FLAT (sem amostras normais)")
        logger.warning("⚠️  Isso aumentará falsos positivos. Recomendado usar PoN real.")

        cmd = f"""
        cnvkit.py batch \
            --targets {self.targets_bed} \
            --fasta {self.ref_genome} \
            --output-reference {output_cnn} \
            --output-dir {self.output_dir}
        """

        try:
            subprocess.run(cmd, shell=True, check=True)
            logger.info(f"✓ Referência flat criada: {output_cnn}")
            return output_cnn
        except subprocess.CalledProcessError as e:
            logger.error(f"✗ Erro ao criar referência flat: {e}")
            return None


def main():
    parser = argparse.ArgumentParser(
        description="Cria Panel of Normals para CNVkit"
    )
    parser.add_argument(
        '--output',
        default='reference/cnvkit_pon.cnn',
        help='Arquivo de saída da referência CNVkit'
    )
    parser.add_argument(
        '--targets',
        required=True,
        help='Arquivo BED com regiões alvo (exons)'
    )
    parser.add_argument(
        '--ref-genome',
        required=True,
        help='Genoma de referência FASTA'
    )
    parser.add_argument(
        '--normal-bams',
        help='Arquivo com lista de BAMs normais (um por linha)'
    )
    parser.add_argument(
        '--download-1000g',
        action='store_true',
        help='Baixar amostras do 1000 Genomes'
    )
    parser.add_argument(
        '--flat-reference',
        action='store_true',
        help='Criar referência flat (sem normais) - menos preciso'
    )
    parser.add_argument(
        '--processes',
        type=int,
        default=8,
        help='Número de processos paralelos'
    )

    args = parser.parse_args()

    builder = PanelOfNormalsBuilder(
        ref_genome=args.ref_genome,
        targets_bed=args.targets
    )

    normal_bams = []

    # Opção 1: Usar BAMs locais
    if args.normal_bams:
        logger.info("Usando BAMs normais locais")
        normal_bams = builder.use_local_normals(args.normal_bams)

    # Opção 2: Baixar do 1000G (opcional, pode ser demorado)
    elif args.download_1000g:
        logger.info("Baixando amostras do 1000 Genomes")
        for sample in builder.normal_samples['1000G'][:3]:  # Limita a 3 para teste
            bam = builder.download_1000g_bam(sample)
            if bam:
                normal_bams.append(bam)

    # Opção 3: Referência flat (fallback)
    elif args.flat_reference:
        logger.warning("Criando referência FLAT (sem controles normais)")
        builder.create_flat_reference(args.output)
        return

    # Se temos BAMs normais, cria PoN
    if normal_bams:
        builder.create_cnvkit_reference(
            normal_bams,
            args.output,
            n_processes=args.processes
        )
    else:
        logger.error("Nenhum BAM normal disponível!")
        logger.info("Opções:")
        logger.info("  1. Use --normal-bams com lista de BAMs locais")
        logger.info("  2. Use --flat-reference (menos preciso)")
        logger.info("  3. Baixe manualmente BAMs do 1000G/GIAB")
        sys.exit(1)


if __name__ == '__main__':
    main()
