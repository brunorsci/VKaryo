#!/bin/bash
##############################################
# Download Reference Resources for Mutect2
# VKaryo Workflow Optimization
##############################################

set -e  # Exit on error

# Reference directory
REF_DIR="./reference"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

echo "============================================"
echo "Downloading Mutect2 Reference Resources"
echo "============================================"

# 1. gnomAD AF-only VCF (for germline filtering)
# Using GATK resource bundle - pre-built AF-only gnomAD for hg38
echo ""
echo "[1/3] Downloading gnomAD AF-only from GATK resource bundle..."
if [ ! -f "af-only-gnomad.hg38.vcf.gz" ]; then
    # Download from GATK resource bundle on Google Cloud
    wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
        -O af-only-gnomad.hg38.vcf.gz

    # Download index
    wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi \
        -O af-only-gnomad.hg38.vcf.gz.tbi

    echo "✓ gnomAD downloaded and indexed"
else
    echo "✓ gnomAD already exists"
fi

# 2. Panel of Normals from GATK Resource Bundle
echo ""
echo "[2/3] Downloading GATK 1000G Panel of Normals..."
if [ ! -f "leukemia_pon.vcf.gz" ]; then
    # Download 1000G PoN
    wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz \
        -O leukemia_pon.vcf.gz

    # Create index
    tabix -p vcf leukemia_pon.vcf.gz

    echo "✓ Panel of Normals downloaded and indexed"
else
    echo "✓ Panel of Normals already exists"
fi

# 3. COSMIC Hematologic Malignancies
echo ""
echo "[3/3] COSMIC Hematologic Database..."
echo ""
echo "IMPORTANT: COSMIC requires registration and authentication."
echo "Please download manually from:"
echo "  https://cancer.sanger.ac.uk/cosmic/download"
echo ""
echo "Required file: CosmicMutantExport_Haem.tsv.gz"
echo "Save as: ${REF_DIR}/cosmic_heme_v99.vcf.gz"
echo ""
echo "Alternatively, convert COSMIC TSV to VCF using:"
echo "  python scripts/cosmic_to_vcf.py CosmicMutantExport_Haem.tsv.gz cosmic_heme_v99.vcf"
echo "  bgzip cosmic_heme_v99.vcf"
echo "  tabix -p vcf cosmic_heme_v99.vcf.gz"

if [ -f "cosmic_heme_v99.vcf.gz" ]; then
    echo "✓ COSMIC Heme already exists"
else
    echo "⚠ COSMIC Heme not found - will need manual download"
    # Create empty placeholder for now
    touch cosmic_heme_v99.vcf.gz.placeholder
fi

echo ""
echo "============================================"
echo "Download Status:"
echo "============================================"
ls -lh af-only-gnomad.hg38.vcf.gz* 2>/dev/null || echo "  ✗ gnomAD missing"
ls -lh leukemia_pon.vcf.gz* 2>/dev/null || echo "  ✗ PoN missing"
ls -lh cosmic_heme_v99.vcf.gz* 2>/dev/null || echo "  ⚠ COSMIC missing (manual download required)"
echo "============================================"
echo ""
echo "Setup complete! Reference resources ready for Mutect2."
echo ""
