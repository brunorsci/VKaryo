#!/bin/bash
# =============================================================
# check_annovar.sh - Verifica instalação local do ANNOVAR
# =============================================================
# Uso: bash scripts/check_annovar.sh [annovar_dir] [buildver]
# Exemplo: bash scripts/check_annovar.sh resources/annovar hg38
#
# Verifica:
# 1. Scripts Perl do ANNOVAR existem e são executáveis
# 2. Diretório humandb existe
# 3. Databases necessários estão presentes
# =============================================================

set -euo pipefail

ANNOVAR_DIR="${1:-resources/annovar}"
BUILDVER="${2:-hg38}"

# Cores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

ERRORS=0

echo "============================================================"
echo " ANNOVAR Installation Check"
echo "============================================================"
echo ""
echo "ANNOVAR directory: $ANNOVAR_DIR"
echo "Build version: $BUILDVER"
echo ""

# 1. Verifica scripts Perl
echo "--- Perl Scripts ---"
for script in table_annovar.pl convert2annovar.pl annotate_variation.pl; do
    if [ -f "$ANNOVAR_DIR/$script" ]; then
        if [ -r "$ANNOVAR_DIR/$script" ]; then
            echo -e "${GREEN}[OK]${NC} $script"
        else
            echo -e "${RED}[FAIL]${NC} $script (exists but not readable)"
            ERRORS=$((ERRORS + 1))
        fi
    else
        echo -e "${RED}[FAIL]${NC} $script (not found)"
        ERRORS=$((ERRORS + 1))
    fi
done
echo ""

# 2. Verifica humandb
echo "--- Database Directory ---"
if [ -d "$ANNOVAR_DIR/humandb" ]; then
    echo -e "${GREEN}[OK]${NC} humandb directory exists"
else
    echo -e "${RED}[FAIL]${NC} humandb directory not found"
    ERRORS=$((ERRORS + 1))
fi
echo ""

# 3. Verifica databases necessários para o workflow
echo "--- Required Databases ($BUILDVER) ---"

# Databases usados no workflow (mutect2_annovar + vardict_annovar)
REQUIRED_DBS=(
    "refGene:g"
    "cytoBand:r"
    "clinvar_20240917:f"
    "cosmic101:f"
    "gnomad41_genome:f"
    "dbnsfp47a:f"
    "avsnp154:f"
    "exac03:f"
    "dbscsnv11:f"
)

for db_entry in "${REQUIRED_DBS[@]}"; do
    db_name="${db_entry%%:*}"
    db_op="${db_entry##*:}"

    # gene-based (g) usa refGene.txt, filter-based (f) usa .txt + .idx
    db_file="${BUILDVER}_${db_name}.txt"

    if [ -f "$ANNOVAR_DIR/humandb/$db_file" ]; then
        size=$(du -sh "$ANNOVAR_DIR/humandb/$db_file" | cut -f1)
        echo -e "${GREEN}[OK]${NC} $db_name ($size)"

        # filter-based precisa de .idx
        if [ "$db_op" = "f" ]; then
            if [ -f "$ANNOVAR_DIR/humandb/${db_file}.idx" ]; then
                echo -e "${GREEN}[OK]${NC}   └── ${db_file}.idx"
            else
                echo -e "${RED}[FAIL]${NC}   └── ${db_file}.idx (missing index)"
                ERRORS=$((ERRORS + 1))
            fi
        fi
    else
        echo -e "${RED}[FAIL]${NC} $db_name (${db_file} not found)"
        ERRORS=$((ERRORS + 1))
    fi
done

# refGene também precisa do refGeneMrna.fa
echo ""
echo "--- Additional Required Files ---"
if [ -f "$ANNOVAR_DIR/humandb/${BUILDVER}_refGeneMrna.fa" ]; then
    echo -e "${GREEN}[OK]${NC} refGeneMrna.fa"
else
    echo -e "${YELLOW}[WARN]${NC} refGeneMrna.fa (optional, needed for mRNA annotation)"
fi

echo ""
echo "============================================================"
if [ $ERRORS -eq 0 ]; then
    echo -e "${GREEN}ANNOVAR installation OK - all checks passed${NC}"
    echo "============================================================"
    exit 0
else
    echo -e "${RED}ANNOVAR installation INCOMPLETE - $ERRORS error(s) found${NC}"
    echo "============================================================"
    exit 1
fi
