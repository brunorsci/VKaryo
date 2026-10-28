#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

# -----------------------------
# Configurações
# -----------------------------
if len(sys.argv) < 4:
    print("Uso: python3 plot_cnv_all_chroms.py <CNR> <CNS> <diretorio_saida>")
    sys.exit(1)

cnr_file = sys.argv[1]
cns_file = sys.argv[2]
outdir = sys.argv[3]

os.makedirs(outdir, exist_ok=True)
pdf_file = os.path.join(outdir, os.path.splitext(os.path.basename(cnr_file))[0] + "_all_chroms.pdf")

# -----------------------------
# Ler dados CNVkit
# -----------------------------
cnr = pd.read_csv(cnr_file, sep="\t", comment="#")
cns = pd.read_csv(cns_file, sep="\t", comment="#")

chrom_col_cnr = "chrom" if "chrom" in cnr.columns else "chromosome"
chrom_col_cns = "chrom" if "chrom" in cns.columns else "chromosome"

chromosomes = sorted(cnr[chrom_col_cnr].unique())

# -----------------------------
# Criar PDF com todos os cromossomos
# -----------------------------
with PdfPages(pdf_file) as pdf:
    for chrom in chromosomes:
        cnr_chrom = cnr[cnr[chrom_col_cnr] == chrom].copy()
        cns_chrom = cns[cns[chrom_col_cns] == chrom].copy()

        if cnr_chrom.empty:
            continue

        fig, axs = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
        
        # Scatter plot
        colors = cnr_chrom["log2"].apply(lambda x: "red" if x > 0 else "blue")
        axs[0].scatter(cnr_chrom["start"]/1e6, cnr_chrom["log2"], c=colors, s=3, label="log2 ratio")
        axs[0].axhline(0, color='black', linestyle='--')
        axs[0].set_ylabel("log2 ratio")
        axs[0].set_title(f"{chrom} - Scatter plot por bin")
        axs[0].legend(loc="upper right")

        # Segment plot
        if not cns_chrom.empty:
            axs[1].hlines(cns_chrom["log2"], cns_chrom["start"]/1e6, cns_chrom["end"]/1e6,
                          colors='orange', linewidth=2)
        axs[1].axhline(0, color='black', linestyle='--')
        axs[1].set_xlabel("Position (Mb)")
        axs[1].set_ylabel("log2 ratio")
        axs[1].set_title(f"{chrom} - Segmentos CNVkit")

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print(f"Todos os plots salvos em PDF: {pdf_file}")

#conda activate cnvkit_env
#python3 scripts/plot_cnv_all_chroms.py \
#    SAMPLE_ID_markduplicates_recal.cnr \
#    SAMPLE_ID_markduplicates_recal.cns \
#    results/plots
