#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
import os

# -----------------------------
# Configurações
# -----------------------------
cnr_file = sys.argv[1] if len(sys.argv) > 1 else "SAMPLE_ID_markduplicates_recal.cnr"
cns_file = sys.argv[2] if len(sys.argv) > 2 else "SAMPLE_ID_markduplicates_recal.cns"
chrom = sys.argv[3] if len(sys.argv) > 3 else "chr7"
cyto_file = sys.argv[4] if len(sys.argv) > 4 else "cytoBand.txt"   # arquivo UCSC

output_file = f"{os.path.splitext(os.path.basename(cnr_file))[0]}_{chrom}_plot.png"

# -----------------------------
# Ler dados CNVkit
# -----------------------------
cnr = pd.read_csv(cnr_file, sep="\t", comment="#")
cns = pd.read_csv(cns_file, sep="\t", comment="#")

# Detecta automaticamente a coluna de cromossomo
chrom_col_cnr = "chrom" if "chrom" in cnr.columns else "chromosome"
chrom_col_cns = "chrom" if "chrom" in cns.columns else "chromosome"

# Filtra pelo cromossomo escolhido
cnr_chrom = cnr[cnr[chrom_col_cnr] == chrom].copy()
cns_chrom = cns[cns[chrom_col_cns] == chrom].copy()

if cnr_chrom.empty:
    raise ValueError(f"Não foram encontrados dados para {chrom} em {cnr_file}")

# -----------------------------
# Função para plotar ideograma
# -----------------------------
def plot_chrom_ideogram(ax, chrom, cyto_file):
    bands = pd.read_csv(cyto_file, sep="\t", header=None,
                        names=["chrom", "start", "end", "name", "gieStain"])
    bands = bands[bands["chrom"] == chrom]

    stain_colors = {
        "gneg": "white",
        "gpos25": "lightgray",
        "gpos50": "gray",
        "gpos75": "darkgray",
        "gpos100": "black",
        "acen": "red",
        "stalk": "brown"
    }

    for _, row in bands.iterrows():
        ax.add_patch(Rectangle(
            (row["start"]/1e6, 0),           # posição inicial em Mb
            (row["end"]-row["start"])/1e6,   # largura
            1,                               # altura
            facecolor=stain_colors.get(row["gieStain"], "white"),
            edgecolor="black"
        ))

    ax.set_xlim(0, bands["end"].max()/1e6)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(f"Ideograma - {chrom}")

# -----------------------------
# Criar figura (3 painéis)
# -----------------------------
fig, axs = plt.subplots(3, 1, figsize=(14, 8), sharex=True,
                        gridspec_kw={"height_ratios": [3, 2, 0.5]})

# Scatter plot por bin
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

# Ideograma
plot_chrom_ideogram(axs[2], chrom, cyto_file)

plt.tight_layout()
plt.savefig(output_file, dpi=300)
print(f"Plot salvo em {output_file}")
plt.show()
