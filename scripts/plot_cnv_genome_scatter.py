#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cariotipo virtual -- todos os 24 cromossomos de UMA amostra, num grid, no
mesmo estilo do plot_cnv_scatter_colored.py (pontos individuais por bin do
CNVkit + segmento colorido por estado: azul=ganho, vermelho=delecao). Feito
para enxergar eventos GRANDES (deleção de braço, monossomia, trissomia) de
relance -- por isso o corte de cor é passado via --gain-threshold/
--loss-threshold e vem de config.yaml:cnv.karyotype_gain_threshold/
karyotype_loss_threshold, um parâmetro PRÓPRIO deste plot (não o
min_log2_ratio_amp/max_log2_ratio_del usado no achado focal por gene em
plot_cnv_scatter_colored.py nem no ideograma) -- o objetivo aqui é alteração
estrutural/cromossômica, não CNV pontual, então o corte certo pode (e deve)
ser diferente. Valor pré-definido = regime "germline samples" do próprio
CNVkit (cnvlib/call.py), ainda pendente de validação contra amostras com
cariótipo anotado -- ver comentário em config.yaml.

Uso:
    python3 plot_cnv_genome_scatter.py \
        --cnr sample.cnr --cns sample.cns \
        --gain-threshold 0.3 --loss-threshold -0.4 \
        --output saida.png [--title "titulo"] \
        [--gene-coords genes.tsv --gene GENE1 --gene GENE2 ...]

--gene/--gene-coords são opcionais -- se passados, destaca os genes de
interesse (mesma lógica do script por gene: rótulo ancorado no segmento real).
"""
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# mesmos comprimentos hg38 já usados em plot_aneuploidy_ideogram.py
CHR_LENGTHS = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
}

GAIN_COLOR = "#1f77b4"
LOSS_COLOR = "#d62728"
NEUTRAL_COLOR = "#888888"
POINT_COLOR = "#555555"
GENE_LINE_COLOR = "#d4a017"
Y_MIN, Y_MAX = -3, 3   # janela mais estreita que o plot por gene -- aqui o alvo e evento GRANDE (braco/cromossomo), nao pico focal


def load_gene_coords(path):
    coords = {}
    with open(path) as f:
        for line in f:
            g, c, s, e = line.rstrip("\n").split("\t")
            coords[g] = (c, int(s), int(e))
    return coords


def segment_color(log2, gain_threshold, loss_threshold):
    if log2 >= gain_threshold:
        return GAIN_COLOR
    if log2 <= loss_threshold:
        return LOSS_COLOR
    return NEUTRAL_COLOR


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cnr", required=True)
    ap.add_argument("--cns", required=True)
    ap.add_argument("--gain-threshold", type=float, required=True)
    ap.add_argument("--loss-threshold", type=float, required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--gene-coords", default=None)
    ap.add_argument("--gene", action="append", default=[])
    args = ap.parse_args()

    cnr_all = pd.read_csv(args.cnr, sep="\t")
    cns_all = pd.read_csv(args.cns, sep="\t")
    gene_coords = load_gene_coords(args.gene_coords) if args.gene_coords else {}
    genes_by_chrom = {}
    for g in args.gene:
        if g in gene_coords:
            c, s, e = gene_coords[g]
            genes_by_chrom.setdefault(c, []).append((g, s, e))

    ncols, nrows = 4, 6  # 24 cromossomos, grid 4x6
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 22))
    axes_flat = axes.flatten()

    for idx, chrom in enumerate(CHR_ORDER):
        ax = axes_flat[idx]
        chr_len_mb = CHR_LENGTHS[chrom] / 1e6

        cnr = cnr_all[cnr_all["chromosome"] == chrom]
        cns = cns_all[cns_all["chromosome"] == chrom]

        if len(cnr) > 0:
            mid = (cnr["start"] + cnr["end"]) / 2 / 1e6
            ax.scatter(mid, cnr["log2"], s=4, c=POINT_COLOR, alpha=0.25, linewidths=0)

        for _, seg in cns.iterrows():
            ax.plot([seg["start"] / 1e6, seg["end"] / 1e6], [seg["log2"], seg["log2"]],
                    color=segment_color(seg["log2"], args.gain_threshold, args.loss_threshold),
                    linewidth=2.5, solid_capstyle="butt", zorder=5)

        for gene, gs, ge in genes_by_chrom.get(chrom, []):
            ax.axvspan(gs / 1e6, ge / 1e6, color=GENE_LINE_COLOR, alpha=0.5, zorder=1)
            ax.text((gs + ge) / 2 / 1e6, Y_MAX - 0.3, gene, ha="center", va="top", fontsize=7,
                    fontweight="bold", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.7))

        ax.axhline(0, color="black", linewidth=0.5, zorder=2)
        ax.set_xlim(0, chr_len_mb)
        ax.set_ylim(Y_MIN, Y_MAX)
        ax.set_title(chrom.replace("chr", "Chr "), fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=7)
        if idx % ncols != 0:
            ax.set_yticklabels([])

    # esconde paineis vazios (24 cromossomos em grid 4x6=24, deve encaixar exato)
    for idx in range(len(CHR_ORDER), nrows * ncols):
        axes_flat[idx].axis("off")

    legend_elements = [
        plt.Line2D([0], [0], color=GAIN_COLOR, lw=3, label=f"Ganho (log2 ≥ {args.gain_threshold})"),
        plt.Line2D([0], [0], color=LOSS_COLOR, lw=3, label=f"Deleção (log2 ≤ {args.loss_threshold})"),
        plt.Line2D([0], [0], color=NEUTRAL_COLOR, lw=3, label="Neutro"),
    ]
    fig.legend(handles=legend_elements, loc="upper right", fontsize=11, bbox_to_anchor=(0.99, 0.995))
    fig.suptitle(args.title or "Cariótipo virtual", fontsize=16, fontweight="bold", y=0.998)
    fig.supxlabel("Posição (Mb)", fontsize=12)
    fig.supylabel("Copy ratio (log2)", fontsize=12)

    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.97])
    plt.savefig(args.output, dpi=130)
    print(f"Salvo: {args.output}")


if __name__ == "__main__":
    main()
