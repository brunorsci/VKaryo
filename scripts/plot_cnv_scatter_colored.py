#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scatter de CNV por cromossomo inteiro, no estilo do `cnvkit.py scatter`, mas
com duas diferencas: (1) o segmento e colorido por estado -- azul = ganho,
vermelho = delecao, cinza = neutro -- em vez de uma cor unica (o cnvkit
nativo nao suporta isso, so uma cor fixa para todos os segmentos); (2) o
rotulo do gene de interesse e ancorado na altura do segmento que de fato
sobrepoe o gene, nao numa posicao fixa do grafico -- evita ambiguidade entre
o agrupamento real de sondas e um ponto isolado de ruido em outra altura.

Uso:
    python3 plot_cnv_scatter_colored.py \
        --cnr sample.cnr --cns sample.cns --chrom chr18 \
        --gene BCL2 --gene-coords gene_coords_all.tsv \
        --output saida.png [--title "titulo customizado"]

Pode receber mais de um --gene (todos precisam estar no mesmo --chrom).
"""
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GAIN_COLOR = "#1f77b4"      # azul -- ganho
LOSS_COLOR = "#d62728"      # vermelho -- delecao
NEUTRAL_COLOR = "#888888"
POINT_COLOR = "#555555"
GENE_LINE_COLOR = "#d4a017"
Y_MIN, Y_MAX = -5, 5        # mesmo corte automatico do cnvkit.py scatter (exclui outliers extremos tipo -20, que sao ruido de contagem baixa, nao sinal linear)
# Sem valor default hardcoded aqui de proposito -- os cortes de ganho/perda
# tem que vir de --gain-threshold/--loss-threshold (config.yaml:
# cnv.min_log2_ratio_amp / cnv.max_log2_ratio_del), pra usar o mesmo padrao
# ja estabelecido no resto do pipeline (plot_aneuploidy_ideogram.py), nao um
# numero inventado so pra este script.


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


def label_y_for_gene(cns_chrom, gene_start, gene_end):
    """Log2 do segmento que sobrepoe o gene (o de maior sobreposicao, se
    houver mais de um) -- usado para ancorar o rotulo na altura certa, em
    vez de uma posicao fixa do grafico."""
    overlapping = cns_chrom[(cns_chrom["start"] < gene_end) & (cns_chrom["end"] > gene_start)]
    if overlapping.empty:
        return None
    overlapping = overlapping.copy()
    overlapping["overlap_bp"] = (
        overlapping[["end"]].min(axis=1).clip(upper=gene_end)
        - overlapping[["start"]].max(axis=1).clip(lower=gene_start)
    )
    best = overlapping.loc[overlapping["overlap_bp"].idxmax()]
    return float(best["log2"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cnr", required=True, help="Arquivo .cnr (log2 por bin)")
    ap.add_argument("--cns", required=True, help="Arquivo .cns (segmentos)")
    ap.add_argument("--chrom", required=True, help="Cromossomo a plotar (ex: chr18)")
    ap.add_argument("--gene", action="append", required=True, help="Gene(s) de interesse a destacar (repetir para mais de um)")
    ap.add_argument("--gene-coords", required=True, help="TSV gene\\tchrom\\tstart\\tend")
    ap.add_argument("--gain-threshold", type=float, required=True,
                     help="Log2 mínimo para colorir segmento como ganho (mesmo valor de config.yaml:cnv.min_log2_ratio_amp)")
    ap.add_argument("--loss-threshold", type=float, required=True,
                     help="Log2 máximo para colorir segmento como deleção (mesmo valor de config.yaml:cnv.max_log2_ratio_del)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    coords = load_gene_coords(args.gene_coords)

    cnr = pd.read_csv(args.cnr, sep="\t")
    cnr = cnr[cnr["chromosome"] == args.chrom].copy()
    cnr["mid"] = (cnr["start"] + cnr["end"]) / 2 / 1e6

    cns = pd.read_csv(args.cns, sep="\t")
    cns_chrom = cns[cns["chromosome"] == args.chrom].copy()

    fig, ax = plt.subplots(figsize=(9, 6))

    # pontos individuais por bin (cinza, sem filtro alem do corte de eixo Y)
    ax.scatter(cnr["mid"], cnr["log2"], s=8, c=POINT_COLOR, alpha=0.3, linewidths=0)

    # segmentos, coloridos por estado
    for _, seg in cns_chrom.iterrows():
        ax.plot([seg["start"] / 1e6, seg["end"] / 1e6], [seg["log2"], seg["log2"]],
                color=segment_color(seg["log2"], args.gain_threshold, args.loss_threshold),
                linewidth=3, solid_capstyle="butt", zorder=5)

    # genes de interesse -- faixa vertical + rotulo ancorado no segmento real
    plotted_genes = []
    for gene in args.gene:
        if gene not in coords:
            continue
        g_chrom, gs, ge = coords[gene]
        if g_chrom != args.chrom:
            continue
        ax.axvspan(gs / 1e6, ge / 1e6, color=GENE_LINE_COLOR, alpha=0.4, zorder=1)

        anchor_y = label_y_for_gene(cns_chrom, gs, ge)
        if anchor_y is None:
            anchor_y = 0.0  # sem segmento sobrepondo (raro) -- cai no neutro como ultimo recurso
        # desloca o texto um pouco acima do segmento, ou abaixo se ja estiver perto do teto
        label_y = anchor_y - 0.5 if anchor_y > Y_MAX - 1 else anchor_y + 0.3
        label_y = max(Y_MIN + 0.2, min(Y_MAX - 0.2, label_y))
        ax.text((gs + ge) / 2 / 1e6, label_y, gene, ha="center", fontsize=9,
                fontweight="bold", zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.75))
        plotted_genes.append(gene)

    ax.axhline(0, color="black", linewidth=0.8, zorder=2)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_xlabel("Posição (Mb)")
    ax.set_ylabel("Copy ratio (log2)")
    ax.set_title(args.title or f"{args.chrom} -- {', '.join(plotted_genes) if plotted_genes else args.chrom}")

    legend_elements = [
        plt.Line2D([0], [0], color=GAIN_COLOR, lw=3, label=f"Ganho (log2 ≥ {args.gain_threshold})"),
        plt.Line2D([0], [0], color=LOSS_COLOR, lw=3, label=f"Deleção (log2 ≤ {args.loss_threshold})"),
        plt.Line2D([0], [0], color=NEUTRAL_COLOR, lw=3, label="Neutro"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Salvo: {args.output}")


if __name__ == "__main__":
    main()
