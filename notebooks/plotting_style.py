"""Shared figure style and builders for the AON snRNA-seq notebooks.

One import point so every notebook produces consistent figures. The conventions follow
Rougier et al. 2014, PLOS Comput Biol; Wong 2010, Nat Methods; and Guha et al. 2022, eLife.
Chartjunk is removed, colour is colourblind-safe and redundant with direct labels, and
continuous scales are perceptually uniform.

Call apply_style() once at the top of a notebook, build figures with the helpers here,
and save with save_fig() to emit a raster PNG.
"""

from __future__ import annotations

import re
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

CELLTYPE_COLORS = {
    "Excitatory": "#0072B2",
    "Inhibitory": "#D55E00",
    "Unassigned": "#999999",
}

DIRECTION_COLORS = {
    "dorsolateral": "#0072B2",
    "ventromedial": "#D55E00",
    "posterior": "#0072B2",
    "anterior": "#D55E00",
}

SEQUENTIAL = "viridis"


def cluster_palette(n):
    """Qualitative colours for Leiden clusters, drawn from the matplotlib tab20 family.

    The muted tab20/tab20b/tab20c ramps stay visually distinct without oversaturation, and
    support up to 60 clusters. Identity is also carried by the on-data cluster numbers, so
    colour is a secondary cue.
    """
    import matplotlib.cm as cm

    base = ([mpl.colors.to_hex(cm.tab20(i)) for i in range(20)]
            + [mpl.colors.to_hex(cm.tab20b(i)) for i in range(20)]
            + [mpl.colors.to_hex(cm.tab20c(i)) for i in range(20)])
    if n > len(base):
        raise ValueError(f"cluster_palette supports up to {len(base)} clusters, got {n}")
    return base[:n]


def apply_style():
    """Override matplotlib defaults with the project figure conventions.

    Sets readable font sizes, hides the top and right spines, and raises the raster
    resolution. Figures are saved as PNG (see save_fig).
    """
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def save_fig(fig, path_stem, close=True):
    """Save a figure as a 300-dpi PNG under the given stem."""
    stem = Path(path_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    fig.savefig(png)
    if close:
        plt.close(fig)
    return png


def wrap(label, width=34):
    """Wrap a long label across lines instead of truncating it."""
    return textwrap.fill(str(label), width=width)


def clean_term(term):
    """Strip Enrichr GO and Reactome accession suffixes from a term label."""
    term = re.sub(r"\s*\(GO:\d+\)$", "", str(term))
    term = re.sub(r"\s*R-HSA-\d+$", "", term)
    return term.strip()


def direction_from_lfc(lfc, pos_label="posterior", neg_label="anterior"):
    """Map a log2 fold change sign to a group label (positive -> pos_label).

    Defaults keep the posterior/anterior robustness figures unchanged. The primary
    dorsolateral figures pass pos_label="dorsolateral", neg_label="ventromedial".
    """
    return np.where(np.asarray(lfc) > 0, pos_label, neg_label)


def volcano(res, fdr=0.05, lfc=1.0, highlight=None, title="Differential expression",
            symbol_col="symbol", pos_label="posterior", neg_label="anterior",
            figsize=(7, 6)):
    """Volcano plot coloured by direction, with labelled thresholds and highlights.

    res needs columns log2FoldChange, padj and a gene-symbol column. Points passing
    both thresholds are coloured by direction (pos_label for positive log2FC, neg_label
    for negative). Everything else is grey. Genes named in highlight are annotated where
    present.
    """
    d = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    d["neglog10padj"] = -np.log10(d["padj"].clip(lower=np.finfo(float).tiny))
    sig = (d["padj"] < fdr) & (d["log2FoldChange"].abs() > lfc)
    dirn = direction_from_lfc(d["log2FoldChange"], pos_label, neg_label)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(d.loc[~sig, "log2FoldChange"], d.loc[~sig, "neglog10padj"],
               s=6, c="0.8", linewidths=0, zorder=1)
    for grp in (pos_label, neg_label):
        m = sig & (dirn == grp)
        ax.scatter(d.loc[m, "log2FoldChange"], d.loc[m, "neglog10padj"],
                   s=10, c=DIRECTION_COLORS[grp], linewidths=0, zorder=2,
                   label=f"up in {grp} ({int(m.sum())})")

    ax.axhline(-np.log10(fdr), ls="--", color="0.55", lw=1)
    for x in (-lfc, lfc):
        ax.axvline(x, ls="--", color="0.55", lw=1)
    # Threshold labels along the bottom axis, clear of the legend and gene labels.
    ax.text(0.985, -np.log10(fdr) + 0.4, f"FDR {fdr}", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=8, color="0.5")
    ax.text(lfc, 0.015, f"|log2FC| = {lfc:g}", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=8, color="0.5")

    if highlight is not None:
        present = d[d[symbol_col].isin(highlight) & sig]
        for _, row in present.iterrows():
            ax.annotate(row[symbol_col], (row["log2FoldChange"], row["neglog10padj"]),
                        fontsize=8, fontstyle="italic",
                        xytext=(4, 2), textcoords="offset points")

    ax.set_xlabel(f"log2 fold change ({pos_label} / {neg_label})")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(title, pad=12)
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="0.85",
              framealpha=0.9)
    return fig


def candidate_effectsize(cand, fdr=0.05, title="Candidate markers in the DE result",
                         pos_label="posterior", neg_label="anterior", figsize=(7, 8)):
    """Horizontal effect-size plot of candidate genes ranked by log2 fold change.

    cand needs columns symbol, log2FoldChange, padj. Points are coloured by direction
    (pos_label for positive log2FC). A filled marker means it passed the FDR cutoff, an
    open marker means it did not.
    """
    d = cand.dropna(subset=["log2FoldChange"]).sort_values("log2FoldChange").copy()
    d["direction"] = direction_from_lfc(d["log2FoldChange"], pos_label, neg_label)
    sigmask = d["padj"].fillna(1.0) < fdr
    y = np.arange(len(d))

    fig, ax = plt.subplots(figsize=figsize)
    ax.axvline(0, color="0.6", lw=1)
    for grp in (pos_label, neg_label):
        for sig in (True, False):
            m = (d["direction"] == grp) & (sigmask if sig else ~sigmask)
            if not m.any():
                continue
            ax.scatter(d.loc[m, "log2FoldChange"], y[m.values], s=55, zorder=3,
                       color=DIRECTION_COLORS[grp],
                       facecolor=DIRECTION_COLORS[grp] if sig else "white",
                       edgecolor=DIRECTION_COLORS[grp], linewidth=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels(d["symbol"], fontstyle="italic", fontsize=9)
    ax.set_xlabel(f"log2 fold change ({pos_label} / {neg_label})")
    ax.set_title(title)

    handles = [
        Line2D([0], [0], marker="o", ls="", color=DIRECTION_COLORS[grp], label=f"{grp}-up")
        for grp in (pos_label, neg_label) if (d["direction"] == grp).any()
    ]
    handles += [
        Line2D([0], [0], marker="o", ls="", markerfacecolor="0.4",
               markeredgecolor="0.4", color="w", label=f"FDR < {fdr}"),
        Line2D([0], [0], marker="o", ls="", markerfacecolor="white",
               markeredgecolor="0.4", color="w", label="not significant"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8)
    return fig


def projection_panel(proj, fdr=0.05, title="Canonical cortical projection markers",
                     pos_label="posterior", neg_label="anterior", figsize=(7, 5)):
    """Diverging bar chart of projection-marker log2 fold changes, coloured by direction.

    proj needs columns symbol and log2FoldChange (optionally padj and a class column).
    Bars clearing FDR < fdr are solid, non-significant bars are hatched and unfilled, so a
    null-effect marker reads differently from a supported one. Makes visible that the markers
    split both ways, i.e. no coherent projection-class signature.
    """
    d = proj.dropna(subset=["log2FoldChange"]).sort_values("log2FoldChange").copy()
    d["direction"] = direction_from_lfc(d["log2FoldChange"], pos_label, neg_label)
    sig = (d["padj"].fillna(1.0) < fdr).to_numpy() if "padj" in d.columns else np.ones(len(d), bool)
    colors = d["direction"].map(DIRECTION_COLORS).to_numpy()
    y = np.arange(len(d))

    fig, ax = plt.subplots(figsize=figsize)
    for yi, (lfc, col, is_sig) in enumerate(zip(d["log2FoldChange"], colors, sig)):
        ax.barh(yi, lfc, color=col if is_sig else "white", edgecolor=col,
                hatch=None if is_sig else "////", linewidth=1.1, zorder=2)
    ax.axvline(0, color="0.5", lw=1)
    ax.set_yticks(y)
    # A real italic font slants the trailing digits too (mathtext \it leaves them upright),
    # so the gene symbols are fully italic. The class tag stays roman on a right-hand twin axis.
    ax.set_yticklabels(d["symbol"], fontstyle="italic", fontsize=9)
    ax.set_xlabel(f"log2 fold change ({pos_label} / {neg_label})")
    ax.set_title(title)

    if "class" in d.columns:
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(y)
        ax2.set_yticklabels(d["class"], fontsize=8, color="0.4")
        ax2.tick_params(length=0)
        ax2.spines["right"].set_visible(False)

    handles = [
        Line2D([0], [0], marker="s", ls="", markerfacecolor="0.4", markeredgecolor="0.4",
               color="w", label=f"FDR < {fdr}"),
        Line2D([0], [0], marker="s", ls="", markerfacecolor="white", markeredgecolor="0.4",
               color="w", label="not significant"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8, frameon=False)
    return fig


def ma_plot(res, fdr=0.05, lfc=1.0, title="MA plot",
            pos_label="posterior", neg_label="anterior", figsize=(7, 5)):
    """MA plot of log2 fold change against mean expression, significant points coloured."""
    d = res.dropna(subset=["baseMean", "log2FoldChange"]).copy()
    d = d[d["baseMean"] > 0]
    sig = (d["padj"].fillna(1.0) < fdr) & (d["log2FoldChange"].abs() > lfc)
    dirn = direction_from_lfc(d["log2FoldChange"], pos_label, neg_label)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(d.loc[~sig, "baseMean"], d.loc[~sig, "log2FoldChange"],
               s=5, c="0.8", linewidths=0)
    for grp in (pos_label, neg_label):
        m = sig & (dirn == grp)
        ax.scatter(d.loc[m, "baseMean"], d.loc[m, "log2FoldChange"],
                   s=7, c=DIRECTION_COLORS[grp], linewidths=0, label=f"up in {grp}")
    ax.axhline(0, color="0.4", lw=1)
    for h in (-lfc, lfc):
        ax.axhline(h, ls="--", color="0.5", lw=0.8)
    ax.text(0.99, lfc, f"log2FC = +/-{lfc:g}", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=8, color="0.5")
    ax.text(0.99, 0.02, f"coloured: FDR < {fdr:g} and |log2FC| > {lfc:g}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="0.5")
    ax.set_xscale("log")
    ax.set_xlabel("mean of normalized counts (log scale)")
    ax.set_ylabel(f"log2 fold change ({pos_label} / {neg_label})")
    ax.set_title(title)
    ax.legend(loc="upper right")
    return fig


def pseudobulk_pca(counts, coldata, title="Pseudobulk PCA", figsize=(7, 6)):
    """PCA of log-CPM pseudobulk profiles, coloured and shaped by group.

    counts is samples-by-genes. coldata is indexed by the same samples with group and
    donor columns. Each donor's two profiles are joined by a line, confirming visually that
    no single donor drives the contrast.
    """
    x = counts.reindex(coldata.index)
    cpm = x.div(x.sum(axis=1), axis=0) * 1e6
    logcpm = np.log1p(cpm)
    logcpm = logcpm.loc[:, logcpm.var(axis=0) > 0]
    z = (logcpm - logcpm.mean(axis=0)) / logcpm.std(axis=0)

    u, s, vt = np.linalg.svd(z.values, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    var = (s ** 2) / (s ** 2).sum() * 100

    fig, ax = plt.subplots(figsize=figsize)
    # Thin lines connect the two pseudobulk profiles of the same donor, so the pairing is
    # visible without per-point labels (which overcrowd at this sample size).
    for donor in pd.unique(coldata["donor"]):
        m = (coldata["donor"] == donor).values
        if m.sum() == 2:
            ax.plot(pcs[m, 0], pcs[m, 1], color="0.8", lw=0.8, zorder=1)
    markers = ["o", "^", "s", "D"]
    for i, grp in enumerate(coldata["group"].unique()):
        m = (coldata["group"] == grp).values
        ax.scatter(pcs[m, 0], pcs[m, 1], s=90, color=DIRECTION_COLORS.get(grp, "0.4"),
                   marker=markers[i % len(markers)], edgecolor="white", linewidth=0.8,
                   label=grp, zorder=3)
    ax.set_xlabel(f"PC1 ({var[0]:.0f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.0f}%)")
    ax.set_title(title, pad=12)
    ax.legend(title="group", loc="best", frameon=True, facecolor="white", edgecolor="0.85")
    return fig


def mean_variance(counts, title="Pseudobulk mean-variance trend", figsize=(7, 5)):
    """Per-gene mean vs variance of log-CPM across pseudobulk samples."""
    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    logcpm = np.log1p(cpm)
    mean = logcpm.mean(axis=0)
    var = logcpm.var(axis=0)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(mean, var, s=5, c="0.6", linewidths=0)
    ax.set_xlabel("mean ln(CPM + 1)")
    ax.set_ylabel("variance of ln(CPM + 1)")
    ax.set_title(title)
    return fig


def marker_dotplot(expr, groups, markers, title="Marker expression", figsize=None):
    """Compact marker dot plot with rows as groups and columns as markers.

    Each dot size is the fraction of cells in the group expressing the marker (count > 0),
    colour is the group mean expression scaled 0-1 within each marker. expr is a
    cells-by-genes DataFrame (log-normalized), groups is a length-n label per cell.
    """
    present = [m for m in markers if m in expr.columns]
    order = list(pd.unique(groups))
    groups = np.asarray(groups)
    frac = np.zeros((len(order), len(present)))
    mean = np.zeros((len(order), len(present)))
    for i, g in enumerate(order):
        gm = groups == g
        for j, mk in enumerate(present):
            vals = np.asarray(expr.loc[gm, mk])
            frac[i, j] = (vals > 0).mean() if gm.any() else 0.0
            mean[i, j] = vals.mean() if gm.any() else 0.0
    scaled = np.zeros_like(mean)
    for j in range(len(present)):
        col = mean[:, j]
        rng = col.max() - col.min()
        scaled[:, j] = (col - col.min()) / rng if rng > 0 else 0.0

    if figsize is None:
        figsize = (0.7 * len(present) + 4.2, 0.55 * len(order) + 2.2)
    fig, ax = plt.subplots(figsize=figsize)
    xs, ys = np.meshgrid(np.arange(len(present)), np.arange(len(order)))
    dot_area_scale, dot_area_min = 340, 10
    sizes = frac.ravel() * dot_area_scale + dot_area_min
    sc = ax.scatter(xs.ravel(), ys.ravel(), s=sizes, c=scaled.ravel(), cmap=SEQUENTIAL,
                    vmin=0, vmax=1, edgecolor="0.3", linewidth=0.4, zorder=3)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels(present, rotation=45, ha="right", fontstyle="italic")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlim(-0.6, len(present) - 0.4)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.invert_yaxis()
    ax.set_axisbelow(True)
    ax.grid(True, color="0.92", lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_color("0.8")
    ax.set_title(title, pad=10)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03)
    cbar.set_label("mean expression (scaled per marker)", fontsize=9)
    handles = [Line2D([0], [0], marker="o", ls="", markerfacecolor="0.6", markeredgecolor="0.3",
                      markersize=np.sqrt(f * dot_area_scale + dot_area_min),
                      label=f"{int(f * 100)}%")
               for f in (0.25, 0.5, 1.0)]
    ax.legend(handles=handles, title="fraction expressing", loc="center left",
              bbox_to_anchor=(1.30, 0.5), frameon=False, labelspacing=2.6,
              handletextpad=2.0, borderpad=0.8)
    fig.subplots_adjust(left=0.2, right=0.78, top=0.86, bottom=0.24)
    return fig


def enrichment_dotplot(anterior_hits, posterior_hits, n_anterior, n_posterior,
                       fdr=0.05, top_n=6, figsize=None):
    """Compact dot plot of over-representation terms, most significant at the top.

    One tightly packed column of terms across all libraries (library shown as a suffix),
    dot size = gene overlap, colour = odds ratio. The tested-but-null direction is stated in
    the caption. Enrichr adjusts p within each library, so ranks compare within a library.
    """
    library_short = {
        "GO_Biological_Process_2021": "GO",
        "KEGG_2019_Mouse": "KEGG",
        "Reactome_2022": "Reactome",
    }
    size_scale = 11

    if anterior_hits is None or not len(anterior_hits):
        fig, ax = plt.subplots(figsize=(7, 3))
        ax.axis("off")
        ax.text(0.5, 0.5, f"No terms reached FDR < {fdr}", ha="center", va="center", fontsize=11)
        return fig

    ant = (anterior_hits.sort_values(["Adjusted P-value", "Combined Score", "Term"],
                                     ascending=[True, False, True])
                        .groupby("Gene_set", observed=True).head(top_n).copy())
    ant["neglog_adjp"] = -np.log10(ant["Adjusted P-value"].clip(lower=np.finfo(float).tiny))
    ant["overlap_count"] = ant["Overlap"].map(lambda o: int(str(o).split("/")[0]))
    ant["lib"] = ant["Gene_set"].map(library_short).fillna("other")
    ant["label"] = [wrap(f"{clean_term(t)}  ({lib})", 38)
                    for t, lib in zip(ant["Term"], ant["lib"])]
    ant = ant.sort_values("neglog_adjp", ascending=True).reset_index(drop=True)

    n = len(ant)
    if figsize is None:
        figsize = (11.5, max(4.2, 0.52 * n + 2.2))
    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(n)
    norm = plt.Normalize(ant["Odds Ratio"].min(), ant["Odds Ratio"].max())
    fdr_line = -np.log10(fdr)
    scatter = ax.scatter(ant["neglog_adjp"], y, s=ant["overlap_count"] * size_scale,
                         c=ant["Odds Ratio"], cmap=SEQUENTIAL, norm=norm, edgecolor="0.25",
                         linewidth=0.4, zorder=3)
    ax.axvline(fdr_line, ls="--", color="0.55", lw=1, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(ant["label"], fontsize=9)
    ax.set_ylim(-1.15, n - 0.2)
    ax.text(fdr_line, -0.98, f"FDR {fdr}", ha="center", va="bottom", fontsize=8, color="0.5")
    ax.set_xlim(0, ant["neglog_adjp"].max() * 1.12)
    ax.set_xlabel("-log10 adjusted p-value", labelpad=8)
    ax.set_title(f"Over-representation of anterior-upregulated genes (n = {n_anterior})", pad=10)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="0.93", lw=0.6)

    # Both keys sit outside the plot on the right, the odds-ratio colour bar in the upper right and
    # the gene-overlap size key stacked below it (smallest dot nearest its title, so the large dot
    # never rides up into the label).
    cax = ax.inset_axes([1.05, 0.55, 0.03, 0.40])
    cbar = fig.colorbar(scatter, cax=cax)
    cbar.set_label("odds ratio", fontsize=9)
    size_ticks = [10, 30, 60]
    handles = [Line2D([0], [0], marker="o", ls="", markerfacecolor="0.6", markeredgecolor="0.25",
                      markersize=np.sqrt(c * size_scale), label=str(c)) for c in size_ticks]
    ax.legend(handles=handles, title="genes in overlap", loc="upper left",
              bbox_to_anchor=(1.02, 0.40), frameon=False, labelspacing=2.4,
              handletextpad=1.4, borderpad=0.8)
    fig.text(0.5, 0.02,
             f"Posterior direction: no terms reached FDR < {fdr} (n = {n_posterior} genes tested). "
             "Adjusted p is Benjamini-Hochberg within each library.",
             ha="center", fontsize=8, color="0.45")
    fig.subplots_adjust(left=0.36, right=0.82, top=0.92, bottom=0.12)
    return fig

