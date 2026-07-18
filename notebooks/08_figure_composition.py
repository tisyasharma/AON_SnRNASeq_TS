#!/usr/bin/env python
"""
08_figure_composition.py — Curated figure set for the AON snRNA-seq project.

Regenerates every figure in results/figures_final/ (main/ + supplement/) from
the committed inputs: data/aon_10x/AON_preprocessed.h5ad and results/tables/*.csv.

These are the hand-composed, presentation-ready figures — distinct from the raw
per-notebook panels in results/figures_raw/ produced by notebooks 01–07. Each
figure is one self-contained function; run this script top-to-bottom to rebuild
the whole set.

    python notebooks/08_figure_composition.py            # all figures
    python notebooks/08_figure_composition.py Fig4_DE     # one figure by name

Requires: scanpy/anndata, pydeseq2 (Fig 4 refits DESeq2 for LFC shrinkage),
matplotlib, colorcet, adjustText (Fig S2), scipy — see environment.yml.

Coverage: this script composes the 6 main figures, the 2 standalone UMAPs, and
supplement Fig S1 (QC) and Fig S2 (DE robustness). Supplement Fig S3–S7 are direct
per-notebook outputs copied into supplement/ (S3 DE diagnostics, S4 projection
markers, S5 candidate signature — from notebooks 04/02; S6 DecontX sensitivity,
S7 correlation-concordance heatmap — from the review analyses), not re-composed here.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- repo-relative paths (script lives in notebooks/) -----------------------
ROOT = Path(__file__).resolve().parent.parent
NB = ROOT / "notebooks"
H5AD = ROOT / "data" / "aon_10x" / "AON_preprocessed.h5ad"
TAB = ROOT / "results" / "tables"
OUT_MAIN = ROOT / "results" / "figures_final" / "main"
OUT_SUPPLEMENT = ROOT / "results" / "figures_final" / "supplement"
OUT_MAIN.mkdir(parents=True, exist_ok=True)
OUT_SUPPLEMENT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(NB))  # for plotting_style
import plotting_style as ps  # noqa: E402  (must follow sys.path setup above)

ps.apply_style()


def _umap_axes(ax, frac=0.16, fs=8):
    """Equal-aspect, tick-free embedding axes with corner UMAP1/UMAP2 arrows.

    Shared by every UMAP panel so the axis styling lives in one place.
    """
    ax.set_aspect("equal", "box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    L = frac * (x1 - x0)
    ax_x = x0 + 0.02 * (x1 - x0)
    ax_y = y0 + 0.02 * (y1 - y0)
    ax.annotate(
        "",
        xy=(ax_x + L, ax_y),
        xytext=(ax_x, ax_y),
        arrowprops=dict(arrowstyle="->", lw=1.0, color="k"),
        annotation_clip=False,
    )
    ax.annotate(
        "",
        xy=(ax_x, ax_y + L),
        xytext=(ax_x, ax_y),
        arrowprops=dict(arrowstyle="->", lw=1.0, color="k"),
        annotation_clip=False,
    )
    ax.text(ax_x + L / 2, ax_y - 0.03 * (y1 - y0), "UMAP1", ha="center", va="top", fontsize=fs)
    ax.text(
        ax_x - 0.03 * (x1 - x0),
        ax_y + L / 2,
        "UMAP2",
        ha="right",
        va="center",
        rotation=90,
        fontsize=fs,
    )


# ============================================================================
# Figure 1 — Single-nucleus atlas (4-panel: lineages, Leiden clusters, markers, composition)
# ============================================================================
def make_fig1_atlas():

    import anndata as ad
    import colorcet as cc
    import matplotlib.patheffects as pe
    import scipy.sparse as sp
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    A = ad.read_h5ad(str(H5AD), backed="r")
    U = np.asarray(A.obsm["X_umap"])
    leiden = np.asarray(A.obs["leiden"]).astype(str)
    lineage = np.asarray(A.obs["lineage"]).astype(str)
    var_names = np.asarray(A.var_names)
    markers = ["Slc17a7", "Slc17a6", "Gad1", "Gad2", "Slc32a1"]
    gidx = [int(np.where(var_names == g)[0][0]) for g in markers]
    Xm = A[:, gidx].X
    Xm = Xm.toarray() if sp.issparse(Xm) else np.asarray(Xm)
    clusters = sorted(set(leiden), key=int)

    PT = 1.0
    STROKE = [pe.withStroke(linewidth=1.6, foreground="white")]
    GLASBEY = cc.glasbey_hv[: len(clusters)]
    CMAP_CL = {c: GLASBEY[i] for i, c in enumerate(clusters)}

    style_umap = _umap_axes

    fig = plt.figure(figsize=(14, 13))
    outer = fig.add_gridspec(
        2, 2, hspace=0.14, wspace=0.10, left=0.055, right=0.965, top=0.905, bottom=0.055
    )
    specs = {"A": outer[0, 0], "B": outer[0, 1], "C": outer[1, 0], "D": outer[1, 1]}

    # A lineages
    axA = fig.add_subplot(specs["A"])
    axA.set_anchor("N")
    for lab in ["Inhibitory", "Excitatory", "Unassigned"]:
        m = lineage == lab
        axA.scatter(
            U[m, 0],
            U[m, 1],
            s=PT,
            color=ps.CELLTYPE_COLORS[lab],
            linewidths=0,
            rasterized=True,
            label=lab,
        )
    style_umap(axA)
    axA.legend(markerscale=6, loc="upper right", handletextpad=0.2, frameon=False)

    # B clusters
    axB = fig.add_subplot(specs["B"])
    axB.set_anchor("N")
    for c in clusters:
        m = leiden == c
        axB.scatter(U[m, 0], U[m, 1], s=PT, color=CMAP_CL[c], linewidths=0, rasterized=True)
        axB.text(
            U[m, 0].mean(),
            U[m, 1].mean(),
            c,
            fontsize=7,
            ha="center",
            va="center",
            color="black",
            path_effects=STROKE,
        )
    style_umap(axB)

    # C dotplot
    axC = fig.add_subplot(specs["C"])
    axC.set_anchor("N")
    lin_order = ["Excitatory", "Inhibitory", "Unassigned"]
    frac = np.zeros((3, 5))
    mean = np.zeros((3, 5))
    for i, lab in enumerate(lin_order):
        m = lineage == lab
        sub = Xm[m]
        frac[i] = (sub > 0).mean(0)
        mean[i] = sub.mean(0)
    sctC = None
    for i in range(3):
        for j in range(5):
            sctC = axC.scatter(
                j,
                2 - i,
                s=frac[i, j] * 300,
                c=mean[i, j],
                cmap="Reds",
                vmin=0,
                vmax=mean.max(),
                edgecolor="grey",
                linewidth=0.4,
            )
    axC.set_xticks(range(5))
    axC.set_xticklabels([f"$\\it{{{g}}}$" for g in markers], fontsize=10)
    axC.set_yticks(range(3))
    axC.set_yticklabels(lin_order[::-1], fontsize=10)
    axC.set_ylim(-0.9, 2.6)
    axC.set_xlim(-0.6, 4.6)
    for s_ in ["top", "right"]:
        axC.spines[s_].set_visible(False)
    caxC = inset_axes(
        axC,
        width="42%",
        height="4.5%",
        loc="lower left",
        bbox_to_anchor=(0.0, -0.13, 1, 1),
        bbox_transform=axC.transAxes,
        borderpad=0,
    )
    cb = fig.colorbar(sctC, cax=caxC, orientation="horizontal")
    cb.set_label("mean expression", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    axC.text(2.35, -0.55, "% expressing:", ha="right", va="center", fontsize=9)
    for xi, fr in zip([2.9, 3.5, 4.1], [0.25, 0.5, 1.0]):
        axC.scatter(xi, -0.55, s=fr * 300, c="lightgrey", edgecolor="grey")
        axC.text(xi, -0.55, f"{int(fr * 100)}", ha="center", va="center", fontsize=6)

    # D composition
    axD = fig.add_subplot(specs["D"])
    axD.set_anchor("N")
    counts = {lin: int((lineage == lin).sum()) for lin in lin_order}
    bars = axD.bar(
        lin_order,
        [counts[lin] for lin in lin_order],
        color=[ps.CELLTYPE_COLORS[lin] for lin in lin_order],
        edgecolor="k",
        linewidth=0.5,
        width=0.6,
    )
    for b, lin in zip(bars, lin_order):
        axD.text(
            b.get_x() + b.get_width() / 2,
            counts[lin] + 300,
            f"{counts[lin]:,}",
            ha="center",
            fontsize=11,
            fontweight="bold",
        )
    axD.set_ylabel("nuclei")
    axD.set_ylim(0, 26000)
    for s_ in ["top", "right"]:
        axD.spines[s_].set_visible(False)

    fig.canvas.draw()
    names = {
        "A": "Cell lineages",
        "B": "Leiden clusters (22)",
        "C": "Lineage markers",
        "D": "Composition",
    }
    for k, spec in specs.items():
        bb = spec.get_position(fig)
        fig.text(
            bb.x0,
            bb.y1 + 0.008,
            f"{k}  {names[k]}",
            ha="left",
            va="bottom",
            fontweight="bold",
            fontsize=13,
        )

    fig.suptitle(
        "Figure 1  |  Single-nucleus atlas of the mouse anterior olfactory nucleus (32,952 nuclei)",
        fontsize=15,
        fontweight="bold",
        y=0.965,
    )
    fig.savefig(str(OUT_MAIN / "Fig1_atlas.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", "Fig1_atlas.png")


# ============================================================================
# Standalone — Leiden cluster UMAP (glasbey_hv, bold labels)
# ============================================================================
def make_umap_leiden_standalone():

    import anndata as ad
    import colorcet as cc
    import matplotlib.patheffects as pe
    from matplotlib.lines import Line2D

    A = ad.read_h5ad(str(H5AD))
    U = A.obsm["X_umap"]
    leiden = np.asarray(A.obs["leiden"]).astype(str)
    clusters = sorted(set(leiden), key=int)

    GLASBEY = cc.glasbey_hv[: len(clusters)]
    CMAP_CL = {c: GLASBEY[i] for i, c in enumerate(clusters)}
    STROKE = [pe.withStroke(linewidth=1.8, foreground="white")]

    style_umap2 = _umap_axes

    fig, ax = plt.subplots(figsize=(9.2, 6.5))
    for c in clusters:
        m = leiden == c
        ax.scatter(U[m, 0], U[m, 1], s=1.0, color=CMAP_CL[c], linewidths=0, rasterized=True)
        ax.text(
            U[m, 0].mean(),
            U[m, 1].mean(),
            c,
            fontsize=7,
            fontweight="bold",
            ha="center",
            va="center",
            color="black",
            path_effects=STROKE,
        )
    style_umap2(ax)
    ax.set_title("Leiden clusters (resolution 0.1, silhouette-selected)", fontweight="bold")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", ls="", mfc=CMAP_CL[c], mec="none", ms=6, label=c)
            for c in clusters
        ],
        title="cluster",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        ncol=2,
        columnspacing=0.8,
        handletextpad=0.2,
        labelspacing=0.35,
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )
    fig.savefig(str(OUT_MAIN / "umap_leiden_clusters_standalone.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("  wrote", "umap_leiden_clusters_standalone.png")


# ============================================================================
# Standalone — cell-type (lineage) UMAP
# ============================================================================
def make_umap_celltypes_standalone():

    import anndata as ad

    A = ad.read_h5ad(str(H5AD))
    U = A.obsm["X_umap"]
    lineage = np.asarray(A.obs["lineage"]).astype(str)

    PT = 1.0

    style_umap = _umap_axes

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for lab in ["Inhibitory", "Excitatory", "Unassigned"]:
        m = lineage == lab
        ax.scatter(
            U[m, 0],
            U[m, 1],
            s=PT,
            c=ps.CELLTYPE_COLORS[lab],
            linewidths=0,
            rasterized=True,
            label=lab,
        )
    style_umap(ax)
    ax.set_title("Cell lineages", fontweight="bold")
    ax.legend(markerscale=6, loc="upper right", handletextpad=0.2, frameon=False)
    fig.savefig(str(OUT_MAIN / "umap_cell_types_standalone.png"), dpi=300, bbox_inches="tight")
    print("  wrote", "umap_cell_types_standalone.png")


# ============================================================================
# Figure 2 — Projection-candidate signature localizes to cluster 6
# ============================================================================
def make_fig2_candidate():

    import anndata as ad
    import scipy.sparse as sp
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    A = ad.read_h5ad(str(H5AD))
    U = A.obsm["X_umap"]
    lineage = np.asarray(A.obs["lineage"]).astype(str)
    var_names = np.asarray(A.var_names)

    CAND = ["Robo2", "Abi3bp", "Gabrg1", "Adcyap1", "Chrm3", "Rprm", "Thrb", "Cntn5"]
    gi = [int(np.where(var_names == g)[0][0]) for g in CAND]
    Xc = A[:, gi].X
    Xc = Xc.toarray() if sp.issparse(Xc) else np.asarray(Xc)
    zc = (Xc - Xc.mean(0)) / (Xc.std(0) + 1e-9)
    sig = zc.mean(1)

    top3 = ["Adcyap1", "Gabrg1", "Abi3bp"]
    gj = [int(np.where(var_names == g)[0][0]) for g in top3]
    Xt = A[:, gj].X
    Xt = Xt.toarray() if sp.issparse(Xt) else np.asarray(Xt)

    PT = 1.0

    style_umap = _umap_axes

    LX, LY = 8.5, -14.4

    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=[1.55, 1, 1],
        height_ratios=[1, 1],
        wspace=0.14,
        hspace=0.16,
        left=0.045,
        right=0.975,
        top=0.90,
        bottom=0.06,
    )
    axb = fig.add_subplot(gs[:, 0])
    order = np.argsort(sig)
    s = axb.scatter(
        U[order, 0],
        U[order, 1],
        c=sig[order],
        s=PT,
        cmap="magma",
        linewidths=0,
        rasterized=True,
        vmin=np.percentile(sig, 1),
        vmax=np.percentile(sig, 99),
    )
    axb.set_title("8-gene projection-candidate signature")
    axb.text(
        LX,
        LY,
        "cluster 6\n(candidate population)",
        fontsize=11,
        ha="center",
        va="center",
        fontweight="bold",
    )
    axb.annotate(
        "",
        xy=(4.8, -0.6),
        xytext=(LX, LY + 1.5),
        arrowprops=dict(
            arrowstyle="-|>",
            color="black",
            lw=1.6,
            mutation_scale=18,
            connectionstyle="arc3,rad=0.45",
        ),
    )
    style_umap(axb)
    axb.set_anchor("N")
    cax = inset_axes(
        axb,
        width="3.5%",
        height="55%",
        loc="center left",
        bbox_to_anchor=(1.02, 0.0, 1, 1),
        bbox_transform=axb.transAxes,
        borderpad=0,
    )
    cb = fig.colorbar(s, cax=cax)
    cb.set_label("signature score (mean z, 8 genes)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    sm = [
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[0, 2]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[1, 2]),
    ]
    for a in sm:
        a.set_anchor("C")
    for j, g in enumerate(top3):
        ax = sm[j]
        v = Xt[:, j]
        o = np.argsort(v)
        ax.scatter(U[o, 0], U[o, 1], c=v[o], s=PT, cmap="viridis", linewidths=0, rasterized=True)
        ax.set_title(f"$\\it{{{g}}}$")
        style_umap(ax)
    axl = sm[3]
    for labn in ["Inhibitory", "Excitatory", "Unassigned"]:
        mm = lineage == labn
        axl.scatter(
            U[mm, 0],
            U[mm, 1],
            s=PT,
            c=ps.CELLTYPE_COLORS[labn],
            linewidths=0,
            rasterized=True,
            label=labn,
        )
    axl.set_title("lineage")
    style_umap(axl)
    axl.legend(
        markerscale=6,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=3,
        handletextpad=0.2,
        columnspacing=1.0,
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(
        "The projection-candidate signature localizes to one excitatory population (cluster 6)",
        fontsize=14,
        fontweight="bold",
        y=0.975,
    )
    fig.canvas.draw()
    fig.savefig(str(OUT_MAIN / "Fig2_candidate_population.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", "Fig2_candidate_population.png")


# ============================================================================
# Figure 3 — Spatial organization (ABC MERFISH: AP violin, AP scatter, dl/vm bars)
# ============================================================================
def make_fig3_spatial():

    import matplotlib.transforms as mtransforms
    from matplotlib.patches import Patch
    from scipy.stats import kruskal

    coords = pd.read_csv(str(TAB / "aon_merfish_coords.csv"))
    T = str(TAB)
    dl = pd.read_csv(f"{T}/spatial_cluster_dl_assignment.csv")

    DIR = ps.DIRECTION_COLORS
    coords["st"] = coords["supertype"].str.extract(r"(Glut_\d)")
    order = list(coords.groupby("st")["z"].median().sort_values().index)
    st_colors = {st: plt.get_cmap("viridis")(i / (len(order) - 1)) for i, st in enumerate(order)}

    fig = plt.figure(figsize=(15, 5.2))
    gs = fig.add_gridspec(1, 3, wspace=0.30)

    axA = fig.add_subplot(gs[0, 0])
    data = [coords.loc[coords["st"] == st, "z"].values for st in order]
    parts = axA.violinplot(data, showmedians=True, widths=0.8)
    for i, b in enumerate(parts["bodies"]):
        b.set_facecolor(st_colors[order[i]])
        b.set_alpha(0.8)
        b.set_edgecolor("grey")
    for k in ["cmedians", "cbars", "cmins", "cmaxes"]:
        parts[k].set_color("grey")
    N = sum(len(d) for d in data)
    H, _ = kruskal(*data)
    k2 = H / ((N**2 - 1) / (N + 1))
    axA.set_xticks(range(1, len(order) + 1))
    axA.set_xticklabels(order, rotation=30, ha="right", fontsize=9)
    axA.set_ylabel("Z: anterior → posterior (CCF)")
    axA.set_title(
        "A  Anteroposterior position by supertype", loc="left", fontweight="bold", fontsize=12
    )
    axA.text(
        0.97,
        0.04,
        f"Kruskal–Wallis $\\epsilon^2$ = {k2:.2f}",
        transform=axA.transAxes,
        va="bottom",
        ha="right",
        fontsize=9,
        color="0.3",
    )
    for s_ in ["top", "right"]:
        axA.spines[s_].set_visible(False)

    axB = fig.add_subplot(gs[0, 1])
    ZTHR = 11.6
    grp = np.where(coords["z"] >= ZTHR, "posterior", "anterior")
    for g, col in [("anterior", DIR["anterior"]), ("posterior", DIR["posterior"])]:
        m = grp == g
        axB.scatter(
            coords.loc[m, "z"],
            coords.loc[m, "y"],
            s=2,
            color=col,
            linewidths=0,
            rasterized=True,
            label=f"{g} (n={int(m.sum()):,})",
        )
    axB.axvline(ZTHR, ls="--", color="0.4", lw=1)
    axB.invert_yaxis()
    axB.set_xlabel("Z: anterior → posterior (CCF)")
    axB.set_ylabel("Y: dorsal → ventral (CCF, inverted)")
    trans = mtransforms.blended_transform_factory(axB.transData, axB.transAxes)
    axB.text(
        ZTHR + 0.35, 0.97, "Z = 11.6", transform=trans, ha="left", va="top", fontsize=8, color="0.4"
    )
    axB.set_title(
        "B  Anteroposterior is the grouping axis", loc="left", fontweight="bold", fontsize=12
    )
    axB.legend(
        markerscale=4, loc="upper left", frameon=True, framealpha=0.85, edgecolor="none", fontsize=8
    )
    for s_ in ["top", "right"]:
        axB.spines[s_].set_visible(False)

    axC = fig.add_subplot(gs[0, 2])
    dl_sorted = dl.sort_values("med_dl_score")
    ys = np.arange(len(dl_sorted))
    cols = [
        DIR["dorsolateral"] if g == "dorsolateral" else DIR["ventromedial"]
        for g in dl_sorted["group"]
    ]
    axC.barh(ys, dl_sorted["med_dl_score"], color=cols, edgecolor="k", linewidth=0.4)
    axC.axvline(0, color="0.5", lw=0.8)
    axC.set_yticks(ys)
    axC.set_yticklabels(dl_sorted["supertype"], fontsize=7)
    axC.set_xlabel("median dorsolateral score")
    axC.set_title(
        "C  Dorsolateral vs ventromedial clusters", loc="left", fontweight="bold", fontsize=12
    )
    axC.legend(
        handles=[
            Patch(color=DIR["dorsolateral"], label="dorsolateral"),
            Patch(color=DIR["ventromedial"], label="ventromedial"),
        ],
        loc="upper left",
        frameon=True,
        framealpha=0.85,
        edgecolor="none",
        fontsize=8,
    )
    for s_ in ["top", "right"]:
        axC.spines[s_].set_visible(False)

    fig.suptitle(
        "Figure 3  |  Spatial organization of AON glutamatergic populations (ABC MERFISH atlas)",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(str(OUT_MAIN / "Fig3_spatial.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", "Fig3_spatial.png")


# ============================================================================
# Figure 4 — Differential expression (volcano + effect sizes, LFC-shrunk)
# ============================================================================
def make_fig4_de():

    from matplotlib.patches import Patch
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    T = str(TAB)
    counts = pd.read_csv(f"{T}/pseudobulk_counts.csv", index_col=0)
    coldata = pd.read_csv(f"{T}/pseudobulk_coldata.csv", index_col=0)
    coldata = coldata.loc[counts.index]
    coldata["group"] = pd.Categorical(coldata["group"], categories=["ventromedial", "dorsolateral"])
    coldata["donor"] = coldata["donor"].astype("category")

    inf = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts.astype(int),
        metadata=coldata,
        design_factors=["donor", "group"],
        ref_level=["group", "ventromedial"],
        inference=inf,
        quiet=True,
    )
    dds.deseq2()

    st = DeseqStats(
        dds, contrast=["group", "dorsolateral", "ventromedial"], inference=inf, quiet=True
    )
    st.summary()
    res = st.results_df.copy()

    de = pd.read_csv(f"{T}/de_dorsolateral_vs_ventromedial_deseq2.csv")
    id2sym = dict(zip(de["gene_identifier"], de["symbol"]))
    sym2id = {v: k for k, v in id2sym.items()}

    shr = pd.read_csv(f"{T}/de_dorsolateral_vs_ventromedial_deseq2_lfcShrink.csv")
    shr = shr.rename(columns={"Unnamed: 0": "gene_identifier"})
    sym2shr = {r["symbol"]: (r["log2FoldChange"], r["lfcSE"], r["padj"]) for _, r in shr.iterrows()}

    DIR = ps.DIRECTION_COLORS

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.24)

    # ---- Panel A: volcano with leader lines ----
    axA = fig.add_subplot(gs[0, 0])
    x = res["log2FoldChange"].values
    y = -np.log10(res["padj"].replace(0, np.nan)).values
    sg = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)
    up_dl = sg & (res["log2FoldChange"] > 0)
    up_vm = sg & (res["log2FoldChange"] < 0)
    axA.scatter(x[~sg], y[~sg], s=3, c="0.8", linewidths=0, rasterized=True)
    axA.scatter(
        x[up_dl],
        y[up_dl],
        s=4,
        color=DIR["dorsolateral"],
        linewidths=0,
        rasterized=True,
        label=f"up dorsolateral ({int(up_dl.sum())})",
    )
    axA.scatter(
        x[up_vm],
        y[up_vm],
        s=4,
        color=DIR["ventromedial"],
        linewidths=0,
        rasterized=True,
        label=f"up ventromedial ({int(up_vm.sum())})",
    )
    label_pos = {"Gabrg1": (3.4, 66), "Abi3bp": (3.6, 45), "Adcyap1": (3.8, 33)}
    for g, (lx, ly) in label_pos.items():
        gid = sym2id[g]
        px = res.loc[gid, "log2FoldChange"]
        py = -np.log10(res.loc[gid, "padj"])
        axA.scatter(
            [px], [py], s=26, facecolors="none", edgecolors="black", linewidths=1.0, zorder=5
        )
        axA.annotate(
            f"$\\it{{{g}}}$",
            xy=(px, py),
            xytext=(lx, ly),
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="center",
            zorder=6,
            arrowprops=dict(arrowstyle="-", color="0.35", lw=0.7, shrinkA=1, shrinkB=3),
        )
    axA.axvline(1, ls="--", color="0.6", lw=0.8)
    axA.axvline(-1, ls="--", color="0.6", lw=0.8)
    axA.axhline(-np.log10(0.05), ls="--", color="0.6", lw=0.8)
    axA.set_xlabel("log2 fold change (dorsolateral / ventromedial)")
    axA.set_ylabel("−log10 adjusted p")
    axA.set_title(
        "A  Differential expression, dorsolateral vs ventromedial",
        loc="left",
        fontweight="bold",
        fontsize=12,
    )
    axA.legend(markerscale=3, loc="upper left", frameon=False, fontsize=8)
    axA.set_xlim(x[np.isfinite(x)].min() - 0.3, 4.6)
    for s_ in ["top", "right"]:
        axA.spines[s_].set_visible(False)

    # ---- Panel B: candidate bars ----
    axB = fig.add_subplot(gs[0, 1])
    CAND = ["Robo2", "Abi3bp", "Gabrg1", "Adcyap1", "Chrm3", "Rprm", "Thrb", "Cntn5"]
    rows = []
    for g in CAND:
        gid = sym2id[g]
        raw = res.loc[gid, "log2FoldChange"]
        slfc, sse, spadj = sym2shr[g]
        rows.append((g, raw, slfc, sse, spadj))
    cs = (
        pd.DataFrame(rows, columns=["symbol", "lfc_raw", "lfc_shrunk", "lfcSE_shrunk", "padj"])
        .sort_values("lfc_shrunk")
        .reset_index(drop=True)
    )
    for i, r_ in cs.iterrows():
        col = DIR["dorsolateral"] if r_["lfc_shrunk"] > 0 else DIR["ventromedial"]
        filled = r_["padj"] < 0.05
        axB.barh(
            i,
            r_["lfc_shrunk"],
            color=col if filled else "white",
            edgecolor=col,
            linewidth=1.4,
            hatch="" if filled else "///",
            zorder=2,
        )
        axB.errorbar(
            r_["lfc_shrunk"], i, xerr=r_["lfcSE_shrunk"], color="0.3", lw=0.8, capsize=2, zorder=3
        )
        axB.scatter(r_["lfc_raw"], i, marker="|", s=90, color="0.45", zorder=4)
    axB.axvline(0, color="0.5", lw=0.8)
    axB.set_yticks(range(len(cs)))
    axB.set_yticklabels([f"$\\it{{{s}}}$" for s in cs["symbol"]], fontsize=10)
    axB.set_xlabel("shrunken log2 fold change (apeglm)")
    axB.set_title("B  Projection-candidate genes", loc="left", fontweight="bold", fontsize=12)
    axB.set_xlim(0, 3.5)
    axB.legend(
        handles=[
            Patch(facecolor="0.4", label="FDR < 0.05"),
            Patch(facecolor="white", edgecolor="0.4", hatch="///", label="not significant"),
            plt.Line2D([0], [0], marker="|", ls="", color="0.45", label="unshrunken LFC"),
        ],
        loc="lower right",
        frameon=False,
        fontsize=8,
    )
    for s_ in ["top", "right"]:
        axB.spines[s_].set_visible(False)

    fig.suptitle(
        "Figure 4  |  Molecular signature of the dorsolateral AON glutamatergic population",
        fontsize=14,
        fontweight="bold",
        y=1.0,
    )
    fig.savefig(str(OUT_MAIN / "Fig4_DE.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", "Fig4_DE.png")


# ============================================================================
# Figure 5 — Pathway enrichment (dorsal/anterior-upregulated, curated)
# ============================================================================
def make_fig5_enrichment():

    import re

    from matplotlib.lines import Line2D

    T = str(TAB) + "/"
    ea = pd.read_csv(T + "enrichment_anterior_up.csv")
    sig = ea[ea["Adjusted P-value"] < 0.05].copy()
    EXCLUDE = [
        "cancer",
        "carcinoma",
        "tumor",
        "breast",
        "pulmonary",
        "valve",
        "heart",
        "cardiac",
        "endothelial",
    ]
    sig["offtopic"] = sig["Term"].str.lower().apply(lambda t: any(k in t for k in EXCLUDE))
    keep = sig[~sig["offtopic"]].copy()
    lib_short = {
        "GO_Biological_Process_2021": "GO",
        "KEGG_2019_Mouse": "KEGG",
        "Reactome_2022": "Reactome",
    }

    def clean(t):
        return re.sub(r"\s*R-HSA-\d+", "", re.sub(r"\s*\(GO:\d+\)", "", t)).strip()

    keep["label"] = keep["Term"].map(clean) + "  (" + keep["Gene_set"].map(lib_short) + ")"
    keep["overlap_n"] = keep["Overlap"].str.split("/").str[0].astype(int)
    keep["neglog10"] = -np.log10(keep["Adjusted P-value"])
    k = keep.sort_values("neglog10").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    y = np.arange(len(k))
    sc = ax.scatter(
        k["neglog10"],
        y,
        s=k["overlap_n"] * 9,
        c=k["Odds Ratio"],
        cmap="viridis",
        edgecolor="black",
        linewidth=0.6,
        zorder=3,
    )
    xmax = k["neglog10"].max()
    ax.set_xlim(1.2, xmax * 1.12)
    fdr = -np.log10(0.05)
    ax.axvline(fdr, ls="--", color="grey", lw=1, zorder=1)
    ax.text(
        fdr + 0.06,
        len(k) / 2 - 0.5,
        "FDR 0.05",
        fontsize=9,
        color="grey",
        ha="left",
        va="center",
        rotation=90,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(k["label"])
    ax.set_xlabel("-log10 adjusted p-value")
    ax.margins(y=0.10)
    ax.set_title("Pathway enrichment, dorsal/anterior-upregulated genes", fontsize=13)
    cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("odds ratio", fontsize=10)
    cb.ax.tick_params(labelsize=10)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="grey",
            markeredgecolor="black",
            markersize=np.sqrt(n * 9),
            label=f"{n}",
        )
        for n in [10, 30, 46]
    ]
    ax.legend(
        handles=handles,
        title="genes in overlap",
        loc="center left",
        bbox_to_anchor=(1.16, 0.5),
        labelspacing=2.2,
        borderpad=1.0,
        handletextpad=1.2,
    )
    fig.savefig(str(OUT_MAIN / "Fig5_enrichment.png"), dpi=300, bbox_inches="tight")
    print("  wrote", "Fig5_enrichment.png")


# ============================================================================
# Figure 6 — Reference concordance via scANVI label transfer
# ============================================================================
def make_fig6_scanvi():

    from matplotlib.patches import Patch

    pred = pd.read_csv(str(TAB / "scanvi_predictions.csv"))

    gc_ = (
        pred.groupby("leiden", observed=True)
        .agg(
            n=("barcode", "size"),
            lineage=("lineage", lambda s: s.value_counts().index[0]),
            conf=("scanvi_maxprob", "mean"),
            top_call=("scanvi_pred", lambda s: s.value_counts().index[0]),
            top_frac=("scanvi_pred", lambda s: s.value_counts(normalize=True).iloc[0]),
        )
        .reset_index()
    )
    gc_["leiden"] = gc_["leiden"].astype(str)

    CT = ps.CELLTYPE_COLORS
    n_hi = int((gc_["conf"] > 0.9).sum())

    fig = plt.figure(figsize=(15, 6.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.05], wspace=0.32)
    axA = fig.add_subplot(gs[0, 0])
    s = gc_.sort_values("conf").reset_index(drop=True)
    ys = np.arange(len(s))
    cols = [CT.get(lin, "#999999") for lin in s["lineage"]]
    axA.barh(ys, s["conf"], color=cols, edgecolor="k", linewidth=0.4)
    axA.set_yticks(ys)
    axA.set_yticklabels([f"cl {c}" for c in s["leiden"]], fontsize=8)
    axA.axvline(0.8, ls="--", color="0.5", lw=0.9)
    # MOVED: right of the dashed line, small gap, left-aligned
    axA.text(0.815, len(s) - 0.3, "0.8", fontsize=7, color="0.5", va="bottom", ha="left")
    axA.set_xlabel("mean scANVI posterior confidence")
    axA.set_xlim(0, 1.45)
    axA.set_title(
        "A  Per-cluster label-transfer confidence", loc="left", fontweight="bold", fontsize=12
    )
    axA.annotate(
        "leiden-4:\nlowest confidence\n(GABAergic; no matched\nAON reference type)",
        xy=(s.iloc[0]["conf"], 0),
        xytext=(1.02, 2.6),
        fontsize=8,
        color=CT["Inhibitory"],
        fontweight="bold",
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="->", color=CT["Inhibitory"], lw=1.2),
    )
    axA.legend(
        handles=[
            Patch(color=CT["Excitatory"], label="Excitatory"),
            Patch(color=CT["Inhibitory"], label="Inhibitory"),
            Patch(color="#999999", label="Unassigned"),
        ],
        loc="upper right",
        frameon=False,
        fontsize=8,
    )
    for sp_ in ["top", "right"]:
        axA.spines[sp_].set_visible(False)

    axB = fig.add_subplot(gs[0, 1])
    l4 = pred[pred["leiden"].astype(str) == "4"]
    vc = l4["scanvi_pred"].value_counts().head(6)[::-1]

    def nt(x):
        return "Gaba" if "Gaba" in x else ("Glut" if "Glut" in x else "other")

    bcol = [
        CT["Inhibitory"]
        if nt(x) == "Gaba"
        else (CT["Excitatory"] if nt(x) == "Glut" else "#999999")
        for x in vc.index
    ]
    yb = np.arange(len(vc))
    axB.barh(yb, vc.values / len(l4) * 100, color=bcol, edgecolor="k", linewidth=0.4)
    axB.set_yticks(yb)
    axB.set_yticklabels([t.replace(" Gaba", "").replace(" Glut", "") for t in vc.index], fontsize=8)
    axB.set_xlabel("% of leiden-4 nuclei")
    axB.set_title(
        "B  leiden-4 assignment is split and GABAergic", loc="left", fontweight="bold", fontsize=12
    )
    axB.set_xlim(0, 58)
    for i, t in enumerate(vc.index):
        if "Barhl2" in t:
            axB.annotate(
                "← correlation method forced\n   ALL of leiden-4 here (glutamatergic)",
                xy=(vc.values[i] / len(l4) * 100, i),
                xytext=(14, i + 0.15),
                fontsize=7.5,
                color=CT["Excitatory"],
                arrowprops=dict(arrowstyle="->", color=CT["Excitatory"], lw=1),
            )
    axB.legend(
        handles=[
            Patch(color=CT["Inhibitory"], label="GABAergic subclass"),
            Patch(color=CT["Excitatory"], label="glutamatergic subclass"),
        ],
        loc="lower right",
        frameon=False,
        fontsize=8,
    )
    for sp_ in ["top", "right"]:
        axB.spines[sp_].set_visible(False)

    fig.suptitle(
        "Figure 6  |  Probabilistic reference mapping (scANVI) flags the one discordant cluster",
        fontsize=14,
        fontweight="bold",
        y=1.0,
    )
    fig.text(
        0.5,
        -0.02,
        f"scANVI Glut/Gaba calls agree with in-house E/I annotation for 94.4% of 32,615 nuclei; "
        f"{n_hi}/22 clusters map at >0.9 mean confidence.",
        ha="center",
        fontsize=9,
        style="italic",
        color="0.3",
    )
    fig.savefig(str(OUT_MAIN / "Fig6_reference_concordance.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", "Fig6_reference_concordance.png")


# ============================================================================
# Figure S1 — Quality control waterfall and distributions
# ============================================================================
def make_figs1_qc():
    pre = pd.read_csv(str(TAB / "prefilter_qc_metrics.csv"))
    ret = pd.read_csv(str(TAB / "retained_qc_metrics.csv"))
    GREY = "#B0B0B0"
    KEEP = ps.CELLTYPE_COLORS["Excitatory"]

    fig = plt.figure(figsize=(14, 7.5))
    gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.26)

    def dist_panel(ax, pre_v, lo, hi, xlabel, title, letter, logx=False):
        data = pre_v.values
        if logx:
            data = data[data > 0]
            bins = np.logspace(np.log10(data.min()), np.log10(data.max()), 60)
            ax.set_xscale("log")
        else:
            bins = np.linspace(data.min(), data.max(), 60)
        ax.hist(data, bins=bins, color=GREY, edgecolor="none")
        if lo is not None:
            ax.axvline(lo, color="k", ls="--", lw=1.2)
        if hi is not None:
            ax.axvline(hi, color="k", ls="--", lw=1.2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("nuclei")
        ax.set_title(f"{letter}  {title}", loc="left", fontsize=12)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)

    axA = fig.add_subplot(gs[0, 0])
    dist_panel(axA, pre.n_genes_by_counts, 500, 6000, "genes per nucleus", "genes detected", "A")
    xr = axA.get_xlim()[1] - axA.get_xlim()[0]
    off = 0.012 * xr
    axA.text(
        500 + off,
        axA.get_ylim()[1] * 0.95,
        "500",
        fontsize=9,
        ha="left",
        va="top",
        rotation=90,
        color="0.3",
    )
    axA.text(
        6000 + off,
        axA.get_ylim()[1] * 0.95,
        "6000",
        fontsize=9,
        ha="left",
        va="top",
        rotation=90,
        color="0.3",
    )

    axB = fig.add_subplot(gs[0, 1])
    dist_panel(
        axB, pre.total_counts, None, None, "total UMI (log scale)", "library size", "B", logx=True
    )

    axC = fig.add_subplot(gs[0, 2])
    d = np.clip(pre.pct_counts_mt.values, 0, 10)
    axC.hist(d, bins=np.linspace(0, 10, 50), color=GREY, edgecolor="none")
    axC.axvline(2, color="k", ls="--", lw=1.2)
    axC.text(
        2 + 0.15, axC.get_ylim()[1] * 0.9, "2%", fontsize=9, ha="left", rotation=90, color="0.3"
    )
    axC.set_xlabel("% mitochondrial counts (clipped at 10)")
    axC.set_ylabel("nuclei")
    axC.set_title("C  mitochondrial fraction", loc="left", fontsize=12)
    for s in ["top", "right"]:
        axC.spines[s].set_visible(False)

    axD = fig.add_subplot(gs[1, 0])
    axD.scatter(
        pre.total_counts,
        pre.n_genes_by_counts,
        s=2,
        c=GREY,
        linewidths=0,
        rasterized=True,
        label=f"all ({len(pre):,})",
    )
    axD.scatter(
        ret.total_counts,
        ret.n_genes_by_counts,
        s=2,
        c=KEEP,
        linewidths=0,
        rasterized=True,
        label=f"retained ({len(ret):,})",
    )
    axD.set_xscale("log")
    axD.set_xlabel("total UMI (log)")
    axD.set_ylabel("genes per nucleus")
    axD.set_title("D  retained population", loc="left", fontsize=12)
    axD.legend(markerscale=4, loc="lower right")
    for s in ["top", "right"]:
        axD.spines[s].set_visible(False)

    axE = fig.add_subplot(gs[1, 1:])
    steps = ["Cell Ranger\nfiltered", "+ gene & mito\nthresholds", "+ doublet\nremoval"]
    counts = [36863, 35028, 32952]
    bars = axE.bar(steps, counts, color=[GREY, GREY, KEEP], edgecolor="k", linewidth=0.5, width=0.6)
    for b, c in zip(bars, counts):
        axE.text(
            b.get_x() + b.get_width() / 2,
            c + 300,
            f"{c:,}",
            ha="center",
            fontsize=11,
            fontweight="bold",
        )
    axE.annotate(f"−{36863 - 35028:,}", xy=(0.5, 35028), ha="center", fontsize=9, color="0.3")
    axE.annotate(f"−{35028 - 32952:,}", xy=(1.5, 32952), ha="center", fontsize=9, color="0.3")
    axE.set_ylabel("nuclei")
    axE.set_ylim(0, 40000)
    axE.set_title("E  nuclei retained through QC", loc="left", fontsize=12)
    for s in ["top", "right"]:
        axE.spines[s].set_visible(False)

    fig.suptitle(
        "Quality control: 36,863 nuclei → 32,952 after filtering and doublet removal",
        fontsize=14,
        fontweight="bold",
        y=0.99,
    )
    fig.savefig(str(OUT_SUPPLEMENT / "FigS1_QC.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", "FigS1_QC.png")


# ============================================================================
# Figure S2 — DE robustness along the anteroposterior axis
# ============================================================================
def make_figs2_de_ap():

    from adjustText import adjust_text

    T = str(TAB)
    ap = pd.read_csv(f"{T}/de_posterior_vs_anterior_deseq2.csv")
    DIR = ps.DIRECTION_COLORS
    x = ap["log2FoldChange"].values
    y = -np.log10(ap["padj"].replace(0, np.nan)).values
    sig = (ap["padj"] < 0.05) & (ap["log2FoldChange"].abs() > 1)
    up_post = sig & (ap["log2FoldChange"] > 0)
    up_ant = sig & (ap["log2FoldChange"] < 0)

    fig, axA = plt.subplots(figsize=(7.5, 6))
    axA.scatter(x[~sig], y[~sig], s=3, c="0.8", linewidths=0, rasterized=True)
    axA.scatter(
        x[up_post],
        y[up_post],
        s=4,
        color=DIR["posterior"],
        linewidths=0,
        rasterized=True,
        label=f"up posterior ({int(up_post.sum())})",
    )
    axA.scatter(
        x[up_ant],
        y[up_ant],
        s=4,
        color=DIR["anterior"],
        linewidths=0,
        rasterized=True,
        label=f"up anterior ({int(up_ant.sum())})",
    )

    lab = ap.dropna(subset=["padj"]).sort_values("padj").head(7).copy()
    lab["yy"] = lab["padj"].apply(lambda p: -np.log10(p) if p > 0 else 300)
    texts = []
    for _, row in lab.iterrows():
        px, py = row["log2FoldChange"], row["yy"]
        axA.scatter(
            [px], [py], s=18, facecolors="none", edgecolors="black", linewidths=0.7, zorder=5
        )
        texts.append(
            axA.text(
                px, py, f"$\\it{{{row['symbol']}}}$", fontsize=8.5, fontweight="bold", zorder=6
            )
        )
    adjust_text(
        texts,
        ax=axA,
        expand=(1.4, 1.8),
        arrowprops=dict(arrowstyle="-", color="0.45", lw=0.6),
        force_text=(0.5, 0.8),
        only_move={"text": "xy", "static": "xy"},
    )
    axA.axvline(1, ls="--", color="0.6", lw=0.8)
    axA.axvline(-1, ls="--", color="0.6", lw=0.8)
    axA.axhline(-np.log10(0.05), ls="--", color="0.6", lw=0.8)
    axA.set_xlim(x[np.isfinite(x)].min() - 0.5, x[np.isfinite(x)].max() + 0.5)
    axA.set_xlabel("log2 fold change (posterior / anterior)")
    axA.set_ylabel("−log10 adjusted p")
    axA.set_title(
        "Figure S2  |  DE robustness: anteroposterior axis",
        loc="left",
        fontweight="bold",
        fontsize=13,
    )
    axA.legend(markerscale=3, loc="lower right", frameon=False, fontsize=9)
    for s_ in ["top", "right"]:
        axA.spines[s_].set_visible(False)
    fig.savefig(
        str(OUT_SUPPLEMENT / "FigS2_DE_robustness_AP_axis.png"), dpi=300, bbox_inches="tight"
    )
    plt.close(fig)
    print("  wrote", "FigS2_DE_robustness_AP_axis.png")


# ============================================================================
# Runner
# ============================================================================
FIGURES = {
    "Fig1_atlas": make_fig1_atlas,
    "umap_leiden_standalone": make_umap_leiden_standalone,
    "umap_celltypes_standalone": make_umap_celltypes_standalone,
    "Fig2_candidate": make_fig2_candidate,
    "Fig3_spatial": make_fig3_spatial,
    "Fig4_DE": make_fig4_de,
    "Fig5_enrichment": make_fig5_enrichment,
    "Fig6_scanvi": make_fig6_scanvi,
    "FigS1_QC": make_figs1_qc,
    "FigS2_DE_AP": make_figs2_de_ap,
}


def main(argv):
    targets = argv[1:] if len(argv) > 1 else list(FIGURES)
    unknown = [t for t in targets if t not in FIGURES]
    if unknown:
        sys.exit(f"unknown figure(s): {unknown}\nchoose from: {list(FIGURES)}")
    for t in targets:
        print(f"[{t}] rendering ...")
        FIGURES[t]()
    print(f"done — {len(targets)} figure(s) written to results/figures_final/")


if __name__ == "__main__":
    main(sys.argv)
