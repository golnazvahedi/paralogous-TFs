#!/usr/bin/env python3
"""
Cytokine-dictionary DEG heatmaps broken out by the FOUR duplication categories (script 06).

Visual style follows iteration14 script 09 (clustered_TF_cytokine_heatmap): per category, a
TWO-PANEL heatmap of the top cytokine-responsive TFs --
  (A) TF x cell type, colour = # distinct cytokines eliciting a response (response breadth);
  (B) TF x cytokine,  colour = # distinct cell types responding (top cytokines);
rows ranked by total response events, annotated with a TF-subfamily colour strip + discrete
integer viridis colormap (0 -> light grey) with the count printed in each non-zero cell.
A compact 4-category x cell-type summary heatmap is emitted alongside.

Source DEGs: results/intermediate/human_cytokine_dict_DEGs.tsv (Oesinghaus 2025 pseudobulk,
gated |Avg_log2FC|>=1 & FDR<0.05 in script 09). Universe = cleaned 4-group catalog
(results/TF_dup_2x2_classification.tsv) intersected with PBMC-detectable genes
(human_pbmc.h5ad var_names, backed) -- same universe rule as scripts 07/10.

Outputs (TAG = cytokine):
  results/figures/TF_dup4_cytokine_DEG_summary.{TAG}.{pdf,png}
        4 categories x cell types: responder rate + mean response breadth
  results/figures/TF_dup4_cytokine_DEG_heatmap.<category>.{TAG}.{pdf,png}
        per-category two-panel breadth heatmap (iteration14 style)
  results/TF_dup4_cytokine_DEG_breadth.{TAG}.tsv          per-TF x cell-type breadth + dup_class4
  results/TF_dup4_cytokine_DEG_category_means.{TAG}.tsv   category x cell-type rate + mean breadth
  results/TF_dup4_cytokine_DEG_ranking.{TAG}.tsv          per-TF response-event ranking + category

Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/10b_cytokine_DEG_heatmaps_by_category.py [TOP_N]
  TOP_N (default 30): max TFs shown per category (ranked by total response events).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
H5AD = RES / "intermediate" / "h5ad" / "human_pbmc.h5ad"
DEGS = RES / "intermediate" / "human_cytokine_dict_DEGs.tsv"

TOP_N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
TOP_CYT = 22
TAG = "cytokine"
OUT_BREADTH = RES / f"TF_dup4_cytokine_DEG_breadth.{TAG}.tsv"
OUT_MEAN = RES / f"TF_dup4_cytokine_DEG_category_means.{TAG}.tsv"
OUT_RANK = RES / f"TF_dup4_cytokine_DEG_ranking.{TAG}.tsv"
FIG_SUM = RES / "figures" / f"TF_dup4_cytokine_DEG_summary.{TAG}"
FIG_HM = RES / "figures" / "TF_dup4_cytokine_DEG_heatmap"   # + .<category>.<TAG>

CT_ORDER = ["CD14_Mono", "CD16_Mono", "Mono", "cDC", "pDC", "Granulocyte",
            "B_cell", "Naive_B_cell", "Intermediate_B_cell", "Plasmablast",
            "CD4_T_cell", "CD4_Naive_T_cell", "CD4_Memory_T_cell", "Treg",
            "CD8_T_cell", "CD8_Naive_T_cell", "CD8_Memory_T_cell", "MAIT", "NKT",
            "NK", "NK_CD56hi", "NK_CD56low", "ILC", "HSPC"]
GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
GLABEL = {"dispersed_ohnolog": "dispersed ohnolog", "clustered_cis_ohnolog": "clustered cis-ohnolog",
          "clustered_SSD_tandem": "clustered SSD-tandem", "dispersed_SSD": "dispersed SSD"}


def log(*a):
    print("[10b]", *a, flush=True)


def main():
    cls = pd.read_csv(CLASS, sep="\t")
    log(f"reading PBMC var_names (backed) from {H5AD} ...")
    det = set(ad.read_h5ad(H5AD, backed="r").var_names)
    u = cls[cls.gene.isin(det)].reset_index(drop=True)
    log(f"cleaned 4-group catalog: {len(cls)}; PBMC-detectable: {len(u)}")
    gene_cat = dict(zip(u.gene, u.dup_class4))
    gene_fam = dict(zip(u.gene, u.TF_subfamily.fillna("other")))

    d = pd.read_csv(DEGS, sep="\t")
    d = d[d.human_symbol.isin(set(u.gene))].copy()
    cts = [c for c in CT_ORDER if c in set(d["sheet"])] + \
          [c for c in sorted(d["sheet"].unique()) if c not in CT_ORDER]
    log(f"cytokine DEG events for catalog TFs: {len(d)} across {len(cts)} cell types, "
        f"{d['Cytokine'].nunique()} cytokines")

    # ---- per-TF response-event ranking (for row ordering + table) ----
    rank = (d.groupby("human_symbol")
              .agg(n_events=("Cytokine", "size"),
                   n_celltypes=("sheet", "nunique"),
                   n_cytokines=("Cytokine", "nunique"),
                   n_up=("direction", lambda s: int((s == "up").sum())),
                   n_down=("direction", lambda s: int((s == "down").sum())))
              .reset_index().rename(columns={"human_symbol": "gene"}))
    rank["dup_class4"] = rank.gene.map(gene_cat)
    rank["TF_subfamily"] = rank.gene.map(gene_fam)
    rank = rank.sort_values("n_events", ascending=False).reset_index(drop=True)
    rank.to_csv(OUT_RANK, sep="\t", index=False)

    # ---- per-TF x cell-type breadth matrix (all catalog TFs, for table + summary) ----
    genes = list(u.gene.values)
    gi = {g: k for k, g in enumerate(genes)}
    ci = {c: k for k, c in enumerate(cts)}
    breadth = np.zeros((len(genes), len(cts)), int)
    for (g, c), sub in d.groupby(["human_symbol", "sheet"]):
        breadth[gi[g], ci[c]] = sub["Cytokine"].nunique()
    bt = u[["gene", "dup_class4", "provenance", "arrangement", "TF_subfamily", "cluster_id"]].copy()
    for j, c in enumerate(cts):
        bt[f"breadth_{c}"] = breadth[:, j]
    bt["total_events"] = breadth.sum(1)
    bt["n_response_celltypes"] = (breadth > 0).sum(1)
    bt["is_responder"] = bt["total_events"] > 0
    bt.to_csv(OUT_BREADTH, sep="\t", index=False)
    log(f"responders: {int(bt.is_responder.sum())}/{len(bt)} = {bt.is_responder.mean():.1%}")

    # ---- category x cell-type means ----
    mrows = []
    for g in GROUPS:
        m = (u.dup_class4 == g).values
        row = dict(dup_class4=g, n=int(m.sum()))
        for j, c in enumerate(cts):
            row[f"rate_{c}"] = float((breadth[m, j] > 0).mean())
            row[f"meanbreadth_{c}"] = float(breadth[m, j].mean())
        mrows.append(row)
    means = pd.DataFrame(mrows)
    means.to_csv(OUT_MEAN, sep="\t", index=False)

    make_summary(means, cts)
    for g in GROUPS:
        make_category_panel(g, d, rank, gene_fam, cts)
    log("done.")


# --------------------------------------------------------------------------- summary
def make_summary(means, cts):
    rate = means[[f"rate_{c}" for c in cts]].values
    brd = means[[f"meanbreadth_{c}" for c in cts]].values
    ylab = [f"{GLABEL[g]}\n(n={int(means.loc[means.dup_class4==g,'n'].iloc[0])})" for g in GROUPS]

    fig, axes = plt.subplots(2, 1, figsize=(12, 6.4))
    ax = axes[0]
    im = ax.imshow(rate, aspect="auto", cmap="Purples", vmin=0, vmax=np.nanmax(rate) * 1.05)
    ax.set_title("(A) Fraction of category that is a cytokine-response DEG, per cell type",
                 fontsize=10, loc="left")
    for i in range(len(GROUPS)):
        for j in range(len(cts)):
            ax.text(j, i, f"{rate[i,j]:.2f}", ha="center", va="center", fontsize=6.3,
                    color="white" if rate[i, j] > np.nanmax(rate) * 0.6 else "black")
    ax.set_yticks(range(len(GROUPS))); ax.set_yticklabels(ylab, fontsize=8)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels([], fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01).set_label("responder rate", fontsize=8)

    ax = axes[1]
    im = ax.imshow(brd, aspect="auto", cmap="Reds", vmin=0, vmax=np.nanmax(brd) * 1.05)
    ax.set_title("(B) Mean cytokine-response breadth (# cytokines / TF) per cell type",
                 fontsize=10, loc="left")
    for i in range(len(GROUPS)):
        for j in range(len(cts)):
            ax.text(j, i, f"{brd[i,j]:.1f}", ha="center", va="center", fontsize=6.3,
                    color="white" if brd[i, j] > np.nanmax(brd) * 0.6 else "black")
    ax.set_yticks(range(len(GROUPS))); ax.set_yticklabels(ylab, fontsize=8)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=90, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01).set_label("mean breadth", fontsize=8)

    fig.suptitle("Cytokine-dictionary DEG response by duplication category "
                 "(Oesinghaus 2025; human PBMC; Lambert TFs, KRAB/HOX/readthrough/pseudo removed)",
                 fontsize=11.5, y=1.0)
    fig.tight_layout()
    FIG_SUM.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(FIG_SUM) + ".pdf", bbox_inches="tight")
    fig.savefig(str(FIG_SUM) + ".png", dpi=190, bbox_inches="tight")
    plt.close(fig)
    log(f"wrote {FIG_SUM}.pdf/.png")


# --------------------------------------------------------------- per-category 2-panel
def _fam_colors(tfs, gene_fam):
    fams = [gene_fam.get(t, "other") for t in tfs]
    uniq = list(dict.fromkeys(fams))
    palette = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(uniq), 1)))
    cmap = {f: palette[i] for i, f in enumerate(uniq)}
    return [cmap[f] for f in fams], cmap


def _draw(ax, mat, title, cbar_label, fam_row_colors, show_ylab=True):
    M = mat.values.astype(float)
    vmax = max(1, int(M.max()))
    base = plt.get_cmap("viridis")
    colors = base(np.linspace(0.12, 1.0, vmax))
    colors = np.vstack([[0.94, 0.94, 0.94, 1.0], colors])      # 0 -> light grey
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, vmax + 1.5, 1), cmap.N)
    im = ax.imshow(M, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index if show_ylab else [], fontsize=8)
    ax.set_title(title, fontsize=9)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = int(M[i, j])
            if v:
                ax.text(j, i, v, ha="center", va="center", fontsize=6,
                        color="white" if v > vmax * 0.45 else "black")
    if show_ylab:
        for i, c in enumerate(fam_row_colors):
            ax.add_patch(plt.Rectangle((-0.9, i - 0.5), 0.35, 1.0, color=c,
                                       clip_on=False, transform=ax.transData))
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.025, pad=0.01, ticks=range(0, vmax + 1))
    cb.set_label(cbar_label, fontsize=7); cb.ax.tick_params(labelsize=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def make_category_panel(g, d, rank, gene_fam, cts):
    sub_rank = rank[rank.dup_class4 == g]
    tfs = sub_rank.gene.head(TOP_N).tolist()
    n_resp = int((rank.dup_class4 == g).sum())
    if not tfs:
        log(f"{g}: no responders, skipping")
        return
    dd = d[d.human_symbol.isin(tfs)]

    # matrix A: TF x cell type (# distinct cytokines), cols ordered by total signal
    matA = (dd.groupby(["human_symbol", "sheet"])["Cytokine"].nunique()
              .unstack(fill_value=0).reindex(index=tfs).fillna(0))
    matA = matA[[c for c in cts if c in matA.columns]]
    matA = matA[matA.sum().sort_values(ascending=False).index]

    # matrix B: TF x cytokine (# distinct cell types), top cytokines
    matB = (dd.groupby(["human_symbol", "Cytokine"])["sheet"].nunique()
              .unstack(fill_value=0).reindex(index=tfs).fillna(0))
    matB = matB[matB.sum().sort_values(ascending=False).index[:TOP_CYT]]

    fam_row_colors, fam_cmap = _fam_colors(tfs, gene_fam)
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(15.5, max(3.2, 0.32 * len(tfs) + 2.2)),
        gridspec_kw={"width_ratios": [matA.shape[1], max(matB.shape[1], 1)]})

    _draw(axA, matA, "(A) Cytokine breadth per cell type\n(colour = # cytokines)",
          "# cytokines", fam_row_colors, show_ylab=True)
    axA.set_xlabel("cell type", fontsize=8)
    axA.set_ylabel(f"{GLABEL[g]} TF  (ranked by total response events)", fontsize=8)
    _draw(axB, matB, f"(B) Cell-type breadth per cytokine\n(colour = # cell types; top {matB.shape[1]})",
          "# cell types", fam_row_colors, show_ylab=False)
    axB.set_xlabel("cytokine", fontsize=8); axB.set_ylabel("")

    handles = [plt.Line2D([0], [0], marker="s", ls="", ms=8, color=c, label=f)
               for f, c in fam_cmap.items()]
    fig.legend(handles=handles, title="TF subfamily", fontsize=7, title_fontsize=8,
               loc="lower center", ncol=min(len(handles), 7),
               bbox_to_anchor=(0.5, -0.04), frameon=False)

    shown = f"top {len(tfs)} of {n_resp}" if n_resp > len(tfs) else f"all {n_resp}"
    fig.suptitle(f"Cytokine-responsive {GLABEL[g]} paralogous TFs "
                 f"(human cytokine dictionary, Oesinghaus 2025; KRAB/HOX/readthrough/pseudo removed)\n"
                 f"{n_resp} responders in category ({shown} shown), ranked by total response events",
                 fontsize=11, y=1.0)
    fig.tight_layout(rect=(0.01, 0.03, 1, 0.96))
    out = f"{FIG_HM}.{g}.{TAG}"
    fig.savefig(out + ".pdf", bbox_inches="tight")
    fig.savefig(out + ".png", dpi=190, bbox_inches="tight")
    plt.close(fig)
    log(f"wrote {out}.pdf/.png  (responders={n_resp}, shown={len(tfs)})")


if __name__ == "__main__":
    main()
