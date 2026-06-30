#!/usr/bin/env python3
"""
Pseudobulk expression heatmaps underlying the z-score DEG definition (script 07b /
iteration14 script 18), broken out by the FOUR duplication categories (script 06).

The z-score DEG call is built on a per-gene PSEUDOBULK matrix: mean log1p(CP10K) per
major_lineage PBMC cell type (one value per cell type), then row z-scored across cell
types; DEG = z > Z_THR in >= 1 cell type. This script renders that matrix as heatmaps so
the expression patterns behind the OR can be seen, grouped by dup_class4.

Universe = cleaned 4-group catalog (results/TF_dup_2x2_classification.tsv) intersected with
PBMC-detectable genes (human_pbmc.h5ad var_names) -- identical to scripts 07/07b.

Z_THR is configurable: `python3 07c_...py [Z_THR]` (default 2.0, the cutoff chosen for the
z>2 OR). Outputs are tagged with the threshold.

Outputs (TAG = zscore_z{Z_THR}):
  results/TF_dup4_pseudobulk_by_category.{TAG}.tsv          per-gene raw + z pseudobulk, dup_class4
  results/TF_dup4_pseudobulk_category_means.{TAG}.tsv       category x cell-type mean-z + DEG-rate
  results/figures/TF_dup4_pseudobulk_summary.{TAG}.{pdf,png}    compact category x cell-type heatmaps
  results/figures/TF_dup4_pseudobulk_heatmap.<category>.{TAG}.{pdf,png}   per-category gene-level z heatmap
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/07c_pseudobulk_heatmaps_by_category.py [Z_THR]
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
H5AD = RES / "intermediate" / "h5ad" / "human_pbmc.h5ad"

# Configurable: `python3 07c_...py [Z_THR] [EXPR_THR]`.
#   EXPR_THR (default 0.0): drop TFs whose TPM-like total (SUM over cell types of mean LINEAR
#   CP10K) is <= EXPR_THR BEFORE z-scoring (iteration14 script 19 low-expression pre-filter).
Z_THR = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
EXPR_THR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
TAG = f"zscore_z{Z_THR:g}" + (f"_e{EXPR_THR:g}" if EXPR_THR > 0 else "")
OUT_PER = RES / f"TF_dup4_pseudobulk_by_category.{TAG}.tsv"
OUT_MEAN = RES / f"TF_dup4_pseudobulk_category_means.{TAG}.tsv"
FIG_SUM = RES / "figures" / f"TF_dup4_pseudobulk_summary.{TAG}"
FIG_HM = RES / "figures" / f"TF_dup4_pseudobulk_heatmap"   # + .<category>.<TAG>

CT_COL = "major_lineage"
CT_ORDER = ["progenitor", "monocyte", "cDC1", "cDC2", "DC_other", "pDC",
            "B", "plasma", "CD4 T", "CD8 T", "other T", "NK", "ILC",
            "platelet", "erythrocyte"]
GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
GLABEL = {"dispersed_ohnolog": "dispersed ohnolog", "clustered_cis_ohnolog": "clustered cis-ohnolog",
          "clustered_SSD_tandem": "clustered SSD-tandem", "dispersed_SSD": "dispersed SSD"}
LABEL_MAX = 130   # show gene labels only for categories with <= this many genes


def log(*a):
    print("[07c]", *a, flush=True)


def pseudobulk(genes):
    """Return (cts, expr, expr_total): expr = (n_genes x n_celltypes) mean log1p(CP10K);
    expr_total = per-gene TPM-like SUM over cell types of mean LINEAR CP10K (script 19)."""
    log(f"reading {H5AD} (backed) ...")
    A = ad.read_h5ad(H5AD, backed="r")
    vn = pd.Index(A.var_names)
    ct = A.obs[CT_COL].astype(str).values
    cts = [c for c in CT_ORDER if c in set(ct)] + [c for c in pd.unique(ct) if c not in CT_ORDER]

    idx = vn.get_indexer(genes)          # all >=0 since universe = catalog ∩ var_names
    X = A[:, idx].to_memory().X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)

    log("computing per-cell library size (chunked) ...")
    n = A.shape[0]
    lib = np.zeros(n)
    step = 20000
    for i in range(0, n, step):
        chunk = A[i:i + step].to_memory().X
        chunk = chunk if sp.issparse(chunk) else sp.csr_matrix(chunk)
        lib[i:i + step] = np.asarray(chunk.sum(1)).ravel()
    lib[lib == 0] = 1.0

    expr = np.zeros((len(genes), len(cts)))      # mean log1p(CP10K)
    lin = np.zeros((len(genes), len(cts)))       # mean linear CP10K
    for j, c in enumerate(cts):
        m = ct == c
        if int(m.sum()) == 0:
            continue
        cp = X[m].multiply(1e4 / lib[m][:, None]).tocsr()
        lin[:, j] = np.asarray(cp.mean(0)).ravel()
        cp.data = np.log1p(cp.data)
        expr[:, j] = np.asarray(cp.mean(0)).ravel()
    return cts, expr, lin.sum(1)


def main():
    cls = pd.read_csv(CLASS, sep="\t")
    det = set(ad.read_h5ad(H5AD, backed="r").var_names)
    u = cls[cls.gene.isin(det)].reset_index(drop=True)
    log(f"4-group catalog: {len(cls)}; PBMC-detectable: {len(u)}")

    genes = u.gene.values
    cts, expr, expr_total = pseudobulk(genes)

    # ---- low-expression pre-filter (TPM-like total of linear CP10K) ----
    u["expr_total_cp10k"] = expr_total
    if EXPR_THR > 0:
        keep = expr_total > EXPR_THR
        log(f"low-expression filter EXPR_THR={EXPR_THR:g}: {int(keep.sum())}/{len(u)} TFs "
            f"survive, {int((~keep).sum())} dropped")
        u = u[keep].reset_index(drop=True)
        expr, expr_total, genes = expr[keep], expr_total[keep], genes[keep]
    log(f"universe after expression filter: {len(u)}")

    mu = expr.mean(1, keepdims=True)
    sd = expr.std(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = np.where(sd > 0, (expr - mu) / sd, 0.0)
    deg = Z > Z_THR

    # ---- per-gene table: raw pseudobulk + z + DEG label ----
    per = u[["gene", "dup_class4", "provenance", "arrangement", "TF_subfamily", "cluster_id",
             "expr_total_cp10k"]].copy()
    for j, c in enumerate(cts):
        per[f"pb_{c}"] = expr[:, j]
    for j, c in enumerate(cts):
        per[f"z_{c}"] = Z[:, j]
    per["n_DEG_celltypes"] = deg.sum(1)
    per["is_DEG"] = per["n_DEG_celltypes"] > 0
    per["peak_celltype"] = [cts[k] for k in Z.argmax(1)]
    per.to_csv(OUT_PER, sep="\t", index=False)
    log(f"DEG (z>{Z_THR:g}) overall: {int(per.is_DEG.sum())}/{len(per)} = {per.is_DEG.mean():.1%}")

    # ---- category x cell-type means: mean z + DEG rate ----
    mrows = []
    for g in GROUPS:
        m = (u.dup_class4 == g).values
        row = dict(dup_class4=g, n=int(m.sum()))
        for j, c in enumerate(cts):
            row[f"meanz_{c}"] = float(Z[m, j].mean())
            row[f"degrate_{c}"] = float(deg[m, j].mean())
        mrows.append(row)
    means = pd.DataFrame(mrows)
    means.to_csv(OUT_MEAN, sep="\t", index=False)

    make_summary(means, cts)
    for g in GROUPS:
        make_category_heatmap(g, u, Z, deg, cts)
    log("done.")


def make_summary(means, cts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    meanz = means[[f"meanz_{c}" for c in cts]].values
    degr = means[[f"degrate_{c}" for c in cts]].values
    ylab = [f"{GLABEL[g]}\n(n={int(means.loc[means.dup_class4==g,'n'].iloc[0])})" for g in GROUPS]

    fig, axes = plt.subplots(2, 1, figsize=(11, 6.4))

    ax = axes[0]
    vmax = np.nanmax(np.abs(meanz))
    im = ax.imshow(meanz, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_title(f"(A) Mean pseudobulk z-score per category x PBMC cell type "
                 f"(row z of mean log1p CP10K)", fontsize=10, loc="left")
    for i in range(len(GROUPS)):
        for j in range(len(cts)):
            ax.text(j, i, f"{meanz[i,j]:.1f}", ha="center", va="center", fontsize=6.5,
                    color="white" if abs(meanz[i, j]) > vmax * 0.6 else "black")
    ax.set_yticks(range(len(GROUPS))); ax.set_yticklabels(ylab, fontsize=8)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels([], fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01).set_label("mean z", fontsize=8)

    ax = axes[1]
    im = ax.imshow(degr, aspect="auto", cmap="Reds", vmin=0, vmax=np.nanmax(degr) * 1.05)
    ax.set_title(f"(B) Fraction of category that is a z-score DEG (z>{Z_THR:g}) per cell type",
                 fontsize=10, loc="left")
    for i in range(len(GROUPS)):
        for j in range(len(cts)):
            ax.text(j, i, f"{degr[i,j]:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white" if degr[i, j] > np.nanmax(degr) * 0.6 else "black")
    ax.set_yticks(range(len(GROUPS))); ax.set_yticklabels(ylab, fontsize=8)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=90, fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01).set_label("DEG rate", fontsize=8)

    efilt = f"; expr>{EXPR_THR:g} CP10K-sum filter, n={int(means.n.sum())}" if EXPR_THR > 0 else ""
    fig.suptitle("Pseudobulk expression behind the z-score DEG, by duplication category "
                 f"(human PBMC; Lambert TFs, KRAB/HOX/readthrough/pseudo removed{efilt})",
                 fontsize=11.5, y=1.0)
    fig.tight_layout()
    FIG_SUM.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(FIG_SUM) + ".pdf", bbox_inches="tight")
    fig.savefig(str(FIG_SUM) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {FIG_SUM}.pdf/.png")


def make_category_heatmap(g, u, Z, deg, cts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m = (u.dup_class4 == g).values
    Zm, dm, genes = Z[m], deg[m], u.gene.values[m]
    order = np.lexsort((-Zm.max(1), Zm.argmax(1)))   # group by peak cell type
    Zm, dm, genes = Zm[order], dm[order], genes[order]
    n = len(genes)
    n_deg = int(dm.any(1).sum())

    h = max(2.6, min(0.16 * n, 26))
    fig, ax = plt.subplots(figsize=(7.2, h))
    im = ax.imshow(Zm, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
    ys, xs = np.where(dm)
    ax.scatter(xs, ys, s=6, facecolors="none", edgecolors="black", linewidths=0.35, marker="s")
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=90, fontsize=8)
    if n <= LABEL_MAX:
        ax.set_yticks(range(n)); ax.set_yticklabels(genes, fontsize=max(3.0, min(6.0, 600 / n)))
    else:
        ax.set_yticks([]); ax.set_ylabel(f"{n} TFs (sorted by peak cell type)", fontsize=9)
    ax.set_title(f"{GLABEL[g]}: pseudobulk z-score across PBMC cell types\n"
                 f"n={n} TFs, {n_deg} are z-score DEG (z>{Z_THR:g}); black squares = DEG cells",
                 fontsize=10, loc="left")
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    cax = fig.add_axes([0.93, 0.4, 0.014, 0.2])
    fig.colorbar(im, cax=cax).set_label("row z-score", fontsize=8)
    out = f"{FIG_HM}.{g}.{TAG}"
    fig.savefig(out + ".pdf", bbox_inches="tight")
    fig.savefig(out + ".png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    log(f"wrote {out}.pdf/.png  (n={n}, DEG={n_deg})")


if __name__ == "__main__":
    main()
