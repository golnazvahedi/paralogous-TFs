#!/usr/bin/env python3
"""
PBMC DEG-enrichment OR across the FOUR duplication categories (script 06) -- VARIANT
that defines DEG by the iteration14 Z-SCORE rule (script 18) instead of the canonical
one-vs-rest scRNA marker used by script 07.

WHY THIS EXISTS
---------------
Script 07 (and its iteration14 predecessor 28) define the DEG/identity set as the
one-vs-rest Wilcoxon PBMC lineage markers (results/intermediate/deg/
human_markers.major_lineage.top.tsv). iteration14 deliberately rejected the z-score DEG
set for the OR ("it calls ~every expressed TF a DEG"). This script honours an explicit
request to redo the PBMC OR with the z-score DEG definition instead, so the two DEG
definitions can be compared head-to-head on identical machinery.

Z-SCORE DEG DEFINITION (verbatim from iteration14 script 18_zscore_DEG_TF_by_celltype.py)
----------------------------------------------------------------------------------------
For every gene in the universe, compute its per-cell-type mean expression
(log1p(CP10K), one value per major_lineage cell type), then z-score that gene's profile
ACROSS cell types. A gene is a DEG in cell type c iff its z-score in c is > Z_THR (=1.0).
A gene is a "DEG TF" (is_DEG) iff it is a DEG in >= 1 cell type. Genes with no
cross-cell-type variance (z undefined) are never DEGs.

UNIVERSE (unchanged from script 07)
-----------------------------------
The cleaned 4-group catalog (results/TF_dup_2x2_classification.tsv) intersected with
PBMC-detectable genes read directly from human_pbmc.h5ad var_names. The z-score is then
computed for exactly those universe genes, so the universe rule matches script 07 byte
for byte; only the DEG label differs.

STATS (identical to script 07)
------------------------------
  (1) per-category DEG rate + ONE-VS-REST Fisher OR (+95% CI, p, BH-FDR);
  (2) the key PAIRWISE contrasts (provenance within arrangement, arrangement within
      provenance, and the two marginals);
  (3) provenance x arrangement FACTORIAL logistic regression (adjusted ORs + interaction);
  (4) per-cell-type DEG rate by category (for the figure).

OUTPUTS (".zscore" suffix -- canonical script-07 outputs are NOT touched)
  results/TF_dup4_DEG_OR.zscore.summary.tsv
  results/TF_dup4_DEG_OR.zscore.by_celltype.tsv
  results/TF_dup4_DEG_OR.zscore.per_TF.tsv
  results/figures/TF_dup4_DEG_OR.zscore.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/07b_human_pbmc_dup4_DEG_OR_zscore.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
H5AD = RES / "intermediate" / "h5ad" / "human_pbmc.h5ad"   # detectable universe + expression

# ---- z-score DEG parameters (verbatim from iteration14 scripts 18/19) ----
# Configurable: `python3 07b_...py [Z_THR] [EXPR_THR]`.
#   Z_THR    (default 1.0) : z-score DEG cutoff (DEG iff z > Z_THR in >=1 cell type).
#   EXPR_THR (default 0.0) : TPM-like low-expression pre-filter. A gene's expression total =
#       SUM over cell types of its mean LINEAR CP10K (pseudobulk, iteration14 script 19).
#       Genes with total <= EXPR_THR are DROPPED from the universe BEFORE z-scoring, so a TF
#       that is silent everywhere cannot reach z>Z_THR from noise. EXPR_THR=0 => no filter
#       (z-score normalises away absolute level, so silent TFs otherwise survive).
# Outputs are tagged with the thresholds so cutoffs coexist without clobbering.
Z_THR = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
EXPR_THR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
TAG = f"zscore_z{Z_THR:g}" + (f"_e{EXPR_THR:g}" if EXPR_THR > 0 else "")
OUT_SUM = RES / f"TF_dup4_DEG_OR.{TAG}.summary.tsv"
OUT_CT = RES / f"TF_dup4_DEG_OR.{TAG}.by_celltype.tsv"
OUT_PER = RES / f"TF_dup4_DEG_OR.{TAG}.per_TF.tsv"
OUT_FIG = RES / "figures" / f"TF_dup4_DEG_OR.{TAG}"
CT_COL = "major_lineage"
CT_ORDER = ["progenitor", "monocyte", "cDC1", "cDC2", "DC_other", "pDC",
            "B", "plasma", "CD4 T", "CD8 T", "other T", "NK", "ILC",
            "platelet", "erythrocyte"]

GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
COLORS = {"dispersed_ohnolog": "#2166ac", "clustered_cis_ohnolog": "#5aae61",
          "clustered_SSD_tandem": "#762a83", "dispersed_SSD": "#bababa"}


def log(*a):
    print("[07b]", *a, flush=True)


def or_ci(a, b, c, d, haldane=True):
    cells = [a, b, c, d]
    used = min(cells) == 0 and haldane
    if used:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = np.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    return orr, np.exp(np.log(orr) - 1.96 * se), np.exp(np.log(orr) + 1.96 * se), used


def contrast(g1, g2, marker, label):
    a = int(np.sum(g1 & marker)); b = int(np.sum(g1 & ~marker))
    c = int(np.sum(g2 & marker)); d = int(np.sum(g2 & ~marker))
    _, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    orr, lo, hi, used = or_ci(a, b, c, d)
    return dict(test=label, OR=orr, ci_lo=lo, ci_hi=hi, p=p,
                note=f"{a}/{a+b} vs {c}/{c+d}", haldane=used)


def zscore_deg(genes):
    """Compute the iteration14 z-score DEG calls for `genes` (a gene-symbol array).

    Returns (cts, Z, deg, expr_total) where Z is (n_genes x n_celltypes) row-z-scored mean
    log1p(CP10K), deg is the boolean (z > Z_THR) matrix, and expr_total is the per-gene
    TPM-like total = SUM over cell types of mean LINEAR CP10K (iteration14 script 19, used
    for the EXPR_THR low-expression pre-filter). Genes absent from the annotation get an
    all-NaN Z row, all-False deg row, and expr_total 0 (caller filters the universe to
    var_names first, so in practice every gene is present)."""
    log(f"reading {H5AD} (backed) for expression ...")
    A = ad.read_h5ad(H5AD, backed="r")
    vn = pd.Index(A.var_names)
    ct = A.obs[CT_COL].astype(str).values
    cts = [c for c in CT_ORDER if c in set(ct)] + [c for c in pd.unique(ct) if c not in CT_ORDER]

    idx = vn.get_indexer(genes)              # -1 for genes not in annotation
    present = idx >= 0
    X = A[:, idx[present]].to_memory().X
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

    expr_present = np.zeros((int(present.sum()), len(cts)))      # mean log1p(CP10K)
    lin_present = np.zeros((int(present.sum()), len(cts)))       # mean linear CP10K
    for j, c in enumerate(cts):
        m = ct == c
        if int(m.sum()) == 0:
            continue
        cp = X[m].multiply(1e4 / lib[m][:, None]).tocsr()        # linear CP10K
        lin_present[:, j] = np.asarray(cp.mean(0)).ravel()
        cp.data = np.log1p(cp.data)
        expr_present[:, j] = np.asarray(cp.mean(0)).ravel()

    # scatter back to full gene set
    expr = np.full((len(genes), len(cts)), np.nan)
    expr[present] = expr_present
    expr_total = np.zeros(len(genes))
    expr_total[present] = lin_present.sum(1)                     # TPM-like total
    mu = np.nanmean(expr, 1, keepdims=True)
    sd = np.nanstd(expr, 1, keepdims=True)
    Z = np.where(sd > 0, (expr - mu) / sd, 0.0)
    Z[~present] = np.nan
    deg = (Z > Z_THR) & present[:, None]
    return cts, Z, deg, expr_total


def main():
    cls = pd.read_csv(CLASS, sep="\t")
    det = set(ad.read_h5ad(H5AD, backed="r").var_names)
    u = cls[cls.gene.isin(det)].reset_index(drop=True)
    log(f"cleaned 4-group catalog: {len(cls)}; PBMC-detectable: {len(u)}")

    genes = u.gene.values
    cts, Z, deg, expr_total = zscore_deg(genes)

    # ---- low-expression pre-filter (TPM-like total of linear CP10K) ----
    u["expr_total_cp10k"] = expr_total
    if EXPR_THR > 0:
        keep = expr_total > EXPR_THR
        log(f"low-expression filter EXPR_THR={EXPR_THR:g} (sum of mean linear CP10K over "
            f"cell types): {int(keep.sum())}/{len(u)} TFs survive, {int((~keep).sum())} dropped")
        u = u[keep].reset_index(drop=True)
        Z, deg, genes = Z[keep], deg[keep], genes[keep]
    log(f"universe after expression filter: {len(u)}")

    # ---- z-score DEG labels (this is the DEG definition swap) ----
    u["n_DEG_celltypes"] = deg.sum(1)
    u["is_DEG"] = u["n_DEG_celltypes"] > 0
    for j, c in enumerate(cts):
        u[f"z_{c}"] = Z[:, j]
    u["peak_celltype"] = [cts[int(np.nanargmax(Z[i]))] if np.isfinite(Z[i]).any() else ""
                          for i in range(len(genes))]
    u.to_csv(OUT_PER, sep="\t", index=False)
    log(f"z-score DEG (z>{Z_THR:g} in >=1 cell type): {int(u.is_DEG.sum())}/{len(u)} TFs "
        f"= {u.is_DEG.mean():.1%}")

    marker = u.is_DEG.values          # <-- the only line that differs in spirit from script 07
    grp = u.dup_class4.values
    prov = u.provenance.values
    arr = u.arrangement.values

    summ = []
    # ---- (1) per-category one-vs-rest ----
    ovr = []
    for g in GROUPS:
        m = grp == g
        r = contrast(m, ~m, marker, f"onevsrest__{g}")
        r["marker_rate"] = float(marker[m].mean())
        r["n"] = int(m.sum())
        ovr.append(r)
    ovr = pd.DataFrame(ovr)
    ovr["fdr"] = multipletests(ovr.p, method="fdr_bh")[1]
    for _, r in ovr.iterrows():
        summ.append(dict(test=r.test, OR=r.OR, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p,
                         fdr=r.fdr, note=f"n={r.n}, DEG_rate={r.marker_rate:.3f} ({r.note})"))

    # ---- (2) key pairwise contrasts ----
    pair_defs = [
        (grp == "clustered_cis_ohnolog", grp == "clustered_SSD_tandem",
         "pair__within_clustered__cis_ohnolog_vs_SSD_tandem"),
        (grp == "dispersed_ohnolog", grp == "dispersed_SSD",
         "pair__within_dispersed__ohnolog_vs_SSD"),
        (grp == "clustered_cis_ohnolog", grp == "dispersed_ohnolog",
         "pair__within_ohnolog__clustered_vs_dispersed"),
        (grp == "clustered_SSD_tandem", grp == "dispersed_SSD",
         "pair__within_SSD__clustered_vs_dispersed"),
        (prov == "ohnolog", prov == "SSD", "marginal__ohnolog_vs_SSD"),
        (arr == "clustered", arr == "dispersed", "marginal__clustered_vs_dispersed"),
    ]
    for g1, g2, lab in pair_defs:
        r = contrast(g1, g2, marker, lab)
        summ.append(dict(test=r["test"], OR=r["OR"], ci_lo=r["ci_lo"], ci_hi=r["ci_hi"],
                         p=r["p"], fdr=np.nan, note=r["note"]))

    # ---- (3) factorial logistic: DEG ~ ohnolog * clustered ----
    try:
        import statsmodels.formula.api as smf
        d = pd.DataFrame({"marker": marker.astype(int),
                          "ohnolog": (prov == "ohnolog").astype(int),
                          "clustered": (arr == "clustered").astype(int)})
        fit = smf.logit("marker ~ ohnolog * clustered", data=d).fit(disp=0)
        for term in ["ohnolog", "clustered", "ohnolog:clustered"]:
            orr = float(np.exp(fit.params[term]))
            lo, hi = np.exp(fit.conf_int().loc[term])
            summ.append(dict(test=f"logit__{term}", OR=orr, ci_lo=float(lo), ci_hi=float(hi),
                             p=float(fit.pvalues[term]), fdr=np.nan,
                             note="adjusted OR (provenance x arrangement factorial)"))
    except Exception as e:
        log("logit failed:", e)

    sdf = pd.DataFrame(summ)
    sdf.to_csv(OUT_SUM, sep="\t", index=False)

    # ---- (4) per-cell-type DEG rate by category ----
    ctrows = []
    for j, c in enumerate(cts):
        mk = deg[:, j]
        row = dict(cell_type=c)
        for g in GROUPS:
            m = grp == g
            row[g] = int(np.sum(m & mk))
            row[g + "_rate"] = float(mk[m].mean()) if m.sum() else np.nan
        ctrows.append(row)
    ctdf = pd.DataFrame(ctrows)
    ctdf.to_csv(OUT_CT, sep="\t", index=False)

    log("per-category one-vs-rest (z-score DEG):")
    print(ovr[["test", "n", "marker_rate", "OR", "ci_lo", "ci_hi", "p", "fdr"]].to_string(index=False))
    log("key contrasts + factorial:")
    print(sdf[sdf.test.str.startswith(("pair", "marginal", "logit"))]
          [["test", "OR", "ci_lo", "ci_hi", "p", "note"]].to_string(index=False))

    make_figure(u, ovr, sdf, ctdf, cts)
    log(f"wrote {OUT_SUM} / {OUT_CT} / {OUT_PER}")


def make_figure(u, ovr, sdf, ctdf, cell_types):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6))

    # A: DEG rate by 4 category
    ax = axes[0]
    o = ovr.set_index("test")
    rates = [o.loc[f"onevsrest__{g}", "marker_rate"] for g in GROUPS]
    ns = [int(o.loc[f"onevsrest__{g}", "n"]) for g in GROUPS]
    ors = [o.loc[f"onevsrest__{g}", "OR"] for g in GROUPS]
    fdr = [o.loc[f"onevsrest__{g}", "fdr"] for g in GROUPS]
    ax.bar(range(4), rates, color=[COLORS[g] for g in GROUPS])
    for i, (rt, n, orr, q) in enumerate(zip(rates, ns, ors, fdr)):
        star = "*" if q < 0.05 else ""
        ax.text(i, rt + 0.004, f"{rt:.2f}{star}\nn={n}\nOR={orr:.2f}", ha="center", fontsize=8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["dispersed\nohnolog", "clustered\ncis-ohnolog",
                        "clustered\nSSD-tandem", "dispersed\nSSD"], fontsize=8.5)
    ax.set_ylabel(f"fraction = z-score DEG (z>{Z_THR:g}, any cell type)")
    ax.set_title("(A) z-score DEG rate by 4 categories\n(one-vs-rest OR; * FDR<0.05)", fontsize=10)
    ax.set_ylim(0, max(rates) * 1.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # B: 2x2 factorial grid (provenance x arrangement), cell = DEG rate
    ax = axes[1]
    grid = np.zeros((2, 2)); ngrid = np.zeros((2, 2), int)
    prov_levels = ["ohnolog", "SSD"]; arr_levels = ["clustered", "dispersed"]
    for i, pv in enumerate(prov_levels):
        for j, ar in enumerate(arr_levels):
            m = (u.provenance == pv) & (u.arrangement == ar)
            grid[i, j] = u.loc[m, "is_DEG"].mean() if m.sum() else np.nan
            ngrid[i, j] = int(m.sum())
    im = ax.imshow(grid, cmap="Greens", vmin=0, vmax=np.nanmax(grid) * 1.05)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{grid[i,j]:.2f}\n(n={ngrid[i,j]})", ha="center", va="center", fontsize=11)
    ax.set_xticks([0, 1]); ax.set_xticklabels(arr_levels)
    ax.set_yticks([0, 1]); ax.set_yticklabels(prov_levels)
    ax.set_xlabel("arrangement"); ax.set_ylabel("provenance")
    ax.set_title("(B) z-score DEG rate: provenance x arrangement", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("DEG rate", fontsize=8)

    # C: forest of key contrasts
    ax = axes[2]
    keys = ["pair__within_clustered__cis_ohnolog_vs_SSD_tandem",
            "pair__within_dispersed__ohnolog_vs_SSD",
            "pair__within_ohnolog__clustered_vs_dispersed",
            "pair__within_SSD__clustered_vs_dispersed",
            "marginal__ohnolog_vs_SSD", "marginal__clustered_vs_dispersed"]
    labels = ["cis-ohno vs SSD-tandem\n(within clustered)",
              "ohnolog vs SSD\n(within dispersed)",
              "clustered vs dispersed\n(within ohnolog)",
              "clustered vs dispersed\n(within SSD)",
              "ohnolog vs SSD (marginal)", "clustered vs dispersed (marginal)"]
    s = sdf.set_index("test")
    y = np.arange(len(keys))[::-1]
    orv = [s.loc[k, "OR"] for k in keys]
    lo = [s.loc[k, "ci_lo"] for k in keys]
    hi = [s.loc[k, "ci_hi"] for k in keys]
    pv = [s.loc[k, "p"] for k in keys]
    colors = ["#2166ac" if v >= 1 else "#b2182b" for v in orv]
    ax.hlines(y, lo, hi, color=colors, lw=1.6)
    ax.scatter(orv, y, c=colors, s=45, edgecolor="k", linewidth=0.5, zorder=3)
    ax.axvline(1, color="grey", ls="--", lw=1)
    ax.set_xscale("log", base=2)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    for yi, v, h, p in zip(y, orv, hi, pv):
        ax.text(h * 1.05, yi, f"OR={v:.2f}{'*' if p<0.05 else ''}", va="center", fontsize=7.5)
    ax.set_xlabel("OR (first vs second), log2")
    ax.set_title("(C) Key contrasts: is it provenance or arrangement?", fontsize=10)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)

    efilt = f"; expr>{EXPR_THR:g} CP10K-sum pre-filter" if EXPR_THR > 0 else ""
    fig.suptitle(f"PBMC z-score DEG (iteration14 def, z>{Z_THR:g}{efilt}) enrichment across the 4 "
                 "duplication categories (Lambert TFs; KRAB/HOX/readthrough/pseudo removed)",
                 fontsize=11.5, y=1.02)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT_FIG}.pdf/.png")


if __name__ == "__main__":
    main()
