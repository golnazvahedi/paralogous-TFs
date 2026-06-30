#!/usr/bin/env python3
"""
DEG (cell-type marker) enrichment OR across the FOUR duplication categories (script 06):
   dispersed_ohnolog | clustered_cis_ohnolog | clustered_SSD_tandem | dispersed_SSD

Marker = CANONICAL one-vs-rest PBMC lineage marker (results/intermediate/deg/
human_markers.major_lineage.top.tsv) of >=1 major_lineage cell type -- the identity-DEG
set.
Universe = the cleaned 4-group catalog (Lambert TFs, KRAB/HOX/readthrough/pseudo removed)
intersected with PBMC-detectable genes (read directly from human_pbmc.h5ad var_names,
backed -- matching the spleen script 08, so all three OR scripts share one universe rule).

Computes:
  (1) per-category marker rate + ONE-VS-REST Fisher OR (+95% CI, p, BH-FDR);
  (2) the key PAIRWISE contrasts (does provenance help within each arrangement, and
      does arrangement matter within each provenance);
  (3) a provenance x arrangement FACTORIAL logistic regression (adjusted ORs +
      interaction) -- the clean decomposition of the 2x2;
  (4) per-cell-type marker rate by category (for the figure).

Outputs:
  results/TF_dup4_DEG_OR.summary.tsv        (one-vs-rest + pairwise + factorial)
  results/TF_dup4_DEG_OR.by_celltype.tsv
  results/TF_dup4_DEG_OR.per_TF.tsv
  results/figures/TF_dup4_DEG_OR.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/07_human_pbmc_dup4_DEG_OR.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
H5AD = RES / "intermediate" / "h5ad" / "human_pbmc.h5ad"   # detectable universe = PBMC var_names
TOP = RES / "intermediate" / "deg" / "human_markers.major_lineage.top.tsv"
OUT_SUM = RES / "TF_dup4_DEG_OR.summary.tsv"
OUT_CT = RES / "TF_dup4_DEG_OR.by_celltype.tsv"
OUT_PER = RES / "TF_dup4_DEG_OR.per_TF.tsv"
OUT_FIG = RES / "figures" / "TF_dup4_DEG_OR"

GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
COLORS = {"dispersed_ohnolog": "#2166ac", "clustered_cis_ohnolog": "#5aae61",
          "clustered_SSD_tandem": "#762a83", "dispersed_SSD": "#bababa"}


def log(*a):
    print("[07]", *a, flush=True)


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


def main():
    cls = pd.read_csv(CLASS, sep="\t")
    log(f"reading PBMC var_names (backed) from {H5AD} ...")
    det = set(ad.read_h5ad(H5AD, backed="r").var_names)
    u = cls[cls.gene.isin(det)].reset_index(drop=True)
    log(f"cleaned 4-group catalog: {len(cls)}; PBMC-detectable: {len(u)}")

    top = pd.read_csv(TOP, sep="\t")
    cell_types = sorted(top.cell_type.unique())
    marker_by_ct = {ct: set(top.loc[top.cell_type == ct, "gene"]) for ct in cell_types}
    any_marker = set(top.gene)

    genes = u.gene.values
    u["is_marker"] = np.isin(genes, list(any_marker))
    u["n_marker_celltypes"] = [sum(g in marker_by_ct[ct] for ct in cell_types) for g in genes]
    u.to_csv(OUT_PER, sep="\t", index=False)

    marker = u.is_marker.values
    grp = u.dup_class4.values

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
                         fdr=r.fdr, note=f"n={r.n}, marker_rate={r.marker_rate:.3f} ({r.note})"))

    prov = u.provenance.values
    arr = u.arrangement.values
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

    # ---- (3) factorial logistic: marker ~ ohnolog * clustered ----
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

    # ---- (4) per-cell-type marker rate by category ----
    ctrows = []
    for ct in cell_types:
        mk = np.isin(genes, list(marker_by_ct[ct]))
        row = dict(cell_type=ct)
        for g in GROUPS:
            m = grp == g
            row[g] = int(np.sum(m & mk))
            row[g + "_rate"] = float(mk[m].mean()) if m.sum() else np.nan
        ctrows.append(row)
    ctdf = pd.DataFrame(ctrows)
    ctdf.to_csv(OUT_CT, sep="\t", index=False)

    log("per-category one-vs-rest:")
    print(ovr[["test", "n", "marker_rate", "OR", "ci_lo", "ci_hi", "p", "fdr"]].to_string(index=False))
    log("key contrasts + factorial:")
    print(sdf[sdf.test.str.startswith(("pair", "marginal", "logit"))]
          [["test", "OR", "ci_lo", "ci_hi", "p", "note"]].to_string(index=False))

    make_figure(u, ovr, sdf, ctdf, cell_types)
    log(f"wrote {OUT_SUM} / {OUT_CT} / {OUT_PER}")


def make_figure(u, ovr, sdf, ctdf, cell_types):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6))

    # A: marker rate by 4 category
    ax = axes[0]
    o = ovr.set_index("test")
    rates = [o.loc[f"onevsrest__{g}", "marker_rate"] for g in GROUPS]
    ns = [int(o.loc[f"onevsrest__{g}", "n"]) for g in GROUPS]
    ors = [o.loc[f"onevsrest__{g}", "OR"] for g in GROUPS]
    fdr = [o.loc[f"onevsrest__{g}", "fdr"] for g in GROUPS]
    bars = ax.bar(range(4), rates, color=[COLORS[g] for g in GROUPS])
    for i, (rt, n, orr, q) in enumerate(zip(rates, ns, ors, fdr)):
        star = "*" if q < 0.05 else ""
        ax.text(i, rt + 0.004, f"{rt:.2f}{star}\nn={n}\nOR={orr:.2f}", ha="center", fontsize=8)
    ax.set_xticks(range(4))
    ax.set_xticklabels(["dispersed\nohnolog", "clustered\ncis-ohnolog",
                        "clustered\nSSD-tandem", "dispersed\nSSD"], fontsize=8.5)
    ax.set_ylabel("fraction = lineage marker (any cell type)")
    ax.set_title("(A) Marker rate by 4 categories\n(one-vs-rest OR; * FDR<0.05)", fontsize=10)
    ax.set_ylim(0, max(rates) * 1.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # B: 2x2 factorial grid (provenance x arrangement), cell = marker rate
    ax = axes[1]
    grid = np.zeros((2, 2)); ngrid = np.zeros((2, 2), int)
    prov_levels = ["ohnolog", "SSD"]; arr_levels = ["clustered", "dispersed"]
    for i, pv in enumerate(prov_levels):
        for j, ar in enumerate(arr_levels):
            m = (u.provenance == pv) & (u.arrangement == ar)
            grid[i, j] = u.loc[m, "is_marker"].mean() if m.sum() else np.nan
            ngrid[i, j] = int(m.sum())
    im = ax.imshow(grid, cmap="Greens", vmin=0, vmax=np.nanmax(grid) * 1.05)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{grid[i,j]:.2f}\n(n={ngrid[i,j]})", ha="center", va="center", fontsize=11)
    ax.set_xticks([0, 1]); ax.set_xticklabels(arr_levels)
    ax.set_yticks([0, 1]); ax.set_yticklabels(prov_levels)
    ax.set_xlabel("arrangement"); ax.set_ylabel("provenance")
    ax.set_title("(B) Marker rate: provenance x arrangement", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("marker rate", fontsize=8)

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
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.suptitle("PBMC lineage-marker DEG enrichment across the 4 duplication categories "
                 "(Lambert TFs; KRAB/HOX/readthrough/pseudo removed)", fontsize=11.5, y=1.02)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT_FIG}.pdf/.png")


if __name__ == "__main__":
    main()
