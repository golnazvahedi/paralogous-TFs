#!/usr/bin/env python3
"""
GWAS (polygenic complex-immune-disease) enrichment across the four duplication categories --
the polygenic counterpart of the monogenic-IEI analysis (script 13).

Curated GWAS genes from script 14 (results/gwas_immune_disease_genes.long.tsv): genome-wide-
significant (P<=5e-8) nearest/overlapping genes for SLE_lupus, rheumatoid_arthritis, allergy,
asthma, type_1_diabetes. We use the PROTEIN-CODING, non-MHC subset and ask, per disease and for
the UNION ("any of 5"), whether the dup_class4 categories are enriched among paralogous TFs.

Matching = canonical HGNC symbol on both sides (catalog gene and GWAS gene normalized through
inputs/hgnc_complete_set.txt symbol/alias/prev), to absorb symbol drift. Universe = the full
classified catalog (results/TF_dup_2x2_classification.tsv); one-vs-rest Fisher OR WITHIN that
universe. Lineage (12), gene-origin window (05d) and the monogenic-IEI flag (13) are merged so
we can directly contrast POLYGENIC (GWAS) vs MONOGENIC (IEI) architecture.

Outputs:
  results/TF_GWAS_disease_genes.tsv        per-TF: catalog annot + per-disease GWAS flags + IEI flag
  results/TF_GWAS_enrichment.summary.tsv   per-category & contrast ORs (per disease + any)
  results/figures/TF_GWAS_disease_enrichment.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/15_gwas_TF_enrichment.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
HGNC = ROOT / "inputs" / "hgnc_complete_set.txt"
GWAS = RES / "gwas_immune_disease_genes.long.tsv"
ALPHABETA = RES / "TF_ohnolog_alpha_beta.tsv"
ORIGIN = RES / "TF_gene_origin_age_2R.full.tsv"
IEI = RES / "TF_IEI_disease_genes.tsv"
OUT_PER = RES / "TF_GWAS_disease_genes.tsv"
OUT_SUM = RES / "TF_GWAS_enrichment.summary.tsv"
OUT_FIG = RES / "figures" / "TF_GWAS_disease_enrichment"

DISEASES = ["SLE_lupus", "rheumatoid_arthritis", "allergy", "asthma", "type_1_diabetes"]
DLABEL = {"SLE_lupus": "SLE/lupus", "rheumatoid_arthritis": "RA", "allergy": "allergy",
          "asthma": "asthma", "type_1_diabetes": "T1D"}
GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
GLABEL = {"dispersed_ohnolog": "dispersed\nohnolog", "clustered_cis_ohnolog": "clustered\ncis-ohnolog",
          "clustered_SSD_tandem": "clustered\nSSD-tandem", "dispersed_SSD": "dispersed\nSSD"}
COL = {"dispersed_ohnolog": "#2166ac", "clustered_cis_ohnolog": "#5aae61",
       "clustered_SSD_tandem": "#762a83", "dispersed_SSD": "#bababa"}


def log(*a):
    print("[15]", *a, flush=True)


def canon_map():
    """UPPER(any symbol/alias/prev) -> canonical HGNC symbol."""
    h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
    m = {}
    for _, r in h.iterrows():
        canon = r["symbol"]
        if not isinstance(canon, str):
            continue
        for col in ("symbol", "alias_symbol", "prev_symbol"):
            v = r.get(col)
            if isinstance(v, str) and v:
                for s in v.split("|"):
                    m.setdefault(s.upper(), canon)
    return m


def or_ci(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return orr, np.exp(np.log(orr) - 1.96 * se), np.exp(np.log(orr) + 1.96 * se)


def contrast(g1, g2, marker, label):
    a, b = int(np.sum(g1 & marker)), int(np.sum(g1 & ~marker))
    c, d = int(np.sum(g2 & marker)), int(np.sum(g2 & ~marker))
    _, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    orr, lo, hi = or_ci(a, b, c, d)
    return dict(test=label, OR=orr, ci_lo=lo, ci_hi=hi, p=p,
                rate1=a / (a + b) if a + b else np.nan, rate2=c / (c + d) if c + d else np.nan,
                note=f"{a}/{a+b} vs {c}/{c+d}")


def main():
    cm = canon_map()
    cls = pd.read_csv(CLASS, sep="\t")
    cls = cls[cls.dup_class4.isin(GROUPS)].reset_index(drop=True)
    cls["canon"] = cls.gene.str.upper().map(lambda s: cm.get(s, s))

    g = pd.read_csv(GWAS, sep="\t")
    g = g[g.is_protein_coding & (~g.is_MHC)].copy()
    g["canon"] = g.gene.str.upper().map(lambda s: cm.get(s, s))
    dis_sets = {d: set(g.loc[g.disease == d, "canon"]) for d in DISEASES}
    any_set = set().union(*dis_sets.values())

    for d in DISEASES:
        cls[f"gwas_{d}"] = cls.canon.isin(dis_sets[d])
    cls["gwas_any"] = cls.canon.isin(any_set)
    cls["n_gwas_diseases"] = cls[[f"gwas_{d}" for d in DISEASES]].sum(1)

    # merge lineage, origin, IEI
    if ALPHABETA.exists():
        cls = cls.merge(pd.read_csv(ALPHABETA, sep="\t")[["gene", "lineage"]], on="gene", how="left")
    if ORIGIN.exists():
        cls = cls.merge(pd.read_csv(ORIGIN, sep="\t")[["gene", "origin_window"]], on="gene", how="left")
    if IEI.exists():
        cls = cls.merge(pd.read_csv(IEI, sep="\t")[["gene", "is_IEI"]], on="gene", how="left")
    cls.to_csv(OUT_PER, sep="\t", index=False)
    log(f"paralogous TFs that are GWAS genes (any of 5, protein-coding non-MHC): "
        f"{int(cls.gwas_any.sum())}/{len(cls)} ({cls.gwas_any.mean():.1%})")

    grp = cls.dup_class4.values
    arr = cls.arrangement.values
    prov = cls.provenance.values
    rows = []

    # per-category one-vs-rest, for each disease and the union
    for col, mname in [(f"gwas_{d}", d) for d in DISEASES] + [("gwas_any", "any")]:
        mk = cls[col].values
        for grp_name in GROUPS:
            m = grp == grp_name
            r = contrast(m, ~m, mk, f"onevsrest__{grp_name}")
            r["modality"] = mname; r["n"] = int(m.sum()); rows.append(r)

    # headline contrasts on the union
    mk = cls.gwas_any.values
    for g1, g2, lab in [
        (arr == "clustered", arr == "dispersed", "buffered_clustered_vs_unbuffered_dispersed"),
        (prov == "ohnolog", prov == "SSD", "ohnolog_vs_SSD"),
        ((cls.lineage == "alpha").values, (cls.lineage == "beta").values, "alpha_vs_beta"),
        ((cls.origin_window == "before_R1").values, (cls.origin_window != "before_R1").values,
         "preR1_ancient_vs_younger"),
    ]:
        r = contrast(g1, g2, mk, lab); r["modality"] = "any"; r["n"] = int(g1.sum() + g2.sum())
        rows.append(r)

    sdf = pd.DataFrame(rows)
    sdf.to_csv(OUT_SUM, sep="\t", index=False)

    log("per-category GWAS rate + one-vs-rest OR (union of 5 diseases):")
    print(sdf[(sdf.modality == "any") & sdf.test.str.startswith("onevsrest")]
          [["test", "n", "rate1", "OR", "ci_lo", "ci_hi", "p"]].to_string(index=False))
    log("headline contrasts (GWAS union):")
    print(sdf[sdf.test.str.contains("vs_|_vs_")][["test", "rate1", "rate2", "OR", "ci_lo", "ci_hi", "p"]]
          .to_string(index=False))

    # monogenic vs polygenic side-by-side per category
    if "is_IEI" in cls.columns:
        log("MONOGENIC (IEI) vs POLYGENIC (GWAS-any) rate by category:")
        for grp_name in GROUPS:
            m = cls.dup_class4 == grp_name
            print(f"  {grp_name}: IEI {cls.loc[m,'is_IEI'].mean():.3f} | GWAS {cls.loc[m,'gwas_any'].mean():.3f}")
        # TFs that are GWAS but NOT monogenic IEI (the polygenic-only paralogous TFs)
        poly = cls[cls.gwas_any & (~cls.is_IEI.fillna(False))]
        log(f"GWAS-only (polygenic, not monogenic-IEI) paralogous TFs: {len(poly)}")
        for grp_name in GROUPS:
            gg = sorted(poly.loc[poly.dup_class4 == grp_name, "gene"])
            print(f"  {grp_name} ({len(gg)}): {', '.join(gg[:30])}{' ...' if len(gg)>30 else ''}")

    make_figure(cls, sdf)
    log(f"wrote {OUT_PER} / {OUT_SUM}")


def make_figure(cls, sdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16.5, 5.2),
                                        gridspec_kw={"width_ratios": [1.05, 1.15, 1.1]})
    s = sdf.set_index(["modality", "test"])

    # (A) GWAS-any rate by category + IEI overlay
    x = np.arange(len(GROUPS))
    rates = [s.loc[("any", f"onevsrest__{g}"), "rate1"] for g in GROUPS]
    ns = [int(s.loc[("any", f"onevsrest__{g}"), "n"]) for g in GROUPS]
    ors = [s.loc[("any", f"onevsrest__{g}"), "OR"] for g in GROUPS]
    ps = [s.loc[("any", f"onevsrest__{g}"), "p"] for g in GROUPS]
    iei = [cls.loc[cls.dup_class4 == g, "is_IEI"].mean() if "is_IEI" in cls else np.nan for g in GROUPS]
    axA.bar(x, rates, color=[COL[g] for g in GROUPS], width=0.7, label="GWAS (any of 5)")
    axA.bar(x, iei, color="black", width=0.3, label="monogenic IEI")
    for i, (rt, n, orr, p) in enumerate(zip(rates, ns, ors, ps)):
        axA.text(i, rt + 0.01, f"{rt:.2f}{'*' if p<0.05 else ''}\nOR={orr:.2f}", ha="center", fontsize=8)
    axA.set_xticks(x); axA.set_xticklabels([GLABEL[g] for g in GROUPS], fontsize=8.5)
    axA.set_ylabel("fraction = disease gene")
    axA.set_title("(A) Polygenic GWAS vs monogenic IEI rate\nby category (one-vs-rest OR; * p<0.05)", fontsize=9.5)
    axA.set_ylim(0, max(rates) * 1.4); axA.legend(fontsize=7.5, frameon=False)
    for sp in ("top", "right"):
        axA.spines[sp].set_visible(False)

    # (B) disease x category rate heatmap
    M = np.array([[s.loc[(d, f"onevsrest__{g}"), "rate1"] for g in GROUPS] for d in DISEASES])
    im = axB.imshow(M, cmap="Reds", aspect="auto", vmin=0, vmax=np.nanmax(M) * 1.05)
    for i in range(len(DISEASES)):
        for j in range(len(GROUPS)):
            axB.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=8.5,
                     color="white" if M[i, j] > np.nanmax(M) * 0.6 else "black")
    axB.set_xticks(range(len(GROUPS))); axB.set_xticklabels([GLABEL[g] for g in GROUPS], fontsize=8)
    axB.set_yticks(range(len(DISEASES))); axB.set_yticklabels([DLABEL[d] for d in DISEASES], fontsize=9)
    axB.set_title("(B) GWAS-gene rate: disease x category", fontsize=9.5)
    fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04).set_label("GWAS rate", fontsize=8)

    # (C) forest of key contrasts (GWAS union)
    keys = ["buffered_clustered_vs_unbuffered_dispersed", "ohnolog_vs_SSD",
            "preR1_ancient_vs_younger", "alpha_vs_beta"]
    labels = ["clustered (buffered) vs\ndispersed (unbuffered)", "ohnolog vs SSD",
              "pre-R1 ancient vs younger", "alpha vs beta"]
    y = np.arange(len(keys))[::-1]
    orv = [s.loc[("any", k), "OR"] for k in keys]; lo = [s.loc[("any", k), "ci_lo"] for k in keys]
    hi = [s.loc[("any", k), "ci_hi"] for k in keys]; pv = [s.loc[("any", k), "p"] for k in keys]
    colors = ["#b2182b" if v >= 1 else "#2166ac" for v in orv]
    axC.hlines(y, lo, hi, color=colors, lw=1.8)
    axC.scatter(orv, y, c=colors, s=46, edgecolor="k", linewidth=0.5, zorder=3)
    axC.axvline(1, color="grey", ls="--", lw=1); axC.set_xscale("log", base=2)
    axC.set_yticks(y); axC.set_yticklabels(labels, fontsize=8.5)
    for yi, v, h, p in zip(y, orv, hi, pv):
        axC.text(h * 1.04, yi, f"OR={v:.2f}{'*' if p<0.05 else ''}", va="center", fontsize=8)
    axC.set_xlabel("odds ratio (log2)")
    axC.set_title("(C) Key GWAS-enrichment contrasts (union)", fontsize=9.5)
    for sp in ("top", "right"):
        axC.spines[sp].set_visible(False)

    fig.suptitle("Polygenic GWAS complex-immune-disease enrichment across the paralogous-TF "
                 "duplication framework (GWAS Catalog P<=5e-8; protein-coding non-MHC; n=1,153 TFs)",
                 fontsize=11.5, y=1.02)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT_FIG}.pdf/.png")


if __name__ == "__main__":
    main()
