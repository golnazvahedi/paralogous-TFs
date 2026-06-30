#!/usr/bin/env python3
"""
Monogenic immune-disease (Inborn Errors of Immunity, IEI) enrichment across the duplication
framework: do the four dup_class4 categories -- and specifically BUFFERED (clustered, cis-
neighbour) vs UNBUFFERED (dispersed/singleton) paralogous TFs -- differ in how often they are
curated monogenic immune-disease genes?

IEI gene set = Genomics England PanelApp, distilled to inputs/iei/panelapp_iei_genes.tsv:
  * panel 398 "Primary immunodeficiency or monogenic inflammatory bowel disease" (v9.18)
  * panel 1075 "Autoinflammatory disorders" (v3.10)
PanelApp confidence: 3=GREEN (diagnostic-grade, our primary set), 2=amber, 1=red. We also flag
GREEN genes with a MONOALLELIC / dominant mode of inheritance -- the haploinsufficiency-
compatible subset most directly relevant to the dosage thesis.

Matching = HGNC ID (catalog symbol -> hgnc_id via inputs/hgnc_complete_set.txt incl. alias/prev;
panel genes carry hgnc_id), with a symbol fallback. Universe = the full classified paralogous-
TF catalog (results/TF_dup_2x2_classification.tsv); enrichment is one-vs-rest WITHIN that
universe (which kind of paralogous TF is most disease-associated), via Fisher exact (OR, CI, p).
The headline contrast is clustered (buffered) vs dispersed (unbuffered). Alpha/beta lineage
(script 12) and gene-origin window (script 05d) are merged for secondary breakdowns.

Outputs:
  results/TF_IEI_disease_genes.tsv        per-TF: catalog annot + lineage + origin + IEI flags + MOI
  results/TF_IEI_enrichment.summary.tsv   all contrasts
  results/figures/TF_IEI_disease_enrichment.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/13_IEI_disease_gene_enrichment.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
IEI = ROOT / "inputs" / "iei" / "panelapp_iei_genes.tsv"
HGNC = ROOT / "inputs" / "hgnc_complete_set.txt"
ALPHABETA = RES / "TF_ohnolog_alpha_beta.tsv"
ORIGIN = RES / "TF_gene_origin_age_2R.full.tsv"
OUT_PER = RES / "TF_IEI_disease_genes.tsv"
OUT_SUM = RES / "TF_IEI_enrichment.summary.tsv"
OUT_FIG = RES / "figures" / "TF_IEI_disease_enrichment"

GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
GLABEL = {"dispersed_ohnolog": "dispersed\nohnolog", "clustered_cis_ohnolog": "clustered\ncis-ohnolog",
          "clustered_SSD_tandem": "clustered\nSSD-tandem", "dispersed_SSD": "dispersed\nSSD"}
COL = {"dispersed_ohnolog": "#2166ac", "clustered_cis_ohnolog": "#5aae61",
       "clustered_SSD_tandem": "#762a83", "dispersed_SSD": "#bababa"}


def log(*a):
    print("[13]", *a, flush=True)


def symbol_to_hgnc():
    h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
    m = {}
    for _, r in h.iterrows():
        hid = r["hgnc_id"]
        for col in ("symbol", "alias_symbol", "prev_symbol"):
            v = r.get(col)
            if isinstance(v, str) and v:
                for s in v.split("|"):
                    m.setdefault(s.upper(), hid)   # current symbol wins (added first per row order)
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
    cls = pd.read_csv(CLASS, sep="\t")
    cls = cls[cls.dup_class4.isin(GROUPS)].reset_index(drop=True)
    s2h = symbol_to_hgnc()
    cls["hgnc_id"] = cls.gene.str.upper().map(s2h)

    iei = pd.read_csv(IEI, sep="\t", dtype={"confidence_level": "Int64"})
    green = set(iei.loc[iei.confidence_level == 3, "hgnc_id"]) - {""}
    green_sym = set(iei.loc[iei.confidence_level == 3, "gene_symbol"].str.upper())
    amber_plus = set(iei.loc[iei.confidence_level >= 2, "hgnc_id"]) - {""}
    mono = set(iei.loc[(iei.confidence_level == 3) &
                       iei.mode_of_inheritance.str.contains("monoallelic", case=False, na=False),
                       "hgnc_id"]) - {""}
    log(f"IEI green genes: {len(green)} | green+amber: {len(amber_plus)} | green monoallelic: {len(mono)}")

    def flag(hidset, symset=None):
        f = cls.hgnc_id.isin(hidset)
        if symset:
            f = f | cls.gene.str.upper().isin(symset)
        return f.values

    cls["is_IEI"] = flag(green, green_sym)
    cls["is_IEI_amberplus"] = flag(amber_plus)
    cls["is_IEI_monoallelic"] = flag(mono)

    # merge lineage + origin window
    if ALPHABETA.exists():
        cls = cls.merge(pd.read_csv(ALPHABETA, sep="\t")[["gene", "lineage"]], on="gene", how="left")
    if ORIGIN.exists():
        cls = cls.merge(pd.read_csv(ORIGIN, sep="\t")[["gene", "origin_window"]], on="gene", how="left")

    cls.to_csv(OUT_PER, sep="\t", index=False)
    log(f"paralogous TFs that are GREEN IEI genes: {int(cls.is_IEI.sum())}/{len(cls)} "
        f"({cls.is_IEI.mean():.1%}); monoallelic subset: {int(cls.is_IEI_monoallelic.sum())}")

    grp = cls.dup_class4.values
    arr = cls.arrangement.values
    prov = cls.provenance.values
    iei_m = cls.is_IEI.values
    mono_m = cls.is_IEI_monoallelic.values

    rows = []
    # per-category one-vs-rest (green; and monoallelic-green)
    for mk, mname in [(iei_m, "IEI_green"), (mono_m, "IEI_green_monoallelic")]:
        for g in GROUPS:
            m = grp == g
            r = contrast(m, ~m, mk, f"onevsrest__{g}"); r["modality"] = mname
            r["n"] = int(m.sum()); rows.append(r)
    # headline contrasts (green)
    for g1, g2, lab in [
        (arr == "clustered", arr == "dispersed", "buffered_clustered_vs_unbuffered_dispersed"),
        (prov == "ohnolog", prov == "SSD", "ohnolog_vs_SSD"),
        ((cls.lineage == "alpha").values, (cls.lineage == "beta").values, "alpha_vs_beta"),
        ((cls.origin_window == "before_R1").values, (cls.origin_window != "before_R1").values,
         "preR1_ancient_vs_younger"),
    ]:
        r = contrast(g1, g2, iei_m, lab); r["modality"] = "IEI_green"; r["n"] = int(g1.sum() + g2.sum())
        rows.append(r)

    sdf = pd.DataFrame(rows)
    sdf.to_csv(OUT_SUM, sep="\t", index=False)
    log("per-category IEI-green rate + one-vs-rest OR:")
    print(sdf[(sdf.modality == "IEI_green") & sdf.test.str.startswith("onevsrest")]
          [["test", "n", "rate1", "OR", "ci_lo", "ci_hi", "p"]].to_string(index=False))
    log("headline contrasts (IEI-green):")
    print(sdf[sdf.test.str.contains("vs")][["test", "rate1", "rate2", "OR", "ci_lo", "ci_hi", "p", "note"]]
          .to_string(index=False))

    # list the actual IEI paralogous TFs per category
    log("IEI-green paralogous TFs by category:")
    for g in GROUPS:
        gg = sorted(cls.loc[(cls.dup_class4 == g) & cls.is_IEI, "gene"])
        print(f"  {g} ({len(gg)}): {', '.join(gg)}")

    make_figure(cls, sdf)
    log(f"wrote {OUT_PER} / {OUT_SUM}")


def make_figure(cls, sdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16, 5.2),
                                        gridspec_kw={"width_ratios": [1.1, 0.8, 1.2]})
    s = sdf[sdf.modality == "IEI_green"].set_index("test")
    sm = sdf[sdf.modality == "IEI_green_monoallelic"].set_index("test")

    # (A) IEI rate by category (green + monoallelic overlay)
    x = np.arange(len(GROUPS))
    rates = [s.loc[f"onevsrest__{g}", "rate1"] for g in GROUPS]
    mrates = [sm.loc[f"onevsrest__{g}", "rate1"] for g in GROUPS]
    ns = [int(s.loc[f"onevsrest__{g}", "n"]) for g in GROUPS]
    ors = [s.loc[f"onevsrest__{g}", "OR"] for g in GROUPS]
    ps = [s.loc[f"onevsrest__{g}", "p"] for g in GROUPS]
    axA.bar(x, rates, color=[COL[g] for g in GROUPS], width=0.7, label="all green IEI")
    axA.bar(x, mrates, color="black", width=0.32, label="monoallelic (dominant)")
    for i, (rt, n, orr, p) in enumerate(zip(rates, ns, ors, ps)):
        axA.text(i, rt + 0.004, f"{rt:.2f}{'*' if p<0.05 else ''}\nn={n}\nOR={orr:.2f}",
                 ha="center", fontsize=8)
    axA.set_xticks(x); axA.set_xticklabels([GLABEL[g] for g in GROUPS], fontsize=8.5)
    axA.set_ylabel("fraction = monogenic IEI gene (PanelApp green)")
    axA.set_title("(A) IEI-gene rate by duplication category\n(one-vs-rest OR; * p<0.05)", fontsize=10)
    axA.set_ylim(0, max(rates) * 1.45); axA.legend(fontsize=7.5, frameon=False)
    for sp in ("top", "right"):
        axA.spines[sp].set_visible(False)

    # (B) provenance x arrangement IEI-rate grid
    grid = np.zeros((2, 2)); ng = np.zeros((2, 2), int)
    for i, pv in enumerate(["ohnolog", "SSD"]):
        for j, ar in enumerate(["clustered", "dispersed"]):
            m = (cls.provenance == pv) & (cls.arrangement == ar)
            grid[i, j] = cls.loc[m, "is_IEI"].mean() if m.sum() else np.nan
            ng[i, j] = int(m.sum())
    im = axB.imshow(grid, cmap="Reds", vmin=0, vmax=np.nanmax(grid) * 1.1)
    for i in range(2):
        for j in range(2):
            axB.text(j, i, f"{grid[i,j]:.2f}\n(n={ng[i,j]})", ha="center", va="center", fontsize=11)
    axB.set_xticks([0, 1]); axB.set_xticklabels(["clustered\n(buffered)", "dispersed\n(unbuffered)"], fontsize=8.5)
    axB.set_yticks([0, 1]); axB.set_yticklabels(["ohnolog", "SSD"])
    axB.set_title("(B) IEI rate: provenance x arrangement", fontsize=10)
    fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04).set_label("IEI rate", fontsize=8)

    # (C) forest of key contrasts
    keys = ["buffered_clustered_vs_unbuffered_dispersed", "ohnolog_vs_SSD",
            "preR1_ancient_vs_younger", "alpha_vs_beta"]
    labels = ["clustered (buffered) vs\ndispersed (unbuffered)", "ohnolog vs SSD",
              "pre-R1 ancient vs younger", "alpha vs beta"]
    y = np.arange(len(keys))[::-1]
    orv = [s.loc[k, "OR"] for k in keys]; lo = [s.loc[k, "ci_lo"] for k in keys]
    hi = [s.loc[k, "ci_hi"] for k in keys]; pv = [s.loc[k, "p"] for k in keys]
    colors = ["#b2182b" if v >= 1 else "#2166ac" for v in orv]
    axC.hlines(y, lo, hi, color=colors, lw=1.8)
    axC.scatter(orv, y, c=colors, s=46, edgecolor="k", linewidth=0.5, zorder=3)
    axC.axvline(1, color="grey", ls="--", lw=1); axC.set_xscale("log", base=2)
    axC.set_yticks(y); axC.set_yticklabels(labels, fontsize=8.5)
    for yi, v, h, p in zip(y, orv, hi, pv):
        axC.text(h * 1.05, yi, f"OR={v:.2f}{'*' if p<0.05 else ''}", va="center", fontsize=8)
    axC.set_xlabel("odds ratio (log2)")
    axC.set_title("(C) Key IEI-enrichment contrasts", fontsize=10)
    for sp in ("top", "right"):
        axC.spines[sp].set_visible(False)

    fig.suptitle("Monogenic immune-disease (IEI) enrichment across the paralogous-TF duplication "
                 "framework (PanelApp green genes; n=1,153 classified TFs)", fontsize=11.5, y=1.02)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT_FIG}.pdf/.png")


if __name__ == "__main__":
    main()
