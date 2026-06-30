#!/usr/bin/env python3
"""
ALPHA vs BETA ohnolog DEG-enrichment, following Zhu et al. (Nature 2026, inputs/Zhu_Nature_2026.pdf).

Zhu et al. show the jawed-vertebrate 2R WGD was an ALLOPOLYPLOIDIZATION fusing two parental
lineages -- ALPHA and BETA -- with strongly ASYMMETRIC retention (alpha-derived ohnologs ~4x
more likely to be kept) and a stronger association of ALPHA ohnologs with cell-type marker
genes (DEGs). They assign alpha/beta via chicken orthology to ancestral paralogons.

Here we reproduce that split on our paralogous-TF catalog and ask the project's question:
does the ALPHA/BETA origin of an ohnolog TF predict cell-type IDENTITY (PBMC steady-state)
or CYTOKINE-RESPONSE DEG enrichment?

Alpha/beta source: Marletaz et al. hagfish paralogons table
(inputs/marletaz_paralogons/Vert_Evt_OGrrA.txt; column `1R` = alpha1/alpha2/beta1/beta2 per
species per orthogroup `FID`). We collapse alpha1/alpha2->alpha, beta1/beta2->beta and assign
each human TF a lineage by:
  (1) CHICKEN orthology by gene symbol (Galgal `Gname`) -- Zhu's stated method; then
  (2) ORTHOGROUP-level label (any species sharing the symbol maps to a single resolved OG)
      for TFs unmapped by (1). Where both exist they agree 100% (n=230).
Lineage is only meaningful for OHNOLOG-provenance TFs (provenance=='ohnolog' in dup_class4).

Universe = ohnolog-provenance TFs ∩ PBMC-detectable genes (human_pbmc.h5ad var_names), same
detectable-universe rule as scripts 07/10. Markers:
  PBMC identity   = results/intermediate/deg/human_markers.major_lineage.top.tsv  (script 07)
  cytokine resp.  = results/intermediate/human_cytokine_dict_DEGs.tsv             (script 09/10)

Stats per modality: alpha marker-rate, beta marker-rate, the ALPHA-vs-BETA Fisher OR (headline,
matched-provenance), and alpha-vs-rest / beta-vs-rest one-vs-rest ORs against the detectable
paralogous-TF background; dispersed_SSD shown as a reference (Zhu: SSD negative).

Outputs:
  results/TF_ohnolog_alpha_beta.tsv                 per-TF lineage + marker/responder flags
  results/TF_alpha_beta_DEG_OR.summary.tsv          all contrasts (both modalities)
  results/figures/TF_alpha_beta_DEG_OR.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/12_alpha_beta_ohnolog_DEG_OR.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
H5AD = RES / "intermediate" / "h5ad" / "human_pbmc.h5ad"
PARA = ROOT / "inputs" / "marletaz_paralogons" / "Vert_Evt_OGrrA.txt"
PBMC_MARK = RES / "intermediate" / "deg" / "human_markers.major_lineage.top.tsv"
CYT_DEGS = RES / "intermediate" / "human_cytokine_dict_DEGs.tsv"
OUT_PER = RES / "TF_ohnolog_alpha_beta.tsv"
OUT_SUM = RES / "TF_alpha_beta_DEG_OR.summary.tsv"
OUT_FIG = RES / "figures" / "TF_alpha_beta_DEG_OR"

LIN = {"alpha1": "alpha", "alpha2": "alpha", "beta1": "beta", "beta2": "beta"}


def log(*a):
    print("[12]", *a, flush=True)


def build_lineage_map():
    """Return dict UPPER(human symbol) -> 'alpha'/'beta' (chicken-primary, OG-union fallback)."""
    d = pd.read_csv(PARA, sep="\t")
    d["lineage"] = d["1R"].map(LIN)
    d["sym"] = d["Gname"].astype(str).str.upper()

    # (1) chicken (Galgal) symbol -> lineage
    g = d[(d.Sp == "Galgal") & d.Gname.notna() & d.lineage.notna()]
    chick = g.groupby("sym")["lineage"].agg(
        lambda s: "alpha" if set(s) == {"alpha"} else "beta" if set(s) == {"beta"} else "ambig")

    # (2) orthogroup-level label + symbol->OG across all species
    ogl = d.dropna(subset=["lineage"]).groupby("FID")["lineage"].agg(
        lambda s: "alpha" if set(s) == {"alpha"} else "beta" if set(s) == {"beta"} else "ambig")
    sym2og = d[d.Gname.notna()].groupby("sym")["FID"].agg(set)

    def og_label(sym):
        labs = {ogl.get(o) for o in sym2og.get(sym, set())}
        labs = {l for l in labs if l in ("alpha", "beta")}
        return next(iter(labs)) if len(labs) == 1 else None

    out = {}
    for sym in set(chick.index) | set(sym2og.index):
        v = chick.get(sym)
        if v not in ("alpha", "beta"):
            v = og_label(sym)
        if v in ("alpha", "beta"):
            out[sym] = v
    log(f"alpha/beta lineage map: {len(out)} chicken symbols "
        f"({sum(v=='alpha' for v in out.values())} alpha / {sum(v=='beta' for v in out.values())} beta)")
    return out


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
    return dict(test=label, OR=orr, ci_lo=lo, ci_hi=hi, p=p, note=f"{a}/{a+b} vs {c}/{c+d}")


def main():
    lin_map = build_lineage_map()
    cls = pd.read_csv(CLASS, sep="\t")
    det = set(ad.read_h5ad(H5AD, backed="r").var_names)
    u = cls[cls.gene.isin(det)].reset_index(drop=True)        # detectable paralogous TFs (all 4 groups)
    u["lineage"] = u.gene.str.upper().map(lin_map)
    u.loc[u.provenance != "ohnolog", "lineage"] = np.nan      # lineage only defined for ohnologs

    # marker flags
    pm = pd.read_csv(PBMC_MARK, sep="\t")
    pbmc_markers = set(pm.gene)
    cy = pd.read_csv(CYT_DEGS, sep="\t")
    cyt_resp = set(cy.human_symbol)
    u["is_pbmc_marker"] = u.gene.isin(pbmc_markers)
    u["is_cytokine_responder"] = u.gene.isin(cyt_resp)

    ohn = u[u.provenance == "ohnolog"]
    log(f"detectable paralogous TFs: {len(u)}; ohnolog-provenance: {len(ohn)}; "
        f"with alpha/beta label: {ohn.lineage.isin(['alpha','beta']).sum()} "
        f"({(ohn.lineage=='alpha').sum()} alpha / {(ohn.lineage=='beta').sum()} beta)")
    u[["gene", "dup_class4", "provenance", "arrangement", "TF_subfamily", "cluster_id",
       "lineage", "is_pbmc_marker", "is_cytokine_responder"]].to_csv(OUT_PER, sep="\t", index=False)

    is_alpha = (u.lineage == "alpha").values
    is_beta = (u.lineage == "beta").values
    is_ssd = (u.dup_class4 == "dispersed_SSD").values
    is_ohn = (u.provenance == "ohnolog").values

    rows = []
    for mod, mname in [("is_pbmc_marker", "PBMC_identity"),
                       ("is_cytokine_responder", "cytokine_response")]:
        mk = u[mod].values
        # rates
        for grp_mask, gname in [(is_alpha, "alpha"), (is_beta, "beta"),
                                (is_ohn, "all_ohnolog"), (is_ssd, "dispersed_SSD_ref")]:
            rows.append(dict(modality=mname, test=f"rate__{gname}", OR=np.nan, ci_lo=np.nan,
                             ci_hi=np.nan, p=np.nan, n=int(grp_mask.sum()),
                             rate=float(mk[grp_mask].mean()) if grp_mask.sum() else np.nan,
                             note=f"{int((mk&grp_mask).sum())}/{int(grp_mask.sum())}"))
        # headline alpha vs beta (matched provenance)
        r = contrast(is_alpha, is_beta, mk, "alpha_vs_beta"); r.update(modality=mname,
            n=int(is_alpha.sum() + is_beta.sum()), rate=np.nan); rows.append(r)
        # one-vs-rest against detectable paralogous-TF background
        r = contrast(is_alpha, ~is_alpha, mk, "alpha_vs_rest"); r.update(modality=mname,
            n=int(is_alpha.sum()), rate=np.nan); rows.append(r)
        r = contrast(is_beta, ~is_beta, mk, "beta_vs_rest"); r.update(modality=mname,
            n=int(is_beta.sum()), rate=np.nan); rows.append(r)

    sdf = pd.DataFrame(rows)
    sdf.to_csv(OUT_SUM, sep="\t", index=False)
    log("summary:")
    print(sdf[["modality", "test", "n", "rate", "OR", "ci_lo", "ci_hi", "p", "note"]]
          .to_string(index=False))
    make_figure(sdf)
    log(f"wrote {OUT_PER} / {OUT_SUM}")


def make_figure(sdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mods = [("PBMC_identity", "(A) PBMC steady-state identity markers"),
            ("cytokine_response", "(B) Cytokine-response DEGs")]
    COL = {"alpha": "#d6604d", "beta": "#4393c3", "all_ohnolog": "#878787",
           "dispersed_SSD_ref": "#d9d9d9"}
    bars = ["alpha", "beta", "all_ohnolog", "dispersed_SSD_ref"]
    blab = {"alpha": "alpha", "beta": "beta", "all_ohnolog": "all\nohnolog",
            "dispersed_SSD_ref": "dispersed\nSSD (ref)"}

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0),
                             gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for ax, (mod, title) in zip(axes[:2], mods):
        s = sdf[sdf.modality == mod].set_index("test")
        rates = [s.loc[f"rate__{b}", "rate"] for b in bars]
        ns = [int(s.loc[f"rate__{b}", "n"]) for b in bars]
        ax.bar(range(len(bars)), rates, color=[COL[b] for b in bars], width=0.7)
        for i, (rt, n) in enumerate(zip(rates, ns)):
            ax.text(i, rt + max(rates) * 0.02, f"{rt:.2f}\nn={n}", ha="center", fontsize=8.5)
        ab = s.loc["alpha_vs_beta"]
        ax.set_title(f"{title}\nalpha-vs-beta OR={ab.OR:.2f} "
                     f"(95% CI {ab.ci_lo:.2f}–{ab.ci_hi:.2f}), p={ab.p:.3f}", fontsize=9.5)
        ax.set_xticks(range(len(bars))); ax.set_xticklabels([blab[b] for b in bars], fontsize=8.5)
        ax.set_ylabel("fraction = marker / responder DEG")
        ax.set_ylim(0, max(rates) * 1.32)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    # forest of ORs
    ax = axes[2]
    keys = []
    for mod, mlab in [("PBMC_identity", "PBMC"), ("cytokine_response", "cytokine")]:
        for t, tl in [("alpha_vs_beta", "alpha vs beta"), ("alpha_vs_rest", "alpha vs rest"),
                      ("beta_vs_rest", "beta vs rest")]:
            keys.append((mod, t, f"{tl}\n[{mlab}]"))
    y = np.arange(len(keys))[::-1]
    s2 = sdf.set_index(["modality", "test"])
    orv = [s2.loc[(m, t), "OR"] for m, t, _ in keys]
    lo = [s2.loc[(m, t), "ci_lo"] for m, t, _ in keys]
    hi = [s2.loc[(m, t), "ci_hi"] for m, t, _ in keys]
    pv = [s2.loc[(m, t), "p"] for m, t, _ in keys]
    colors = ["#b2182b" if v >= 1 else "#2166ac" for v in orv]
    ax.hlines(y, lo, hi, color=colors, lw=1.8)
    ax.scatter(orv, y, c=colors, s=46, edgecolor="k", linewidth=0.5, zorder=3)
    ax.axvline(1, color="grey", ls="--", lw=1)
    ax.set_xscale("log", base=2)
    ax.set_yticks(y); ax.set_yticklabels([k[2] for k in keys], fontsize=8)
    for yi, v, h, p in zip(y, orv, hi, pv):
        ax.text(h * 1.05, yi, f"OR={v:.2f}{'*' if p < 0.05 else ''}", va="center", fontsize=7.5)
    ax.set_xlabel("odds ratio (log2)")
    ax.set_title("(C) Alpha/beta DEG-enrichment ORs", fontsize=9.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.suptitle("Alpha vs beta ohnolog DEG enrichment in human PBMC "
                 "(Zhu et al. 2026 allopolyploid lineages; Marletaz paralogons; Lambert TFs)",
                 fontsize=11.5, y=1.02)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT_FIG}.pdf/.png")


if __name__ == "__main__":
    main()
