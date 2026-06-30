#!/usr/bin/env python3
"""
Standalone OR forest plots for the z-score-DEG enrichment OR results (script 07b).

07b embeds the ORs in a 3-panel figure; this renders a clean, publication-style forest
plot dedicated to the OR values, read straight from the 07b summary TSV.

Default target = the expression-filtered z>=2 run (TAG="zscore_z2_e1"). Override with
  `python3 07d_plot_OR_zscore.py [TAG]`   e.g. zscore_z2  (no expr filter) or zscore_z1_e1.

Reads  : results/TF_dup4_DEG_OR.<TAG>.summary.tsv
Writes : results/figures/TF_dup4_DEG_OR_forest.<TAG>.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/07d_plot_OR_zscore.py [TAG]
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
TAG = sys.argv[1] if len(sys.argv) > 1 else "zscore_z2_e1"
SUM = RES / f"TF_dup4_DEG_OR.{TAG}.summary.tsv"
OUT = RES / "figures" / f"TF_dup4_DEG_OR_forest.{TAG}"

GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
GCOLOR = {"dispersed_ohnolog": "#2166ac", "clustered_cis_ohnolog": "#5aae61",
          "clustered_SSD_tandem": "#762a83", "dispersed_SSD": "#bababa"}


def log(*a):
    print("[07d]", *a, flush=True)


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else ""


CLIP_LO, CLIP_HI = 2.0 ** -5, 2.0 ** 5    # OR/CI beyond this are "unstable" (quasi-separation)


def forest(ax, rows, title, color_by_group=False):
    """rows: list of (label, OR, lo, hi, p, fdr_or_None, color). Drawn top->bottom.

    ORs/CIs beyond [CLIP_LO, CLIP_HI] (e.g. logit quasi-separation when ~every gene is a
    DEG) are clamped to the axis edge, drawn with an arrow marker, and flagged 'unstable'
    so a degenerate threshold renders honestly instead of exploding the axis."""
    y = np.arange(len(rows))[::-1]
    for (lab, orr, lo, hi, p, fdr, col), yi in zip(rows, y):
        unstable = not (CLIP_LO <= orr <= CLIP_HI) or not np.isfinite(hi) or hi > CLIP_HI * 4
        cl = lambda v: min(max(v, CLIP_LO), CLIP_HI) if np.isfinite(v) else CLIP_HI
        clo, chi, corr = cl(lo), cl(hi), cl(orr)
        ax.hlines(yi, clo, chi, color=col, lw=2.0, zorder=2)
        ax.plot([clo, clo], [yi - 0.14, yi + 0.14], color=col, lw=1.4)
        ax.plot([chi, chi], [yi - 0.14, yi + 0.14], color=col, lw=1.4)
        ax.scatter([corr], [yi], color=col, s=70, edgecolor="k", linewidth=0.6, zorder=3,
                   marker=(">" if orr > CLIP_HI else "<" if orr < CLIP_LO else "o"))
        sig = stars(fdr) if fdr is not None and np.isfinite(fdr) else stars(p)
        tag = f"  OR={orr:.2f}{sig}" + ("  (unstable)" if unstable else "")
        ax.text(chi * 1.06, yi, tag, va="center", fontsize=8.5)
    ax.axvline(1.0, color="grey", ls="--", lw=1)
    ax.set_xscale("log", base=2)
    ax.set_xlim(CLIP_LO / 1.3, CLIP_HI * 2.4)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("odds ratio (log2 scale)")
    ax.set_title(title, fontsize=10.5, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    if not SUM.exists():
        sys.exit(f"missing {SUM} -- run 07b for TAG={TAG} first "
                 f"(e.g. scripts/07b_human_pbmc_dup4_DEG_OR_zscore.py 2 1)")
    s = pd.read_csv(SUM, sep="\t").set_index("test")
    log(f"loaded {SUM}")

    # ---- panel A: per-category one-vs-rest ----
    a_rows = []
    for g in GROUPS:
        r = s.loc[f"onevsrest__{g}"]
        note = str(r["note"])
        rate = note.split("DEG_rate=")[1].split(" ")[0] if "DEG_rate=" in note else "?"
        n = note.split("n=")[1].split(",")[0] if "n=" in note else "?"
        lab = f"{g.replace('_', ' ')}\n(n={n}, DEG rate={rate})"
        a_rows.append((lab, r.OR, r.ci_lo, r.ci_hi, r.p, r.fdr, GCOLOR[g]))

    # ---- panel B: key contrasts + factorial ----
    cdefs = [
        ("ohnolog vs SSD  (within dispersed)", "pair__within_dispersed__ohnolog_vs_SSD"),
        ("cis-ohnolog vs SSD-tandem  (within clustered)", "pair__within_clustered__cis_ohnolog_vs_SSD_tandem"),
        ("clustered vs dispersed  (within ohnolog)", "pair__within_ohnolog__clustered_vs_dispersed"),
        ("clustered vs dispersed  (within SSD)", "pair__within_SSD__clustered_vs_dispersed"),
        ("ohnolog vs SSD  (marginal)", "marginal__ohnolog_vs_SSD"),
        ("clustered vs dispersed  (marginal)", "marginal__clustered_vs_dispersed"),
        ("ohnolog  (factorial logit, adj.)", "logit__ohnolog"),
        ("clustered  (factorial logit, adj.)", "logit__clustered"),
        ("ohnolog x clustered  (interaction)", "logit__ohnolog:clustered"),
    ]
    b_rows = []
    for lab, key in cdefs:
        if key not in s.index:
            continue
        r = s.loc[key]
        col = "#2166ac" if r.OR >= 1 else "#b2182b"
        b_rows.append((lab, r.OR, r.ci_lo, r.ci_hi, r.p, r.get("fdr", np.nan), col))

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.4), gridspec_kw={"width_ratios": [1, 1.25]})
    forest(axes[0], a_rows,
           "(A) Per-category enrichment vs rest\n(one-vs-rest Fisher OR; * by BH-FDR)")
    forest(axes[1], b_rows,
           "(B) Key contrasts & provenance x arrangement factorial\n(* by raw p)")

    # human-readable tag: "zscore_z2" -> "z>2, no expr filter"; "zscore_z2_e1" -> "z>2, expr>1 CP10K-sum filter"
    zpart = "z>1" if "_z1" in TAG else "z>2" if "_z2" in TAG else TAG
    if "_e" in TAG:
        ethr = TAG.split("_e")[1]
        tagnote = f"{zpart}, expr>{ethr} CP10K-sum filter"
    else:
        tagnote = f"{zpart}, no expr filter"
    fig.suptitle("PBMC z-score-DEG enrichment OR across the 4 duplication categories "
                 f"[{tagnote}]", fontsize=12, y=1.03)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT}.pdf / .png")


if __name__ == "__main__":
    main()
