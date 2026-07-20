#!/usr/bin/env python3
"""
08_e2g_reproducibility.py -- Aggregate the per-cell-type ENCODE-rE2G group comparisons
(script 07 outputs) into one reproducibility table + figure, testing whether the
clustered_cis_ohnolog vs clustered_SSD_tandem enhancer-sharing signature reproduces
across a panel of lymphoid cell types.

Panel (all ATAC-based ENCODE-rE2G, extended model, hg38):
  Treg        (ENCFF393GIF)  results/e2g/                 [T cell]
  memory CD4  (ENCFF440AKC)  results/e2g/memoryCD4/       [T cell]
  CD8         (ENCFF993UOE)  results/e2g/CD8/             [T cell]
  Jurkat      (ENCFF705FWU)  results/e2g/Jurkat/          [T cell (T-ALL line)]
  NK          (ENCFF359ZSL)  results/e2g/NK/              [lymphoid control]
  B cell      (ENCFF186ABX)  results/e2g/Bcell/           [lymphoid control]

Key metrics carried across (Mann-Whitney U, two-sided, ohnolog vs SSD-tandem):
  obs_exp_window, n_distinct_elem_total, n_shared_enhancer_any, shared_frac

Outputs:
  results/e2g/e2g_reproducibility_summary.tsv
  results/figures/e2g/e2g_reproducibility.{pdf,png}
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
C_OHNO, C_SSD = "#0072B2", "#D55E00"

# (label, subdir, accession, celltype_class)
PANEL = [
    ("Treg",       "",          "ENCFF393GIF", "T cell"),
    ("memory CD4", "memoryCD4", "ENCFF440AKC", "T cell"),
    ("CD8",        "CD8",       "ENCFF993UOE", "T cell"),
    ("Jurkat",     "Jurkat",    "ENCFF705FWU", "T cell (line)"),
    ("NK",         "NK",        "ENCFF359ZSL", "lymphoid ctrl"),
    ("B cell",     "Bcell",     "ENCFF186ABX", "lymphoid ctrl"),
]
KEY = ["obs_exp_window", "n_distinct_elem_total", "n_shared_enhancer_any", "shared_frac"]


def load_comp(subdir):
    d = ROOT / "results" / "e2g" / subdir if subdir else ROOT / "results" / "e2g"
    rows = {}
    with open(d / "e2g_group_comparison.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[r["metric"]] = r
    return rows


def main():
    summary = []
    for label, subdir, acc, ctype in PANEL:
        comp = load_comp(subdir)
        rec = dict(cell_type=label, accession=acc, celltype_class=ctype)
        for k in KEY:
            r = comp[k]
            rec[f"{k}__med_ohno"] = float(r["median_ohnolog"])
            rec[f"{k}__med_SSD"] = float(r["median_SSD"])
            rec[f"{k}__p"] = float(r["p_value"])
            rec[f"{k}__rrb"] = float(r["rank_biserial"])
        summary.append(rec)

    out = ROOT / "results" / "e2g"
    with open(out / "e2g_reproducibility_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(summary)

    labels = [s["cell_type"] for s in summary]
    x = np.arange(len(labels))

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))

    # --- panel A: shared enhancers, median ohno vs SSD, grouped bars ---
    ax = axes[0]
    med_o = [s["n_shared_enhancer_any__med_ohno"] for s in summary]
    med_s = [s["n_shared_enhancer_any__med_SSD"] for s in summary]
    ax.bar(x - 0.19, med_o, 0.38, color=C_OHNO, label="cis-ohnolog", edgecolor="black", lw=0.4)
    ax.bar(x + 0.19, med_s, 0.38, color=C_SSD, label="SSD-tandem", edgecolor="black", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("median shared enhancers / pair")
    ax.set_title("(a) Enhancers shared by both paralogs", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    # --- panel B: -log10 p across cell types for each metric ---
    ax = axes[1]
    marks = {"obs_exp_window": "o", "n_distinct_elem_total": "s",
             "n_shared_enhancer_any": "D", "shared_frac": "^"}
    names = {"obs_exp_window": "obs/exp density", "n_distinct_elem_total": "distinct enhancers",
             "n_shared_enhancer_any": "shared enhancers", "shared_frac": "shared fraction"}
    for k in KEY:
        y = [-np.log10(max(s[f"{k}__p"], 1e-300)) for s in summary]
        ax.plot(x, y, marks[k] + "-", ms=6, lw=1.0, label=names[k])
    ax.axhline(-np.log10(0.05), color="0.5", ls="--", lw=0.8)
    ax.text(len(labels) - 1, -np.log10(0.05) + 0.05, "p=0.05", fontsize=7, color="0.4", ha="right", va="bottom")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel(r"$-\log_{10}\, p$ (MWU, ohnolog vs SSD)")
    ax.set_title("(b) Significance across the panel", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    # --- panel C: rank-biserial effect size (shared enhancers) forest ---
    ax = axes[2]
    rrb = [s["n_shared_enhancer_any__rrb"] for s in summary]
    colors = [("#222" if s["celltype_class"].startswith("T") else "#888") for s in summary]
    ax.barh(x, rrb, color=colors, edgecolor="black", lw=0.4)
    ax.set_yticks(x); ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_xlabel("rank-biserial effect size\n(shared enhancers, SSD > ohnolog)")
    ax.set_title("(c) Effect size (shared enhancers)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=[Patch(color="#222", label="T cell"), Patch(color="#888", label="lymphoid control")],
              frameon=False, fontsize=8, loc="lower right")

    fig.suptitle("ENCODE-rE2G reproducibility across 6 lymphoid cell types: SSD-tandem paralogs share enhancers, ohnologs do not",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(ROOT / "results" / "figures" / "e2g" / "e2g_reproducibility.pdf", bbox_inches="tight")
    fig.savefig(ROOT / "results" / "figures" / "e2g" / "e2g_reproducibility.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # console
    print(f"{'cell_type':12s} {'shared_ohno':>11} {'shared_SSD':>10} {'p(shared)':>11} {'rrb':>6}")
    for s in summary:
        print(f"{s['cell_type']:12s} {s['n_shared_enhancer_any__med_ohno']:>11} "
              f"{s['n_shared_enhancer_any__med_SSD']:>10} {s['n_shared_enhancer_any__p']:>11.3g} "
              f"{s['n_shared_enhancer_any__rrb']:>6}")
    print(f"wrote {out}/e2g_reproducibility_summary.tsv and figures/e2g/e2g_reproducibility.*")


if __name__ == "__main__":
    main()
