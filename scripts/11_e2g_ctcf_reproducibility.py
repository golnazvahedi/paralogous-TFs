#!/usr/bin/env python3
"""
11_e2g_ctcf_reproducibility.py -- Aggregate the per-cell-type E2G-CTCF group comparisons
(script 10 outputs) across the 6-cell lymphoid panel into one table + figure.

Panel (ATAC-based ENCODE-rE2G, extended model, hg38; CTCF = per-element feature col 14):
  Treg (ENCFF393GIF) | memory CD4 (ENCFF440AKC) | CD8 (ENCFF993UOE) |
  Jurkat (ENCFF705FWU) | NK (ENCFF359ZSL) | B cell (ENCFF186ABX)

Question: does CTCF signal at the paralog-pair regulatory elements differ between
clustered_cis_ohnolog and clustered_SSD_tandem, and is any difference reproducible?

Outputs:
  results/e2g/e2g_ctcf_reproducibility_summary.tsv
  results/figures/e2g/e2g_ctcf_reproducibility.{pdf,png}
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
C_OHNO, C_SSD = "#0072B2", "#D55E00"

PANEL = [
    ("Treg",       "",          "ENCFF393GIF"),
    ("memory CD4", "memoryCD4", "ENCFF440AKC"),
    ("CD8",        "CD8",       "ENCFF993UOE"),
    ("Jurkat",     "Jurkat",    "ENCFF705FWU"),
    ("NK",         "NK",        "ENCFF359ZSL"),
    ("B cell",     "Bcell",     "ENCFF186ABX"),
]
KEY = ["ctcf_mean_pair_enh", "ctcf_mean_pair_prom", "ctcf_mean_shared_enh", "ctcf_mean_window"]
NAMES = {"ctcf_mean_pair_enh": "pair-linked enhancers",
         "ctcf_mean_pair_prom": "pair promoters",
         "ctcf_mean_shared_enh": "shared enhancers",
         "ctcf_mean_window": "+/-100 kb window"}


def load_comp(subdir):
    d = ROOT / "results" / "e2g" / subdir if subdir else ROOT / "results" / "e2g"
    rows = {}
    with open(d / "e2g_ctcf_group_comparison.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows[r["metric"]] = r
    return rows


def main():
    summary = []
    for label, subdir, acc in PANEL:
        comp = load_comp(subdir)
        rec = dict(cell_type=label, accession=acc)
        for k in KEY:
            r = comp[k]
            rec[f"{k}__med_ohno"] = float(r["median_ohnolog"]) if r["median_ohnolog"] not in ("nan", "") else float("nan")
            rec[f"{k}__med_SSD"] = float(r["median_SSD"]) if r["median_SSD"] not in ("nan", "") else float("nan")
            rec[f"{k}__p"] = float(r["p_value"]) if r["p_value"] not in ("nan", "") else float("nan")
            rec[f"{k}__rrb"] = float(r["rank_biserial"]) if r["rank_biserial"] not in ("nan", "") else float("nan")
        summary.append(rec)

    out = ROOT / "results" / "e2g"
    with open(out / "e2g_ctcf_reproducibility_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(summary)

    labels = [s["cell_type"] for s in summary]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    # panel A: median CTCF at pair-linked enhancers, ohnolog vs SSD
    ax = axes[0]
    med_o = [s["ctcf_mean_pair_enh__med_ohno"] for s in summary]
    med_s = [s["ctcf_mean_pair_enh__med_SSD"] for s in summary]
    ax.bar(x - 0.19, med_o, 0.38, color=C_OHNO, label="cis-ohnolog", edgecolor="black", lw=0.4)
    ax.bar(x + 0.19, med_s, 0.38, color=C_SSD, label="SSD-tandem", edgecolor="black", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("median CTCF signal at pair-linked enhancers")
    ax.set_title("(a) CTCF at pair-linked enhancers", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    # panel B: -log10 p across cell types for each CTCF metric
    ax = axes[1]
    marks = {"ctcf_mean_pair_enh": "o", "ctcf_mean_pair_prom": "s",
             "ctcf_mean_shared_enh": "D", "ctcf_mean_window": "^"}
    for k in KEY:
        y = [-np.log10(max(s[f"{k}__p"], 1e-300)) if s[f"{k}__p"] == s[f"{k}__p"] else np.nan for s in summary]
        ax.plot(x, y, marks[k] + "-", ms=6, lw=1.0, label=NAMES[k])
    ax.axhline(-np.log10(0.05), color="0.5", ls="--", lw=0.8)
    ax.text(len(labels) - 1, -np.log10(0.05) + 0.03, "p=0.05", fontsize=7, color="0.4", ha="right", va="bottom")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel(r"$-\log_{10}\, p$ (MWU, ohnolog vs SSD)")
    ax.set_title("(b) Significance across the panel", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("E2G CTCF at paralog-pair elements: no reproducible ohnolog-vs-SSD difference across 6 lymphoid cell types",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(out.parent / "figures" / "e2g" / "e2g_ctcf_reproducibility.pdf", bbox_inches="tight")
    fig.savefig(out.parent / "figures" / "e2g" / "e2g_ctcf_reproducibility.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"{'cell_type':12s} {'enh_ohno':>9} {'enh_SSD':>8} {'p(enh)':>9} {'rrb':>7}")
    for s in summary:
        print(f"{s['cell_type']:12s} {s['ctcf_mean_pair_enh__med_ohno']:>9.3f} "
              f"{s['ctcf_mean_pair_enh__med_SSD']:>8.3f} {s['ctcf_mean_pair_enh__p']:>9.3g} "
              f"{s['ctcf_mean_pair_enh__rrb']:>7.3f}")
    print(f"wrote {out}/e2g_ctcf_reproducibility_summary.tsv and figures/e2g/e2g_ctcf_reproducibility.*")


if __name__ == "__main__":
    main()
