#!/usr/bin/env python3
"""
09_integrate_ctcf_atac_e2g.py -- Integrate the three assays profiled over the SAME 77
clustered TF-paralog loci (47 clustered_cis_ohnolog / 30 clustered_SSD_tandem) into a
single headline figure + a combined cross-assay summary table.

Three complementary layers of regulatory architecture, all ohnolog vs SSD-tandem:
  1. CTCF  (insulation)   : CTCF peaks BETWEEN the two paralog members
                            [results/ctcf/cd4/ctcf_per_locus_summary.tsv]  (CD4 T cell, matched lineage)
  2. ATAC  (accessibility): distance from each member TSS to nearest ATAC peak
                            [results/ctcf/atac/ctcf_per_locus_summary.tsv] (Treg)
  3. E2G   (wiring)       : enhancer elements shared by BOTH paralogs
                            [results/e2g/e2g_pair_promoter_links.tsv]      (Treg)

Outputs:
  results/figures/integrated_three_assay_headline.{pdf,png}
  results/integrated_cross_assay_summary.tsv
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
G1, G2 = "clustered_cis_ohnolog", "clustered_SSD_tandem"
COL = {G1: "#0072B2", G2: "#D55E00"}
LAB = {G1: "cis-ohnolog", G2: "SSD-tandem"}


def load(path, key):
    out = {G1: [], G2: []}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            g = r["dup_class4"]
            if g in out and r.get(key) not in (None, "", "nan"):
                try:
                    out[g].append(float(r[key]))
                except ValueError:
                    pass
    return out


def mwu(a, b):
    a = [x for x in a if x == x]; b = [x for x in b if x == x]
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    rrb = 1.0 - (2.0 * U) / (len(a) * len(b))
    return U, p, rrb


def box_points(ax, data, ylab, title, log=False):
    vals = [data[G1], data[G2]]
    labels = [f"{LAB[G1]}\n(n={len(data[G1])})", f"{LAB[G2]}\n(n={len(data[G2])})"]
    colors = [COL[G1], COL[G2]]
    bp = ax.boxplot(vals, widths=0.55, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=1.4))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.35)
    rng = np.random.default_rng(0)
    for i, (v, c) in enumerate(zip(vals, colors), start=1):
        vv = np.array(v, float)
        if log:
            vv = np.clip(vv, 1.0, None)
        x = rng.normal(i, 0.06, size=len(vv))
        ax.scatter(x, vv, s=20, color=c, edgecolor="black", lw=0.3, zorder=3, alpha=0.85)
    if log:
        ax.set_yscale("log")
    ax.set_xticks([1, 2]); ax.set_xticklabels(labels)
    ax.set_ylabel(ylab); ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    U, p, rrb = mwu(vals[0], vals[1])
    ax.text(0.5, 0.97, f"MWU p={p:.2g}", transform=ax.transAxes, ha="center", va="top", fontsize=8.5)
    return U, p, rrb


def main():
    ctcf = load(ROOT / "results" / "ctcf" / "cd4" / "ctcf_per_locus_summary.tsv", "n_peaks_between_members")
    atac = load(ROOT / "results" / "ctcf" / "atac" / "ctcf_per_locus_summary.tsv", "nearest_peak_to_tss_mean")
    e2g = load(ROOT / "results" / "e2g" / "e2g_pair_promoter_links.tsv", "n_shared_enhancer_any")

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.4))
    r1 = box_points(axes[0], ctcf, "CTCF peaks between members",
                    "(a) CTCF insulation\n(CD4 T cell)")
    r2 = box_points(axes[1], atac, "nearest ATAC peak to TSS (bp)",
                    "(b) ATAC accessibility\n(Treg)", log=True)
    r3 = box_points(axes[2], e2g, "shared enhancers per pair",
                    "(c) E2G enhancer sharing\n(Treg)")
    fig.suptitle("Three assays, one locus set: WGD (ohnolog) paralogs are wired for independence; "
                 "tandem SSD paralogs for co-regulation", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(ROOT / "results" / "figures" / "integrated_three_assay_headline.pdf", bbox_inches="tight")
    fig.savefig(ROOT / "results" / "figures" / "integrated_three_assay_headline.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # combined summary table
    rows = [
        ("CTCF insulation (CD4 T)", "CTCF peaks between the two members", ctcf, r1, "ohnolog > SSD"),
        ("ATAC accessibility (Treg)", "Nearest ATAC peak to member TSS (bp, mean)", atac, r2, "SSD closer (more open)"),
        ("E2G wiring (Treg)", "Enhancers shared by both paralogs", e2g, r3, "SSD > ohnolog"),
    ]
    with open(ROOT / "results" / "integrated_cross_assay_summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["assay", "headline_metric", "median_ohnolog", "median_SSD",
                    "U", "p_value", "rank_biserial", "direction"])
        for assay, metric, data, (U, p, rrb), direction in rows:
            w.writerow([assay, metric, round(np.median(data[G1]), 3), round(np.median(data[G2]), 3),
                        U, f"{p:.4g}", round(rrb, 4), direction])

    for assay, metric, data, (U, p, rrb), direction in rows:
        print(f"{assay:26s} med_o={np.median(data[G1]):>9.2f} med_s={np.median(data[G2]):>9.2f} p={p:.4g} ({direction})")
    print("wrote results/figures/integrated_three_assay_headline.* and results/integrated_cross_assay_summary.tsv")


if __name__ == "__main__":
    main()
