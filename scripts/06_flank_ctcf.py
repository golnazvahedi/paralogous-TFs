#!/usr/bin/env python3
"""
06_flank_ctcf.py -- CTCF enrichment in the two 100 kb FLANKS of each clustered locus,
reported separately by group.

Regions per locus (genomic coordinates; boundaries = outermost member TSS/TES = the
locus body edges cluster_start / cluster_end):
  - LEFT  (upstream)   flank: [win_start, cluster_start]   ~ "-100 kb -> locus start"
  - RIGHT (downstream) flank: [cluster_end, win_end]        ~ "locus end -> +100 kb"
(win_start/win_end are already clipped to chromosome bounds.)

Per flank: CTCF summit count, length (kb), peaks/kb, and observed/expected vs the
genome-wide baseline (44,217 summits / 3.03 Gb = 0.01459 peaks/kb).

Comparisons (Mann-Whitney U, two-sided):
  * clustered_cis_ohnolog vs clustered_SSD_tandem, for the LEFT flank and RIGHT flank.
  * combined-flank peaks/kb by group (matches iteration19 CLAUDE.md caveat wording).
Also within-group LEFT vs RIGHT (Wilcoxon signed-rank, paired) to test flank asymmetry.

Outputs:
  results/ctcf/flank_ctcf_per_locus.tsv
  results/ctcf/flank_ctcf_group_comparison.tsv
  results/figures/flank_ctcf_by_group.{pdf,png}
"""
from pathlib import Path
import csv, gzip
import numpy as np
from bisect import bisect_left
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
G1, G2 = "clustered_cis_ohnolog", "clustered_SSD_tandem"
COL = {G1: "#0072B2", G2: "#D55E00"}
LAB = {G1: "cis-ohnolog", G2: "SSD-tandem"}


def load_summits(peaks_path):
    summ = defaultdict(list)
    with gzip.open(peaks_path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] not in PRIMARY:
                continue
            s, e = int(p[1]), int(p[2])
            off = int(p[9]) if len(p) > 9 and p[9] not in (".", "-1") else (e - s) // 2
            if off < 0:
                off = (e - s) // 2
            summ[p[0]].append(s + off)
    return {c: sorted(v) for c, v in summ.items()}


def n_summits(summits, chrom, a, b):
    if chrom not in summits or b <= a:
        return 0
    s = summits[chrom]
    return bisect_left(s, b) - bisect_left(s, a)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--peaks", default=str(ROOT / "inputs" / "encode" / "ENCFF796WRU.bed.gz"))
    ap.add_argument("--subdir", default="", help="write outputs under results/ctcf/<subdir> and results/figures/<subdir>")
    ap.add_argument("--label", default="GM12878 CTCF, ENCFF796WRU", help="cell-type label for figure title")
    args = ap.parse_args()
    ctcf_dir = ROOT / "results" / "ctcf" / args.subdir if args.subdir else ROOT / "results" / "ctcf"
    fd = ROOT / "results" / "figures" / args.subdir if args.subdir else ROOT / "results" / "figures"
    ctcf_dir.mkdir(parents=True, exist_ok=True)
    fd.mkdir(parents=True, exist_ok=True)
    CELL = args.label

    summits = load_summits(args.peaks)
    total = sum(len(v) for v in summits.values())
    sizes = {}
    with open(ROOT / "inputs" / "hg38.chrom.sizes") as f:
        for line in f:
            c, s = line.rstrip("\n").split("\t")[:2]
            sizes[c] = int(s)
    total_bp = sum(sizes[c] for c in summits if c in sizes)
    exp_rate_kb = total / total_bp * 1000  # peaks per kb genome-wide

    rows = []
    with open(ROOT / "results" / "clustered_loci_windows.pm100kb.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            chrom = r["chrom"]
            ws, we = int(r["win_start"]), int(r["win_end"])
            cs, ce = int(r["cluster_start"]), int(r["cluster_end"])
            lstart, lend = ws, cs            # left flank
            rstart, rend = ce, we            # right flank
            llen = lend - lstart
            rlen = rend - rstart
            ln = n_summits(summits, chrom, lstart, lend)
            rn = n_summits(summits, chrom, rstart, rend)
            l_per_kb = ln / (llen / 1000) if llen > 0 else 0
            r_per_kb = rn / (rlen / 1000) if rlen > 0 else 0
            rows.append({
                "cluster_id": r["cluster_id"], "dup_class4": r["dup_class4"],
                "family": r["family"], "chrom": chrom, "n_members": r["n_members"],
                "members": r["members"], "clipped": r["clipped"],
                "left_flank_len_kb": round(llen / 1000, 2),
                "left_ctcf_peaks": ln,
                "left_peaks_per_kb": round(l_per_kb, 5),
                "left_obs_exp": round(l_per_kb / exp_rate_kb, 4) if exp_rate_kb else 0,
                "right_flank_len_kb": round(rlen / 1000, 2),
                "right_ctcf_peaks": rn,
                "right_peaks_per_kb": round(r_per_kb, 5),
                "right_obs_exp": round(r_per_kb / exp_rate_kb, 4) if exp_rate_kb else 0,
            })

    outp = ctcf_dir / "flank_ctcf_per_locus.tsv"
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    def vals(g, k):
        return np.array([float(r[k]) for r in rows if r["dup_class4"] == g], float)

    # ---- between-group comparisons per flank ----
    comp = []
    for key, label in [("left_peaks_per_kb", "LEFT flank CTCF peaks/kb (-100 kb -> locus start)"),
                       ("left_obs_exp", "LEFT flank observed/expected"),
                       ("right_peaks_per_kb", "RIGHT flank CTCF peaks/kb (locus end -> +100 kb)"),
                       ("right_obs_exp", "RIGHT flank observed/expected")]:
        a, b = vals(G1, key), vals(G2, key)
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rb = 1.0 - (2.0 * u) / (len(a) * len(b))
        comp.append({
            "comparison": "ohnolog_vs_SSD", "metric": key, "label": label,
            "n_ohnolog": len(a), "n_SSD": len(b),
            "median_ohnolog": round(float(np.median(a)), 4),
            "median_SSD": round(float(np.median(b)), 4),
            "U": round(float(u), 1), "p_value": f"{p:.4g}",
            "rank_biserial_r": round(float(rb), 4),
            "test": "Mann-Whitney U (two-sided), clustered_cis_ohnolog vs clustered_SSD_tandem",
        })

    # ---- within-group LEFT vs RIGHT (paired Wilcoxon signed-rank) ----
    for g in (G1, G2):
        l = vals(g, "left_peaks_per_kb")
        rr = vals(g, "right_peaks_per_kb")
        try:
            w_stat, p = wilcoxon(l, rr)
        except ValueError:
            w_stat, p = float("nan"), 1.0
        comp.append({
            "comparison": f"{LAB[g]}_left_vs_right", "metric": "left_vs_right_peaks_per_kb",
            "label": f"{LAB[g]}: LEFT vs RIGHT flank peaks/kb (paired)",
            "n_ohnolog": len(l) if g == G1 else "", "n_SSD": len(rr) if g == G2 else "",
            "median_ohnolog": round(float(np.median(l)), 4),
            "median_SSD": round(float(np.median(rr)), 4),
            "U": round(float(w_stat), 1) if not np.isnan(w_stat) else "",
            "p_value": f"{p:.4g}", "rank_biserial_r": "",
            "test": "Wilcoxon signed-rank (paired), LEFT vs RIGHT flank within group",
        })

    cp = ctcf_dir / "flank_ctcf_group_comparison.tsv"
    with open(cp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comp[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(comp)

    # ---- figure ----
    rng = np.random.RandomState(20260718)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
    panels = [("left_peaks_per_kb", "LEFT flank\n(-100 kb -> locus start)"),
              ("right_peaks_per_kb", "RIGHT flank\n(locus end -> +100 kb)")]
    for ax, (key, ttl) in zip(axes, panels):
        a, b = vals(G1, key), vals(G2, key)
        bp = ax.boxplot([a, b], widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.4))
        for patch, g in zip(bp["boxes"], [G1, G2]):
            patch.set_facecolor(COL[g]); patch.set_alpha(0.35)
        for xi, (arr, g) in enumerate(zip([a, b], [G1, G2]), start=1):
            j = (rng.rand(len(arr)) - 0.5) * 0.28
            ax.scatter(np.full(len(arr), xi) + j, arr, s=20, color=COL[g],
                       alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
        row = next(c for c in comp if c["metric"] == key)
        ax.axhline(exp_rate_kb, color="0.5", ls="--", lw=0.9)
        ax.set_title(ttl, fontsize=10)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"{LAB[G1]}\n(n={len(a)})", f"{LAB[G2]}\n(n={len(b)})"], fontsize=9)
        ax.text(0.5, 0.98, f"MWU p={row['p_value']}\nr_rb={row['rank_biserial_r']}",
                transform=ax.transAxes, ha="center", va="top", fontsize=8.5, color="0.2")
    axes[0].set_ylabel("CTCF peaks/kb")
    axes[1].text(0.98, exp_rate_kb, " genome-wide expected", transform=axes[1].get_yaxis_transform(),
                 ha="right", va="bottom", fontsize=7.5, color="0.4")
    fig.suptitle("Per-flank CTCF enrichment by clustered duplication class\n"
                 f"({CELL}; dashed = genome-wide expected rate)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(fd / "flank_ctcf_by_group.pdf")
    fig.savefig(fd / "flank_ctcf_by_group.png", dpi=200)
    plt.close(fig)

    # ---- console ----
    print(f"genome-wide expected CTCF rate: {exp_rate_kb:.5f} peaks/kb\n")
    print("=== group medians (peaks/kb | obs/exp) ===")
    for g in (G1, G2):
        print(f"{LAB[g]:12s} (n={int((np.array([r['dup_class4'] for r in rows])==g).sum())}): "
              f"LEFT {np.median(vals(g,'left_peaks_per_kb')):.4f} ({np.median(vals(g,'left_obs_exp')):.2f}x)  "
              f"RIGHT {np.median(vals(g,'right_peaks_per_kb')):.4f} ({np.median(vals(g,'right_obs_exp')):.2f}x)")
    print("\n=== comparisons ===")
    for c in comp:
        print(f"{c['comparison']:26s} {c['metric']:26s} "
              f"med_o={c['median_ohnolog']:>8} med_s={c['median_SSD']:>8} "
              f"p={c['p_value']:>9} r_rb={c['rank_biserial_r']}")
    print(f"\nwrote {outp}\n      {cp}\n      {fd}/flank_ctcf_by_group.png")


if __name__ == "__main__":
    main()
