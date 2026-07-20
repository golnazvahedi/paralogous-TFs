#!/usr/bin/env python3
"""
05_paired_promoters_ctcf.py -- Restrict to PAIRED clusters (exactly 2 members) and
report, per pair:
  - the two promoter (TSS) coordinates and orientation
  - promoter-promoter distance  = |TSS_A - TSS_B|  (bp)
  - number of CTCF peaks between the two promoters (summits strictly inside [minTSS,maxTSS))
  - CTCF peaks/kb between promoters (distance-normalized)

Also compares the two clustered classes (clustered_cis_ohnolog vs clustered_SSD_tandem)
on promoter distance and between-promoter CTCF (Mann-Whitney U, two-sided).

Inputs : results/clustered_loci_windows.pm100kb.tsv  (has member_tss, member_strands)
         inputs/encode/ENCFF796WRU.bed.gz            (CTCF peaks)
Outputs: results/ctcf/paired_promoters_ctcf.tsv
         results/ctcf/paired_promoters_group_comparison.tsv
         results/figures/paired_promoters_ctcf.{pdf,png}
"""
from pathlib import Path
import csv, gzip
import numpy as np
from bisect import bisect_left
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

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
            start, end = int(p[1]), int(p[2])
            off = int(p[9]) if len(p) > 9 and p[9] not in (".", "-1") else (end - start) // 2
            if off < 0:
                off = (end - start) // 2
            summ[p[0]].append(start + off)
    return {c: sorted(v) for c, v in summ.items()}


def count_between(summits, chrom, a, b):
    if chrom not in summits:
        return 0
    s = summits[chrom]
    lo, hi = min(a, b), max(a, b)
    return bisect_left(s, hi) - bisect_left(s, lo)


def orientation(strands, tss):
    """Classify the 2-gene arrangement from strands + TSS order.

    Convention: '+' gene transcribes toward higher coords (TSS/head at its low end),
    '-' gene transcribes toward lower coords (TSS/head at its high end). Order the two
    genes by TSS coordinate into left (l) and right (r):
      - same strand                 -> tandem     (-> -> or <- <-)
      - left '-' , right '+'        -> head-to-head (divergent: 5' heads face each other,
                                       transcription points away, promoters in the middle)
      - left '+' , right '-'        -> tail-to-tail (convergent: 3' tails meet in the middle,
                                       transcription points toward each other)
    """
    if strands[0] == strands[1]:
        return "tandem"
    (t_l, s_l), (t_r, s_r) = sorted(zip(tss, strands))
    if s_l == "-" and s_r == "+":
        return "head-to-head"
    if s_l == "+" and s_r == "-":
        return "tail-to-tail"
    return "tandem"


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
    rows = []
    with open(ROOT / "results" / "clustered_loci_windows.pm100kb.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["n_members"] != "2":
                continue
            tss = [int(x) for x in r["member_tss"].split(",")]
            strands = r["member_strands"].split(",")
            members = r["members"].split(",")
            chrom = r["chrom"]
            dist = abs(tss[0] - tss[1])
            n_between = count_between(summits, chrom, tss[0], tss[1])
            per_kb = round(n_between / (dist / 1000), 5) if dist > 0 else ""
            rows.append({
                "cluster_id": r["cluster_id"],
                "dup_class4": r["dup_class4"],
                "family": r["family"],
                "chrom": chrom,
                "geneA": members[0], "geneB": members[1],
                "strandA": strands[0], "strandB": strands[1],
                "tssA": tss[0], "tssB": tss[1],
                "orientation": orientation(strands, tss),
                "promoter_distance_bp": dist,
                "n_ctcf_between_promoters": n_between,
                "ctcf_between_per_kb": per_kb,
            })

    rows.sort(key=lambda x: (x["dup_class4"], x["promoter_distance_bp"]))
    cols = list(rows[0].keys())
    outp = ctcf_dir / "paired_promoters_ctcf.tsv"
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    # ---- group comparison ----
    def vals(g, k):
        return np.array([float(r[k]) for r in rows if r["dup_class4"] == g and r[k] != ""], float)

    comp = []
    for key, label in [("promoter_distance_bp", "Promoter-promoter distance (bp)"),
                       ("n_ctcf_between_promoters", "CTCF peaks between promoters (count)"),
                       ("ctcf_between_per_kb", "CTCF peaks/kb between promoters")]:
        a, b = vals(G1, key), vals(G2, key)
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rb = 1.0 - (2.0 * u) / (len(a) * len(b))
        comp.append({
            "metric": key, "label": label,
            "n_ohnolog": len(a), "n_SSD": len(b),
            "median_ohnolog": round(float(np.median(a)), 4),
            "median_SSD": round(float(np.median(b)), 4),
            "U": round(float(u), 1), "p_value": f"{p:.4g}",
            "rank_biserial_r": round(float(rb), 4),
            "test": "Mann-Whitney U (two-sided), clustered_cis_ohnolog vs clustered_SSD_tandem; 2-member clusters only",
        })
    cp = ctcf_dir / "paired_promoters_group_comparison.tsv"
    with open(cp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comp[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(comp)

    # ---- figure ----
    rng = np.random.RandomState(20260718)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3))
    panels = [("promoter_distance_bp", "Promoter-promoter distance", "bp", True),
              ("n_ctcf_between_promoters", "CTCF peaks between promoters", "count", False),
              ("ctcf_between_per_kb", "CTCF peaks/kb between promoters", "peaks/kb", False)]
    for ax, (key, ttl, ylab, logy) in zip(axes, panels):
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
        if logy:
            ax.set_yscale("log")
        ax.set_title(ttl, fontsize=10)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"{LAB[G1]}\n(n={len(a)})", f"{LAB[G2]}\n(n={len(b)})"], fontsize=9)
        ax.text(0.5, 0.98, f"MWU p={row['p_value']}\nr_rb={row['rank_biserial_r']}",
                transform=ax.transAxes, ha="center", va="top", fontsize=8.5, color="0.2")
    fig.suptitle("Paired (2-member) clustered TF loci: promoter distance and inter-promoter CTCF\n"
                 f"({CELL}; Mann-Whitney U, two-sided)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(fd / "paired_promoters_ctcf.pdf")
    fig.savefig(fd / "paired_promoters_ctcf.png", dpi=200)
    plt.close(fig)

    # ---- console ----
    from collections import Counter
    cc = Counter(r["dup_class4"] for r in rows)
    print(f"paired (2-member) clusters: {dict(cc)}")
    print(f"\norientation breakdown:")
    oc = Counter((r["dup_class4"], r["orientation"]) for r in rows)
    for k in sorted(oc):
        print(f"  {k[0]:24s} {k[1]:14s} {oc[k]}")
    print("\n=== group comparison (2-member clusters) ===")
    for c in comp:
        print(f"{c['metric']:28s} med_ohno={c['median_ohnolog']:>10} "
              f"med_SSD={c['median_SSD']:>10} U={c['U']:>7} p={c['p_value']:>9} r_rb={c['rank_biserial_r']}")
    print(f"\nwrote {outp}\n      {cp}\n      {fd}/paired_promoters_ctcf.png")


if __name__ == "__main__":
    main()
