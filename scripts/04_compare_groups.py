#!/usr/bin/env python3
"""
04_compare_groups.py -- Compare CTCF architecture between the two CLUSTERED groups:
  clustered_cis_ohnolog  vs  clustered_SSD_tandem.

Tests (named explicitly in outputs per project reproducibility rule):
  * Mann-Whitney U (two-sided) per per-locus continuous metric.
    Report U, p, n per group, medians, rank-biserial effect size.
  * Per-bin metagene MWU (density_per_kb) with BH-FDR across bins to localize
    where along the window the two groups differ.
  * Mean CTCF metagene profiles per group with bootstrap 95% CI bands.

Figures:
  results/figures/ctcf_metagene_clustered_groups.{pdf,png}
  results/figures/ctcf_peaks_per_kb_by_group.{pdf,png}
  results/figures/ctcf_perlocus_heatmap.{pdf,png}   (per-locus binned CTCF)

Tables:
  results/ctcf/ctcf_group_comparison.tsv       (per-locus metric MWU)
  results/ctcf/ctcf_metagene_perbin_mwu.tsv     (per-bin MWU + BH-FDR)
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
G1 = "clustered_cis_ohnolog"
G2 = "clustered_SSD_tandem"
# colorblind-safe
C1 = "#0072B2"   # ohnolog (blue)
C2 = "#D55E00"   # SSD tandem (vermillion)
COL = {G1: C1, G2: C2}
LAB = {G1: "clustered cis-ohnolog", G2: "clustered SSD-tandem"}

METRICS = [
    ("peaks_per_kb_window", "CTCF peaks/kb (+/-100 kb window)"),
    ("peaks_per_kb_body", "CTCF peaks/kb (locus body)"),
    ("peaks_per_kb_flank", "CTCF peaks/kb (flanks)"),
    ("obs_exp_window", "Observed/expected CTCF (window)"),
    ("nearest_peak_to_tss_min", "Nearest CTCF peak to a member TSS (bp, min)"),
    ("nearest_peak_to_tss_mean", "Nearest CTCF peak to member TSS (bp, mean)"),
    ("n_peaks_between_members", "CTCF peaks between paralog members"),
    ("inter_member_tss_dist_bp", "Inter-member TSS distance (bp)"),
    ("peaks_between_per_kb", "CTCF peaks/kb between members (dist-normalized)"),
]


def rank_biserial(u, n1, n2):
    """Rank-biserial effect size from Mann-Whitney U (group1 orientation)."""
    return 1.0 - (2.0 * u) / (n1 * n2)


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def load_summary(ctcf_dir):
    rows = []
    with open(ctcf_dir / "ctcf_per_locus_summary.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    return rows


def load_metagene(ctcf_dir):
    rows = []
    with open(ctcf_dir / "ctcf_metagene.matrix.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    return rows


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--subdir", default="", help="read/write under results/ctcf/<subdir> and results/figures/<subdir>")
    ap.add_argument("--label", default="GM12878 CTCF, ENCFF796WRU", help="cell-type label for figure titles")
    args = ap.parse_args()
    ctcf_dir = ROOT / "results" / "ctcf" / args.subdir if args.subdir else ROOT / "results" / "ctcf"
    fig_dir = ROOT / "results" / "figures" / args.subdir if args.subdir else ROOT / "results" / "figures"
    ctcf_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    CELL = args.label

    summ = load_summary(ctcf_dir)
    g1 = [r for r in summ if r["dup_class4"] == G1]
    g2 = [r for r in summ if r["dup_class4"] == G2]
    n1, n2 = len(g1), len(g2)

    # ---------- per-locus metric MWU ----------
    comp_rows = []
    for key, label in METRICS:
        a = np.array([float(r[key]) for r in g1 if r[key] != ""], float)
        b = np.array([float(r[key]) for r in g2 if r[key] != ""], float)
        if len(a) < 3 or len(b) < 3:
            continue
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rb = rank_biserial(u, len(a), len(b))
        comp_rows.append({
            "metric": key, "label": label,
            "n_ohnolog": len(a), "n_SSD": len(b),
            "median_ohnolog": round(float(np.median(a)), 5),
            "median_SSD": round(float(np.median(b)), 5),
            "mean_ohnolog": round(float(np.mean(a)), 5),
            "mean_SSD": round(float(np.mean(b)), 5),
            "U": round(float(u), 2), "p_value": f"{p:.4g}",
            "rank_biserial_r": round(float(rb), 4),
            "test": "Mann-Whitney U (two-sided), clustered_cis_ohnolog vs clustered_SSD_tandem",
        })
    cp = ctcf_dir / "ctcf_group_comparison.tsv"
    with open(cp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comp_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(comp_rows)

    # ---------- metagene per-group profile + bootstrap CI ----------
    meta = load_metagene(ctcf_dir)
    # index -> rel_pos, region
    idx_meta = {}
    for r in meta:
        i = int(r["meta_index"])
        idx_meta.setdefault(i, (float(r["rel_pos"]), r["region"]))
    indices = sorted(idx_meta)

    def group_matrix(group):
        # per-locus vector of density_per_kb across bins
        by_locus = {}
        for r in meta:
            if r["dup_class4"] != group:
                continue
            by_locus.setdefault(r["cluster_id"], {})[int(r["meta_index"])] = float(r["density_per_kb"])
        mat = []
        for cid, d in by_locus.items():
            mat.append([d.get(i, 0.0) for i in indices])
        return np.array(mat, float)

    M1 = group_matrix(G1)
    M2 = group_matrix(G2)

    rng = np.random.RandomState(20260718)

    def boot_ci(mat, nboot=2000):
        means = mat.mean(axis=0)
        n = mat.shape[0]
        boot = np.empty((nboot, mat.shape[1]))
        for b in range(nboot):
            idx = rng.randint(0, n, n)
            boot[b] = mat[idx].mean(axis=0)
        lo = np.percentile(boot, 2.5, axis=0)
        hi = np.percentile(boot, 97.5, axis=0)
        return means, lo, hi

    m1, lo1, hi1 = boot_ci(M1)
    m2, lo2, hi2 = boot_ci(M2)

    rel = np.array([idx_meta[i][0] for i in indices])

    # ---------- per-bin MWU + BH-FDR ----------
    perbin_rows = []
    pvals = []
    for j, i in enumerate(indices):
        a = M1[:, j]; b = M2[:, j]
        try:
            u, p = mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            u, p = np.nan, 1.0
        pvals.append(p)
        perbin_rows.append({
            "meta_index": i, "rel_pos": idx_meta[i][0], "region": idx_meta[i][1],
            "mean_ohnolog": round(float(a.mean()), 5), "mean_SSD": round(float(b.mean()), 5),
            "U": round(float(u), 2) if not np.isnan(u) else "", "p_value": f"{p:.4g}",
        })
    q = bh_fdr(pvals)
    for r, qq in zip(perbin_rows, q):
        r["q_value_BH"] = f"{qq:.4g}"
    pbp = ctcf_dir / "ctcf_metagene_perbin_mwu.tsv"
    with open(pbp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(perbin_rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(perbin_rows)

    # fig_dir already set from args at top of main()

    # ---------- FIG 1: metagene ----------
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for m, lo, hi, g in [(m1, lo1, hi1, G1), (m2, lo2, hi2, G2)]:
        ax.plot(rel, m, color=COL[g], lw=2, label=f"{LAB[g]} (n={M1.shape[0] if g==G1 else M2.shape[0]})")
        ax.fill_between(rel, lo, hi, color=COL[g], alpha=0.18, linewidth=0)
    ax.axvspan(0, 1, color="0.85", alpha=0.4, zorder=0)
    ax.axvline(0, color="0.4", ls=":", lw=0.9)
    ax.axvline(1, color="0.4", ls=":", lw=0.9)
    # mark FDR-significant bins
    sig = [rel[j] for j in range(len(indices)) if float(q[j]) < 0.05]
    ymax = max(hi1.max(), hi2.max())
    for x in sig:
        ax.plot(x, ymax * 1.04, marker="v", color="black", ms=5)
    ax.set_xticks([-1, 0, 1, 2])
    ax.set_xticklabels(["-100 kb", "locus start", "locus end", "+100 kb"])
    ax.set_ylabel("CTCF peak density (summits per kb)")
    ax.set_xlabel("relative position across extended clustered locus")
    ax.set_title("CTCF metagene profile: clustered ohnolog vs SSD-tandem loci\n"
                 f"({CELL}; shaded = locus body; bands = bootstrap 95% CI)",
                 fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    if sig:
        ax.text(0.99, 0.02, "▼ per-bin MWU BH-FDR<0.05", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.5, color="0.3")
    fig.tight_layout()
    fig.savefig(fig_dir / "ctcf_metagene_clustered_groups.pdf")
    fig.savefig(fig_dir / "ctcf_metagene_clustered_groups.png", dpi=200)
    plt.close(fig)

    # ---------- FIG 2: per-locus peaks/kb by group (box+points) ----------
    panel_metrics = ["peaks_per_kb_window", "peaks_per_kb_body", "peaks_per_kb_flank", "obs_exp_window"]
    fig, axes = plt.subplots(1, 4, figsize=(13, 4.2))
    for ax, key in zip(axes, panel_metrics):
        a = np.array([float(r[key]) for r in g1 if r[key] != ""], float)
        b = np.array([float(r[key]) for r in g2 if r[key] != ""], float)
        data = [a, b]
        bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.4))
        for patch, g in zip(bp["boxes"], [G1, G2]):
            patch.set_facecolor(COL[g]); patch.set_alpha(0.35)
        for xi, (arr, g) in enumerate(zip(data, [G1, G2]), start=1):
            jitter = (rng.rand(len(arr)) - 0.5) * 0.28
            ax.scatter(np.full(len(arr), xi) + jitter, arr, s=16, color=COL[g],
                       alpha=0.8, edgecolor="white", linewidth=0.4, zorder=3)
        row = next((c for c in comp_rows if c["metric"] == key), None)
        pv = row["p_value"] if row else "NA"
        rb = row["rank_biserial_r"] if row else "NA"
        ax.set_title(dict(METRICS)[key], fontsize=9)
        ax.set_xticks([1, 2]); ax.set_xticklabels(["ohnolog", "SSD-tandem"], fontsize=8.5)
        ax.text(0.5, 0.97, f"MWU p={pv}\nr_rb={rb}", transform=ax.transAxes,
                ha="center", va="top", fontsize=8, color="0.25")
        if key == "obs_exp_window":
            ax.axhline(1.0, color="0.5", ls="--", lw=0.8)
    axes[0].set_ylabel("value")
    fig.suptitle("Per-locus CTCF metrics by clustered duplication class "
                 f"(ohnolog n={n1}, SSD-tandem n={n2}; Mann-Whitney U, two-sided)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(fig_dir / "ctcf_peaks_per_kb_by_group.pdf")
    fig.savefig(fig_dir / "ctcf_peaks_per_kb_by_group.png", dpi=200)
    plt.close(fig)

    # ---------- FIG 3: per-locus heatmap ----------
    order_ids = [r["cluster_id"] for r in g1] + [r["cluster_id"] for r in g2]
    dens = {}
    for r in meta:
        dens.setdefault(r["cluster_id"], {})[int(r["meta_index"])] = float(r["density_per_kb"])
    H = np.array([[dens[cid].get(i, 0.0) for i in indices] for cid in order_ids], float)
    fig, ax = plt.subplots(figsize=(9, 10))
    vmax = np.percentile(H, 98)
    im = ax.imshow(H, aspect="auto", cmap="magma", vmin=0, vmax=vmax, interpolation="nearest")
    ax.axhline(n1 - 0.5, color="cyan", lw=1.6)
    # x ticks at rel positions -1,0,1,2
    xt = [np.argmin(np.abs(rel - v)) for v in [-1, 0, 1, 2]]
    ax.set_xticks(xt); ax.set_xticklabels(["-100kb", "start", "end", "+100kb"], fontsize=9)
    ax.set_yticks([n1 / 2, n1 + n2 / 2])
    ax.set_yticklabels([f"cis-ohnolog\n(n={n1})", f"SSD-tandem\n(n={n2})"], fontsize=10)
    for j in [0, 1]:  # body boundaries
        ax.axvline(xt[1] - 0.5, color="white", ls=":", lw=0.7)
        ax.axvline(xt[2] + 0.5, color="white", ls=":", lw=0.7)
    ax.set_title("Per-locus CTCF density (metagene bins), rows split by class", fontsize=11)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("CTCF summits per kb", fontsize=9)
    fig.tight_layout()
    fig.savefig(fig_dir / "ctcf_perlocus_heatmap.pdf")
    fig.savefig(fig_dir / "ctcf_perlocus_heatmap.png", dpi=200)
    plt.close(fig)

    # ---------- FIG 4: CTCF between paralog members (headline metric) ----------
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.2))
    for ax, key, ttl in zip(
            axes,
            ["n_peaks_between_members", "peaks_between_per_kb"],
            ["CTCF peaks between\nparalog members (raw count)",
             "CTCF peaks/kb between members\n(inter-member-distance normalized)"]):
        a = np.array([float(r[key]) for r in g1 if r[key] != ""], float)
        b = np.array([float(r[key]) for r in g2 if r[key] != ""], float)
        bp = ax.boxplot([a, b], widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.4))
        for patch, g in zip(bp["boxes"], [G1, G2]):
            patch.set_facecolor(COL[g]); patch.set_alpha(0.35)
        for xi, (arr, g) in enumerate(zip([a, b], [G1, G2]), start=1):
            jitter = (rng.rand(len(arr)) - 0.5) * 0.28
            ax.scatter(np.full(len(arr), xi) + jitter, arr, s=18, color=COL[g],
                       alpha=0.8, edgecolor="white", linewidth=0.4, zorder=3)
        row = next((c for c in comp_rows if c["metric"] == key), None)
        ax.set_title(ttl, fontsize=9.5)
        ax.set_xticks([1, 2]); ax.set_xticklabels(["ohnolog", "SSD-tandem"], fontsize=9)
        ax.text(0.5, 0.97, f"MWU p={row['p_value']}\nr_rb={row['rank_biserial_r']}",
                transform=ax.transAxes, ha="center", va="top", fontsize=8.5, color="0.2")
    axes[0].set_ylabel("count / density")
    fig.suptitle("Insulator separation of paralog members\n"
                 f"(ohnolog n={n1}, SSD-tandem n={n2}; {CELL})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(fig_dir / "ctcf_between_members.pdf")
    fig.savefig(fig_dir / "ctcf_between_members.png", dpi=200)
    plt.close(fig)

    # ---------- console summary ----------
    print("=== Per-locus metric comparison (Mann-Whitney U, two-sided) ===")
    print(f"{'metric':28s} {'med_ohno':>10s} {'med_SSD':>10s} {'U':>8s} {'p':>10s} {'r_rb':>7s}")
    for r in comp_rows:
        print(f"{r['metric']:28s} {r['median_ohnolog']:>10} {r['median_SSD']:>10} "
              f"{r['U']:>8} {r['p_value']:>10} {r['rank_biserial_r']:>7}")
    nsig = sum(1 for qq in q if qq < 0.05)
    print(f"\nper-bin metagene MWU: {nsig}/{len(indices)} bins BH-FDR<0.05")
    print(f"wrote {cp}\n      {pbp}\n      figures in {fig_dir}")


if __name__ == "__main__":
    main()
