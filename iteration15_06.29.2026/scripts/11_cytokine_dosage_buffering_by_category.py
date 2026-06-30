#!/usr/bin/env python3
"""
DOSAGE BUFFERING of paralogous TFs across the FOUR duplication categories (script 06),
in cytokine-response space. Dosage buffering != co-regulation: co-induction AMPLIFIES the
family total, buffering STABILISES it by letting members compensate. Raw correlation is
dominated by shared cytokine drive, so we use the two script-12 signatures, here per
within-family PAIR so each pair can be assigned a category:

  TEST 1  buffering index  B = Var(x+y) / (Var(x)+Var(y))  over union-active conditions.
            B < 1 -> sub-additive (members CANCEL = BUFFER); B > 1 -> super-additive
            (co-AMPLIFY); B = 1 -> independent.
  TEST 2  common-mode-removed PARTIAL correlation: regress each member's response on the
            common drive (mean response of all 4-category TFs per condition), correlate the
            residuals. partial_r < 0 = compensation beyond shared drive (the clean dosage
            signal).

Pairs (within-DBD-family pair construction): within-DBD-family, >=MIN_COND union conditions,
non-zero variance. Tagged by pair_category (shared dup_class4) and proximity
(NEIGHBOR=same cluster / DISTAL). Benchmarked against random CROSS-family pairs (null).

Comparisons: per-category median B / frac buffered (B<1) / median partial_r / frac
compensating (partial_r<0); Kruskal-Wallis across 4; proximity-matched cis-ohnolog vs
SSD-tandem (neighbors); exemplars (STAT1/STAT4, SAND).

Outputs:
  results/cytokine_dosage_buffering_by_category.pairs.tsv
  results/cytokine_dosage_buffering_by_category.summary.tsv
  results/figures/cytokine_dosage_buffering_by_category.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/11_cytokine_dosage_buffering_by_category.py
"""
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CSV = ROOT / "inputs" / "human_cytokine_dict" / "human_cytokine_dict_mini.csv"
CLASS = RES / "TF_dup_2x2_classification.tsv"
OUT_PAIRS = RES / "cytokine_dosage_buffering_by_category.pairs.tsv"
OUT_SUM = RES / "cytokine_dosage_buffering_by_category.summary.tsv"
OUT_FIG = RES / "figures" / "cytokine_dosage_buffering_by_category"

MIN_COND = 8
N_NULL = 3000
GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
COLORS = {"dispersed_ohnolog": "#2166ac", "clustered_cis_ohnolog": "#5aae61",
          "clustered_SSD_tandem": "#762a83", "dispersed_SSD": "#bababa"}
RNG = np.random.default_rng(0)


def log(*a):
    print("[11]", *a, flush=True)


def main():
    d = pd.read_csv(CSV, usecols=["gene", "log_fc", "celltype", "cytokine"])
    d["cond"] = d.celltype.astype(str) + "|" + d.cytokine.astype(str)
    gv = d.groupby(["gene", "cond"])["log_fc"].mean().unstack()
    genes_all = list(gv.index)
    gidx = {g: i for i, g in enumerate(genes_all)}
    A0 = gv.fillna(0.0).values
    MASK = gv.notna().values
    log(f"pseudobulk matrix: {A0.shape[0]} genes x {A0.shape[1]} conditions")

    cls = pd.read_csv(CLASS, sep="\t")
    present = set(gv.index)
    u = cls[cls.gene.isin(present)].copy()
    gene2cat = dict(zip(u.gene, u.dup_class4))
    gene2clu = dict(zip(u.gene, u.cluster_id))
    gene2fam = dict(zip(u.gene, u.TF_subfamily))
    log(f"4-group TFs with a profile: {len(u)}")

    # common mode = mean response of all 4-category TFs per condition (0-filled)
    cc_rows = [gidx[g] for g in u.gene]
    common = A0[cc_rows].mean(0)

    def metrics(a, b):
        ia, ib = gidx[a], gidx[b]
        m = MASK[ia] | MASK[ib]
        n = int(m.sum())
        if n < MIN_COND:
            return None
        x, y, cm = A0[ia][m], A0[ib][m], common[m]
        if np.std(x) == 0 or np.std(y) == 0:
            return None
        vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
        B = np.var(x + y, ddof=1) / (vx + vy) if (vx + vy) > 0 else np.nan
        mr = stats.pearsonr(x, y)[0]
        if np.std(cm) > 0:
            rx = x - np.polyval(np.polyfit(cm, x, 1), cm)
            ry = y - np.polyval(np.polyfit(cm, y, 1), cm)
            pr = stats.pearsonr(rx, ry)[0] if (np.std(rx) and np.std(ry)) else np.nan
        else:
            pr = mr
        return n, B, mr, pr

    # ---- within-family pairs ----
    rows = []
    for fam, sub in u.groupby("TF_subfamily"):
        gs = sorted(sub.gene)
        for a, b in combinations(gs, 2):
            mt = metrics(a, b)
            if mt is None:
                continue
            n, B, mr, pr = mt
            ca, cb = gene2cat[a], gene2cat[b]
            cla, clb = gene2clu.get(a), gene2clu.get(b)
            rows.append(dict(gene_a=a, gene_b=b, family=fam, n_conditions=n,
                             buffering_index_B=B, marginal_r=mr, partial_r=pr,
                             pair_category=(ca if ca == cb else "mixed"),
                             proximity=("NEIGHBOR" if (pd.notna(cla) and cla == clb) else "DISTAL")))
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT_PAIRS, sep="\t", index=False)
    log(f"scored within-family pairs: {len(pairs)}")

    # ---- null: random CROSS-family pairs ----
    fam_of = gene2fam
    glist = list(u.gene)
    null_B, null_pr = [], []
    tries = 0
    while len(null_B) < N_NULL and tries < N_NULL * 20:
        tries += 1
        a, b = RNG.choice(glist, 2, replace=False)
        if fam_of[a] == fam_of[b]:
            continue
        mt = metrics(a, b)
        if mt is None:
            continue
        _, B, _, pr = mt
        if np.isfinite(B):
            null_B.append(B)
        if np.isfinite(pr):
            null_pr.append(pr)
    null_B = np.array(null_B); null_pr = np.array(null_pr)
    log(f"null cross-family pairs: {len(null_B)} (median B={np.median(null_B):.3f}, "
        f"median partial_r={np.median(null_pr):.3f})")

    # ---- per-category summary ----
    summ = []
    for g in GROUPS:
        s = pairs[pairs.pair_category == g]
        n = len(s)
        # buffered fraction vs null; one-sided MWU that partial_r is LOWER (more compensating) than null
        pr_vals = s.partial_r.dropna().values
        try:
            _, p_pr = stats.mannwhitneyu(pr_vals, null_pr, alternative="less") if len(pr_vals) >= 3 else (np.nan, np.nan)
        except Exception:
            p_pr = np.nan
        try:
            _, p_B = stats.mannwhitneyu(s.buffering_index_B.dropna().values, null_B, alternative="two-sided") if n >= 3 else (np.nan, np.nan)
        except Exception:
            p_B = np.nan
        summ.append(dict(panel="per_category", group=g, n_pairs=n,
                         median_B=round(s.buffering_index_B.median(), 3) if n else np.nan,
                         frac_buffered_B_lt1=round((s.buffering_index_B < 1).mean(), 3) if n else np.nan,
                         median_marginal_r=round(s.marginal_r.median(), 3) if n else np.nan,
                         median_partial_r=round(s.partial_r.median(), 3) if n else np.nan,
                         frac_compensating_pr_lt0=round((s.partial_r < 0).mean(), 3) if n else np.nan,
                         p_partial_r_vs_null=p_pr, p_B_vs_null=p_B))
    summ.append(dict(panel="null", group="random_cross_family", n_pairs=len(null_B),
                     median_B=round(float(np.median(null_B)), 3), frac_buffered_B_lt1=round(float(np.mean(null_B < 1)), 3),
                     median_marginal_r=np.nan, median_partial_r=round(float(np.median(null_pr)), 3),
                     frac_compensating_pr_lt0=round(float(np.mean(null_pr < 0)), 3),
                     p_partial_r_vs_null=np.nan, p_B_vs_null=np.nan))

    # Kruskal-Wallis on B across the 4 categories
    arrs = [pairs.loc[pairs.pair_category == g, "buffering_index_B"].dropna().values for g in GROUPS]
    arrs = [a for a in arrs if len(a) >= 3]
    if len(arrs) >= 2:
        H, p = stats.kruskal(*arrs)
        summ.append(dict(panel="test", group="KruskalWallis_B_4cat", n_pairs=sum(len(a) for a in arrs),
                         median_B=np.nan, frac_buffered_B_lt1=np.nan, median_marginal_r=np.nan,
                         median_partial_r=np.nan, frac_compensating_pr_lt0=np.nan,
                         p_partial_r_vs_null=np.nan, p_B_vs_null=p))

    # proximity-matched neighbor contrast: cis-ohnolog vs SSD-tandem (B and partial_r)
    nb = pairs.proximity == "NEIGHBOR"
    cis = pairs.loc[nb & (pairs.pair_category == "clustered_cis_ohnolog")]
    sst = pairs.loc[nb & (pairs.pair_category == "clustered_SSD_tandem")]
    for metric in ["buffering_index_B", "partial_r"]:
        a = cis[metric].dropna().values; b = sst[metric].dropna().values
        p = stats.mannwhitneyu(a, b, alternative="two-sided")[1] if (len(a) >= 3 and len(b) >= 3) else np.nan
        summ.append(dict(panel="neighbor_cis_vs_SSD", group=metric, n_pairs=len(a) + len(b),
                         median_B=round(np.median(a), 3) if metric == "buffering_index_B" and len(a) else np.nan,
                         frac_buffered_B_lt1=np.nan,
                         median_marginal_r=np.nan,
                         median_partial_r=round(np.median(a), 3) if metric == "partial_r" and len(a) else np.nan,
                         frac_compensating_pr_lt0=np.nan, p_partial_r_vs_null=np.nan, p_B_vs_null=p,
                         ))
        summ[-1]["note"] = f"cis median={np.median(a):.2f} (n={len(a)}) vs SSD-tandem median={np.median(b):.2f} (n={len(b)})"

    sdf = pd.DataFrame(summ)
    sdf.to_csv(OUT_SUM, sep="\t", index=False)
    log("per-category dosage buffering:")
    print(sdf[sdf.panel.isin(["per_category", "null"])]
          [["group", "n_pairs", "median_B", "frac_buffered_B_lt1", "median_marginal_r",
            "median_partial_r", "frac_compensating_pr_lt0", "p_partial_r_vs_null"]].to_string(index=False))
    log("neighbor cis-ohnolog vs SSD-tandem:")
    print(sdf[sdf.panel == "neighbor_cis_vs_SSD"][["group", "note", "p_B_vs_null"]].to_string(index=False))

    # exemplars
    ex = pairs[((pairs.gene_a == "STAT1") & (pairs.gene_b == "STAT4")) |
               ((pairs.gene_a == "STAT4") & (pairs.gene_b == "STAT1")) |
               (pairs.gene_a.isin(["SP100", "SP110", "SP140"]) & pairs.gene_b.isin(["SP110", "SP140", "SP140L"]))]
    if len(ex):
        log("exemplar pairs:")
        print(ex[["gene_a", "gene_b", "pair_category", "buffering_index_B", "marginal_r", "partial_r"]].to_string(index=False))

    make_figure(pairs, null_B, null_pr)
    log(f"wrote {OUT_PAIRS} / {OUT_SUM}")


def make_figure(pairs, null_B, null_pr):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.4))

    # A: log2 B by category + null
    ax = axes[0]
    data = [np.log2(pairs.loc[pairs.pair_category == g, "buffering_index_B"].dropna().values) for g in GROUPS]
    data.append(np.log2(null_B[null_B > 0]))
    labels = ["disp\nohnolog", "clust\ncis-ohno", "clust\nSSD-tand", "disp\nSSD", "null\n(x-fam)"]
    cols = [COLORS[g] for g in GROUPS] + ["#ffffff"]
    bp = ax.boxplot(data, vert=True, widths=0.6, patch_artist=True, showfliers=False)
    for box, c in zip(bp["boxes"], cols):
        box.set(facecolor=c, alpha=0.75)
    for i, arr in enumerate(data, 1):
        if len(arr):
            ax.text(i, ax.get_ylim()[1] if False else np.percentile(arr, 95) + 0.1,
                    f"n={len(arr)}", ha="center", fontsize=7)
    ax.axhline(0, color="grey", ls="--", lw=1)
    ax.set_xticks(range(1, 6)); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("log2 buffering index B  (<0 = buffer, >0 = amplify)")
    ax.set_title("(A) Buffering index by category", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # B: partial_r by category + null
    ax = axes[1]
    data = [pairs.loc[pairs.pair_category == g, "partial_r"].dropna().values for g in GROUPS]
    data.append(null_pr)
    bp = ax.boxplot(data, vert=True, widths=0.6, patch_artist=True, showfliers=False)
    for box, c in zip(bp["boxes"], cols):
        box.set(facecolor=c, alpha=0.75)
    ax.axhline(0, color="grey", ls="--", lw=1)
    ax.set_xticks(range(1, 6)); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("common-mode-removed partial r  (<0 = compensation)")
    ax.set_title("(B) Compensation (partial r) by category", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # C: marginal -> partial r shift, neighbors of the two clustered categories
    ax = axes[2]
    nb = pairs.proximity == "NEIGHBOR"
    for g, c in [("clustered_cis_ohnolog", "#5aae61"), ("clustered_SSD_tandem", "#762a83")]:
        s = pairs.loc[nb & (pairs.pair_category == g)]
        for _, r in s.iterrows():
            ax.plot([0, 1], [r.marginal_r, r.partial_r], color=c, alpha=0.5, lw=0.9)
        ax.scatter([0] * len(s), s.marginal_r, color=c, s=18, zorder=3)
        ax.scatter([1] * len(s), s.partial_r, color=c, s=18, zorder=3,
                   label=f"{g.replace('clustered_','')} (n={len(s)})")
    ax.axhline(0, color="grey", ls="--", lw=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["marginal r", "partial r\n(common-mode removed)"], fontsize=8)
    ax.set_ylabel("pair correlation")
    ax.set_title("(C) Neighbor pairs: shared-drive vs compensation", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.suptitle("Cytokine-response DOSAGE BUFFERING across the 4 duplication categories "
                 "(Oesinghaus pseudobulk; B<1 / partial r<0 = buffering)", fontsize=11.5, y=1.02)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT_FIG}.pdf/.png")


if __name__ == "__main__":
    main()
