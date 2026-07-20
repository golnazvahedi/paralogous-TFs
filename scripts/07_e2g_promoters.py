#!/usr/bin/env python3
"""
07_e2g_promoters.py -- ENCODE-rE2G (ATAC-based enhancer->gene) analysis of the two
CLUSTERED paralog groups, focused on the GENE-PAIR PROMOTERS.

Data: ENCODE-rE2G thresholded element-gene links for regulatory CD4 T cells (Treg),
      annotation ENCSR759OVP, file ENCFF393GIF (bed3+; multiomic ATAC+H3K27ac+HiC+RNA).
Each row = a regulatory ELEMENT (chr,start,end) predicted to regulate TargetGene, with
an ABC/rE2G Score, DistanceToTSS, per-element ATAC accessibility (`Open`), and class
(promoter/lncRNA/intergenic/genic).

Column layout (0-indexed) of ENCFF393GIF:
  0 chr  1 start  2 end  3 name  4 class  5 TargetGene  6 EnsemblID  7 TargetGeneTSS
  8 CellType  9 Score  10 DistanceToTSS  11 H3K27ac  12 Open(ATAC)  13 H3K4me  14 CTCF ...

Two analyses (only the two clustered classes, per iteration19 CLAUDE.md):
  (1) PER-LOCUS E2G ELEMENT DENSITY (accessibility mirror of the CTCF density analysis):
      unique elements deduped by (chr,start,end); element midpoints counted in the
      +/-100 kb window / locus body / flanks; elements/kb and observed/expected vs a
      genome-wide baseline (n_unique_elements / total primary bp).
  (2) GENE-PAIR PROMOTER REGULATORY WIRING (the core question):
      for each cluster, collect E2G links whose TargetGene is a cluster member. Report,
      per member: n links and n distinct enhancer elements landing on that promoter, and
      the mean element ATAC accessibility. Per cluster: the number of enhancer elements
      SHARED between two (or more) paralog members -- a shared-regulatory-input / co-
      regulation proxy -- plus a window-restricted shared count (elements inside +/-100 kb).

Group comparison: Mann-Whitney U (two-sided), clustered_cis_ohnolog vs clustered_SSD_tandem,
reported with U, p, n per group, medians, and rank-biserial effect size.

Outputs:
  results/e2g/e2g_per_locus_summary.tsv
  results/e2g/e2g_pair_promoter_links.tsv
  results/e2g/e2g_group_comparison.tsv
  results/figures/e2g/e2g_element_density_by_group.{pdf,png}
  results/figures/e2g/e2g_shared_enhancers_by_group.{pdf,png}
  results/figures/e2g/e2g_promoter_links_by_group.{pdf,png}
  log/07_e2g_promoters.log
"""
from pathlib import Path
import csv, gzip
from bisect import bisect_left, bisect_right
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
G1, G2 = "clustered_cis_ohnolog", "clustered_SSD_tandem"
COL = {G1: "#0072B2", G2: "#D55E00"}
LAB = {G1: "cis-ohnolog", G2: "SSD-tandem"}
PAD = 100_000
E2G = ROOT / "inputs" / "encode" / "ENCFF393GIF.bed.gz"       # default: Treg
CELL = "Treg ENCODE-rE2G, ENCFF393GIF"


def rank_biserial(u, n1, n2):
    return 1.0 - (2.0 * u) / (n1 * n2) if n1 and n2 else float("nan")


def mwu(a, b):
    a = [x for x in a if x == x]
    b = [x for x in b if x == x]
    if len(a) < 3 or len(b) < 3:
        return dict(U=float("nan"), p=float("nan"), n1=len(a), n2=len(b),
                    med1=np.median(a) if a else float("nan"),
                    med2=np.median(b) if b else float("nan"), r_rb=float("nan"))
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    return dict(U=U, p=p, n1=len(a), n2=len(b),
                med1=float(np.median(a)), med2=float(np.median(b)),
                r_rb=rank_biserial(U, len(a), len(b)))


def load_e2g(path):
    """Return:
       elem_mid[chrom] -> sorted list of unique element midpoints (for density)
       gene_links[gene] -> list of (chrom, estart, eend, score, open_atac, klass, dist)
       n_unique_elements, total considered links
    """
    elem_seen = defaultdict(set)          # chrom -> set of (start,end)
    gene_links = defaultdict(list)
    n_links = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            chrom = p[0]
            if chrom not in PRIMARY:
                continue
            start, end = int(p[1]), int(p[2])
            klass = p[4]
            gene = p[5]
            try:
                score = float(p[9])
            except ValueError:
                score = float("nan")
            try:
                openv = float(p[12])
            except (ValueError, IndexError):
                openv = float("nan")
            try:
                dist = int(p[10])
            except (ValueError, IndexError):
                dist = -1
            elem_seen[chrom].add((start, end))
            gene_links[gene].append((chrom, start, end, score, openv, klass, dist))
            n_links += 1
    elem_mid = {}
    n_unique = 0
    for c, s in elem_seen.items():
        mids = sorted((a + b) // 2 for (a, b) in s)
        elem_mid[c] = mids
        n_unique += len(mids)
    return elem_mid, gene_links, n_unique, n_links


def count_mid(mids, chrom, a, b):
    if chrom not in mids or b <= a:
        return 0
    m = mids[chrom]
    return bisect_left(m, b) - bisect_left(m, a)


def main():
    import argparse
    global E2G, CELL
    ap = argparse.ArgumentParser()
    ap.add_argument("--e2g", default=str(E2G), help="ENCODE-rE2G thresholded element-gene links bed.gz")
    ap.add_argument("--subdir", default="", help="write outputs under results/e2g/<subdir> and results/figures/e2g/<subdir>")
    ap.add_argument("--label", default=CELL, help="cell-type label for figure titles")
    args = ap.parse_args()
    E2G = Path(args.e2g)
    CELL = args.label
    sub = args.subdir

    # ---- load windows (the 77 clustered loci) ----
    loci = []
    with open(ROOT / "results" / "clustered_loci_windows.pm100kb.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["dup_class4"] not in (G1, G2):
                continue
            loci.append(row)

    elem_mid, gene_links, n_unique, n_links = load_e2g(E2G)
    total_bp = 0
    # genome-wide baseline element rate over primary chroms present in the file
    for c, m in elem_mid.items():
        pass
    # use hg38 chrom sizes for baseline denominator (primary chroms in the file)
    sizes = {}
    with open(ROOT / "inputs" / "hg38.chrom.sizes") as f:
        for line in f:
            cc, ss = line.split()[:2]
            if cc in PRIMARY:
                sizes[cc] = int(ss)
    total_bp = sum(sizes[c] for c in elem_mid if c in sizes)
    base_rate_kb = n_unique / (total_bp / 1000.0)   # elements per kb genome-wide

    # ---------- (1) per-locus element density ----------
    perloc = []
    for r in loci:
        chrom = r["chrom"]
        ws, we = int(r["win_start"]), int(r["win_end"])
        cs, ce = int(r["cluster_start"]), int(r["cluster_end"])
        win_kb = (we - ws) / 1000.0
        body_kb = max(ce - cs, 1) / 1000.0
        flank_kb = ((cs - ws) + (we - ce)) / 1000.0
        n_win = count_mid(elem_mid, chrom, ws, we)
        n_body = count_mid(elem_mid, chrom, cs, ce)
        n_flank = n_win - n_body
        perloc.append(dict(
            cluster_id=r["cluster_id"], dup_class4=r["dup_class4"], chrom=chrom,
            n_members=r["n_members"], members=r["members"],
            elem_win=n_win, elem_body=n_body, elem_flank=n_flank,
            elem_per_kb_window=n_win / win_kb,
            elem_per_kb_body=n_body / body_kb,
            elem_per_kb_flank=(n_flank / flank_kb) if flank_kb > 0 else float("nan"),
            obs_exp_window=(n_win / win_kb) / base_rate_kb,
        ))

    # ---------- (2) gene-pair promoter regulatory wiring ----------
    pair_rows = []
    for r in loci:
        chrom = r["chrom"]
        ws, we = int(r["win_start"]), int(r["win_end"])
        members = [m for m in r["members"].split(",") if m]
        # links per member; distinct elements per member; accessibility
        per_member = {}
        elem_to_members = defaultdict(set)   # (chrom,s,e) -> set(members) with a link
        elem_open = {}
        for gene in members:
            lk = gene_links.get(gene, [])
            # promoter-class self element(s) vs distal enhancer elements
            dist_elems = set()
            prom_elems = set()
            opens = []
            for (c, s, e, sc, ov, klass, dist) in lk:
                key = (c, s, e)
                elem_to_members[key].add(gene)
                if ov == ov:
                    elem_open[key] = ov
                if klass == "promoter":
                    prom_elems.add(key)
                else:
                    dist_elems.add(key)
                    if ov == ov:
                        opens.append(ov)
            per_member[gene] = dict(
                n_links=len(lk),
                n_distinct_elem=len({(c, s, e) for (c, s, e, *_ ) in lk}),
                n_enhancer_elem=len(dist_elems),
                n_promoter_elem=len(prom_elems),
                mean_open_enh=float(np.mean(opens)) if opens else float("nan"),
            )
        # shared enhancer elements: linked to >=2 members (any-distance), and window-restricted
        shared_any = 0
        shared_in_win = 0
        for key, gset in elem_to_members.items():
            if len(gset) >= 2:
                shared_any += 1
                c, s, e = key
                mid = (s + e) // 2
                if c == chrom and ws <= mid < we:
                    shared_in_win += 1
        total_distinct = len(elem_to_members)
        shared_frac = (shared_any / total_distinct) if total_distinct else float("nan")

        # aggregate across members (report the pair as a unit; also min/max)
        n_links_tot = sum(per_member[m]["n_links"] for m in members)
        enh_counts = [per_member[m]["n_enhancer_elem"] for m in members]
        opens_all = [per_member[m]["mean_open_enh"] for m in members if per_member[m]["mean_open_enh"] == per_member[m]["mean_open_enh"]]
        pair_rows.append(dict(
            cluster_id=r["cluster_id"], dup_class4=r["dup_class4"], chrom=chrom,
            n_members=len(members), members=r["members"],
            n_links_total=n_links_tot,
            n_distinct_elem_total=total_distinct,
            n_enhancer_elem_min=min(enh_counts) if enh_counts else 0,
            n_enhancer_elem_max=max(enh_counts) if enh_counts else 0,
            n_enhancer_elem_sum=sum(enh_counts),
            n_shared_enhancer_any=shared_any,
            n_shared_enhancer_in_window=shared_in_win,
            shared_frac=shared_frac,
            mean_open_enh=float(np.mean(opens_all)) if opens_all else float("nan"),
            per_member=";".join(f"{m}:links={per_member[m]['n_links']},enh={per_member[m]['n_enhancer_elem']}" for m in members),
        ))

    # ---------- write tables ----------
    out = ROOT / "results" / "e2g" / sub if sub else ROOT / "results" / "e2g"
    out.mkdir(parents=True, exist_ok=True)
    figd = ROOT / "results" / "figures" / "e2g" / sub if sub else ROOT / "results" / "figures" / "e2g"
    figd.mkdir(parents=True, exist_ok=True)

    with open(out / "e2g_per_locus_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(perloc[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(perloc)
    with open(out / "e2g_pair_promoter_links.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pair_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)

    # ---------- group comparison ----------
    def by_group(rows, key):
        g1 = [float(r[key]) for r in rows if r["dup_class4"] == G1 and r[key] == r[key]]
        g2 = [float(r[key]) for r in rows if r["dup_class4"] == G2 and r[key] == r[key]]
        return g1, g2

    comp = []
    metrics = [
        ("per_locus", "elem_per_kb_window", perloc, "E2G elements/kb (+/-100 kb window)"),
        ("per_locus", "elem_per_kb_body", perloc, "E2G elements/kb (locus body)"),
        ("per_locus", "elem_per_kb_flank", perloc, "E2G elements/kb (flanks)"),
        ("per_locus", "obs_exp_window", perloc, "E2G observed/expected (window)"),
        ("pair", "n_links_total", pair_rows, "E2G links to pair promoters (total)"),
        ("pair", "n_distinct_elem_total", pair_rows, "Distinct enhancer elements on pair promoters"),
        ("pair", "n_enhancer_elem_sum", pair_rows, "Enhancer elements summed over members"),
        ("pair", "n_shared_enhancer_any", pair_rows, "Shared enhancer elements (both paralogs)"),
        ("pair", "n_shared_enhancer_in_window", pair_rows, "Shared enhancers within +/-100 kb"),
        ("pair", "shared_frac", pair_rows, "Shared-enhancer fraction"),
        ("pair", "mean_open_enh", pair_rows, "Mean enhancer ATAC accessibility (Open)"),
    ]
    for scope, key, rows, label in metrics:
        g1, g2 = by_group(rows, key)
        m = mwu(g1, g2)
        comp.append(dict(scope=scope, metric=key, label=label,
                         median_ohnolog=round(m["med1"], 5), median_SSD=round(m["med2"], 5),
                         n_ohnolog=m["n1"], n_SSD=m["n2"],
                         U=m["U"], p_value=m["p"], rank_biserial=round(m["r_rb"], 4) if m["r_rb"] == m["r_rb"] else float("nan")))
    with open(out / "e2g_group_comparison.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comp[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(comp)

    # ---------- figures ----------
    def box_points(ax, key, rows, ylab, title):
        data, labels, colors = [], [], []
        for g in (G1, G2):
            vals = [float(r[key]) for r in rows if r["dup_class4"] == g and r[key] == r[key]]
            data.append(vals); labels.append(f"{LAB[g]}\n(n={len(vals)})"); colors.append(COL[g])
        bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color="black", lw=1.4))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.35)
        rng = np.random.default_rng(0)
        for i, (vals, c) in enumerate(zip(data, colors), start=1):
            x = rng.normal(i, 0.06, size=len(vals))
            ax.scatter(x, vals, s=18, color=c, edgecolor="black", lw=0.3, zorder=3, alpha=0.85)
        ax.set_xticks([1, 2]); ax.set_xticklabels(labels)
        ax.set_ylabel(ylab); ax.set_title(title, fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)

    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    box_points(ax, "elem_per_kb_window", perloc, "E2G elements / kb",
               f"E2G element density in +/-100 kb window\n{CELL}")
    m = [c for c in comp if c["metric"] == "elem_per_kb_window"][0]
    ax.text(0.5, 0.97, f"MWU p={m['p_value']:.3g}", transform=ax.transAxes, ha="center", va="top", fontsize=8)
    fig.tight_layout(); fig.savefig(figd / "e2g_element_density_by_group.pdf"); fig.savefig(figd / "e2g_element_density_by_group.png", dpi=200); plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    box_points(ax, "n_shared_enhancer_any", pair_rows, "shared enhancer elements",
               f"Enhancers shared by BOTH paralogs (E2G)\n{CELL}")
    m = [c for c in comp if c["metric"] == "n_shared_enhancer_any"][0]
    ax.text(0.5, 0.97, f"MWU p={m['p_value']:.3g}", transform=ax.transAxes, ha="center", va="top", fontsize=8)
    fig.tight_layout(); fig.savefig(figd / "e2g_shared_enhancers_by_group.pdf"); fig.savefig(figd / "e2g_shared_enhancers_by_group.png", dpi=200); plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    box_points(ax, "n_distinct_elem_total", pair_rows, "distinct enhancer elements",
               f"Distinct E2G enhancers on pair promoters\n{CELL}")
    m = [c for c in comp if c["metric"] == "n_distinct_elem_total"][0]
    ax.text(0.5, 0.97, f"MWU p={m['p_value']:.3g}", transform=ax.transAxes, ha="center", va="top", fontsize=8)
    fig.tight_layout(); fig.savefig(figd / "e2g_promoter_links_by_group.pdf"); fig.savefig(figd / "e2g_promoter_links_by_group.png", dpi=200); plt.close(fig)

    # ---------- log ----------
    logname = f"07_e2g_promoters.{sub}.log" if sub else "07_e2g_promoters.log"
    with open(ROOT / "log" / logname, "w") as f:
        f.write(f"E2G file: {E2G}\n")
        f.write(f"n links (primary chroms): {n_links}\n")
        f.write(f"n unique elements: {n_unique}\n")
        f.write(f"genome-wide baseline: {base_rate_kb:.5f} elements/kb\n")
        f.write(f"n loci: {len(loci)} (ohnolog={sum(1 for r in loci if r['dup_class4']==G1)}, "
                f"SSD={sum(1 for r in loci if r['dup_class4']==G2)})\n")
        f.write("\ngroup comparison (Mann-Whitney U, clustered_cis_ohnolog vs clustered_SSD_tandem):\n")
        for c in comp:
            f.write(f"  {c['metric']:32s} med_ohno={c['median_ohnolog']:>10} med_SSD={c['median_SSD']:>10} "
                    f"U={c['U']} p={c['p_value']:.4g} r_rb={c['rank_biserial']}\n")

    # console
    print(f"n unique E2G elements: {n_unique}  baseline {base_rate_kb:.5f} elem/kb")
    for c in comp:
        print(f"  {c['metric']:32s} med_o={c['median_ohnolog']:>9} med_s={c['median_SSD']:>9} "
              f"U={c['U']} p={c['p_value']:.4g} r={c['rank_biserial']}")
    print(f"wrote {out}/  and figures in {figd}/")


if __name__ == "__main__":
    main()
