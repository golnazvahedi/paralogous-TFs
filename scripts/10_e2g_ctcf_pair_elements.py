#!/usr/bin/env python3
"""
10_e2g_ctcf_pair_elements.py -- CTCF signal carried in the ENCODE-rE2G element-gene
links, read at the PARALOG-PAIR regulatory elements of the two CLUSTERED groups.

Motivation: the main CTCF architecture analysis (scripts 03-06) used a GM12878
(B-lymphoblastoid) IDR narrowPeak set. The ENCODE-rE2G files carry a per-element CTCF
ChIP-seq SIGNAL feature (column index 14, "CTCF") that is matched to each T-cell subset
the rE2G model was built in. This script reads that column at the regulatory elements
that rE2G links to the paralog pair, giving a T-cell-subset CTCF readout to compare
against the GM12878 result.

Data / columns (0-indexed) of the ENCODE-rE2G thresholded links (same layout as script 07):
  0 chr  1 start  2 end  3 name  4 class  5 TargetGene  6 EnsemblID  7 TargetGeneTSS
  8 CellType  9 Score  10 DistanceToTSS  11 H3K27ac  12 Open(ATAC)  13 H3K4me  14 CTCF ...

CTCF is an ELEMENT feature: it is (near-)constant across the multiple element->gene links
of a given element, so it is deduped by element key (chrom,start,end).

Per-locus CTCF metrics (only the two clustered classes):
  - ctcf_mean_pair_enh  : mean CTCF over distinct ENHANCER-class elements rE2G links to a
                          pair member (the "pair-linked enhancer" set)
  - ctcf_max_pair_enh   : max CTCF over that set (strongest insulator-ish element)
  - ctcf_mean_pair_prom : mean CTCF over PROMOTER-class elements at the pair TSSs
  - ctcf_mean_shared_enh: mean CTCF over enhancer elements linked to BOTH paralogs
                          (the shared-regulatory-input set; NaN when a locus shares none)
  - ctcf_mean_window    : mean CTCF over ALL unique E2G elements whose midpoint falls in the
                          +/-100 kb window (locus-level, target-gene-agnostic)

Group comparison: Mann-Whitney U (two-sided), clustered_cis_ohnolog vs clustered_SSD_tandem,
reported with U, p, n per group, medians, and rank-biserial effect size.

Outputs (mirrors 07's --subdir convention so it is drop-in for the panel):
  results/e2g/<subdir>/e2g_ctcf_pair_elements.tsv
  results/e2g/<subdir>/e2g_ctcf_group_comparison.tsv
  results/figures/e2g/<subdir>/e2g_ctcf_pair_enhancers_by_group.{pdf,png}
  log/10_e2g_ctcf_pair_elements.<subdir>.log
"""
from pathlib import Path
import csv, gzip, argparse
from bisect import bisect_left
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
G1, G2 = "clustered_cis_ohnolog", "clustered_SSD_tandem"
COL = {G1: "#0072B2", G2: "#D55E00"}
LAB = {G1: "cis-ohnolog", G2: "SSD-tandem"}
CTCF_IDX = 14


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
       elem_ctcf[chrom] -> (sorted midpoints array, parallel CTCF array) for unique elements
       gene_links[gene] -> list of (chrom, start, end, klass, ctcf)
    CTCF (element feature) deduped per element key; if an element appears with several
    CTCF values (should not happen), keep the max.
    """
    elem_ctcf_map = defaultdict(dict)     # chrom -> {(start,end): ctcf}
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
                ctcf = float(p[CTCF_IDX])
            except (ValueError, IndexError):
                ctcf = float("nan")
            key = (start, end)
            prev = elem_ctcf_map[chrom].get(key)
            if prev is None or (ctcf == ctcf and ctcf > prev):
                elem_ctcf_map[chrom][key] = ctcf
            gene_links[gene].append((chrom, start, end, klass, ctcf))
            n_links += 1
    elem_ctcf = {}
    n_unique = 0
    for c, d in elem_ctcf_map.items():
        items = sorted(((a + b) // 2, v) for (a, b), v in d.items())
        mids = np.array([m for m, _ in items])
        vals = np.array([v for _, v in items], dtype=float)
        elem_ctcf[c] = (mids, vals)
        n_unique += len(mids)
    return elem_ctcf, gene_links, n_unique, n_links


def window_ctcf(elem_ctcf, chrom, a, b):
    """mean CTCF over unique elements with midpoint in [a,b); count too."""
    if chrom not in elem_ctcf or b <= a:
        return float("nan"), 0
    mids, vals = elem_ctcf[chrom]
    lo = bisect_left(mids, a)
    hi = bisect_left(mids, b)
    if hi <= lo:
        return float("nan"), 0
    seg = vals[lo:hi]
    seg = seg[~np.isnan(seg)]
    return (float(np.mean(seg)) if seg.size else float("nan")), int(hi - lo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--e2g", default=str(ROOT / "inputs" / "encode" / "ENCFF393GIF.bed.gz"))
    ap.add_argument("--subdir", default="")
    ap.add_argument("--label", default="Treg ENCODE-rE2G, ENCFF393GIF")
    args = ap.parse_args()
    sub, cell = args.subdir, args.label

    loci = []
    with open(ROOT / "results" / "clustered_loci_windows.pm100kb.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["dup_class4"] in (G1, G2):
                loci.append(row)

    elem_ctcf, gene_links, n_unique, n_links = load_e2g(Path(args.e2g))

    perloc = []
    for r in loci:
        chrom = r["chrom"]
        ws, we = int(r["win_start"]), int(r["win_end"])
        members = [m for m in r["members"].split(",") if m]

        elem_to_members = defaultdict(set)     # (chrom,s,e) -> members linked
        elem_class = {}                        # (chrom,s,e) -> klass (promoter if any link is promoter)
        elem_ctcfv = {}                        # (chrom,s,e) -> ctcf
        for gene in members:
            for (c, s, e, klass, ctcf) in gene_links.get(gene, []):
                key = (c, s, e)
                elem_to_members[key].add(gene)
                elem_ctcfv[key] = ctcf
                if klass == "promoter":
                    elem_class[key] = "promoter"
                else:
                    elem_class.setdefault(key, "enhancer")

        enh_ctcf, prom_ctcf, shared_ctcf = [], [], []
        for key, klass in elem_class.items():
            v = elem_ctcfv[key]
            if klass == "promoter":
                if v == v:
                    prom_ctcf.append(v)
            else:
                if v == v:
                    enh_ctcf.append(v)
                if len(elem_to_members[key]) >= 2 and v == v:
                    shared_ctcf.append(v)

        win_mean, n_win = window_ctcf(elem_ctcf, chrom, ws, we)
        perloc.append(dict(
            cluster_id=r["cluster_id"], dup_class4=r["dup_class4"], chrom=chrom,
            n_members=len(members), members=r["members"],
            n_pair_enh_elem=len(enh_ctcf),
            ctcf_mean_pair_enh=float(np.mean(enh_ctcf)) if enh_ctcf else float("nan"),
            ctcf_max_pair_enh=float(np.max(enh_ctcf)) if enh_ctcf else float("nan"),
            n_pair_prom_elem=len(prom_ctcf),
            ctcf_mean_pair_prom=float(np.mean(prom_ctcf)) if prom_ctcf else float("nan"),
            n_shared_enh=len(shared_ctcf),
            ctcf_mean_shared_enh=float(np.mean(shared_ctcf)) if shared_ctcf else float("nan"),
            n_win_elem=n_win,
            ctcf_mean_window=win_mean,
        ))

    out = ROOT / "results" / "e2g" / sub if sub else ROOT / "results" / "e2g"
    out.mkdir(parents=True, exist_ok=True)
    figd = ROOT / "results" / "figures" / "e2g" / sub if sub else ROOT / "results" / "figures" / "e2g"
    figd.mkdir(parents=True, exist_ok=True)

    with open(out / "e2g_ctcf_pair_elements.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(perloc[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(perloc)

    metrics = [
        ("ctcf_mean_pair_enh", "Mean CTCF signal at pair-linked enhancers"),
        ("ctcf_max_pair_enh", "Max CTCF signal at pair-linked enhancers"),
        ("ctcf_mean_pair_prom", "Mean CTCF signal at pair promoter elements"),
        ("ctcf_mean_shared_enh", "Mean CTCF at enhancers shared by both paralogs"),
        ("ctcf_mean_window", "Mean CTCF over E2G elements in +/-100 kb window"),
    ]
    comp = []
    for key, label in metrics:
        g1 = [float(r[key]) for r in perloc if r["dup_class4"] == G1 and r[key] == r[key]]
        g2 = [float(r[key]) for r in perloc if r["dup_class4"] == G2 and r[key] == r[key]]
        m = mwu(g1, g2)
        comp.append(dict(metric=key, label=label,
                         median_ohnolog=round(m["med1"], 5) if m["med1"] == m["med1"] else float("nan"),
                         median_SSD=round(m["med2"], 5) if m["med2"] == m["med2"] else float("nan"),
                         n_ohnolog=m["n1"], n_SSD=m["n2"], U=m["U"], p_value=m["p"],
                         rank_biserial=round(m["r_rb"], 4) if m["r_rb"] == m["r_rb"] else float("nan")))
    with open(out / "e2g_ctcf_group_comparison.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comp[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(comp)

    # figure: pair-linked-enhancer CTCF by group (box + points)
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    data, labels, colors = [], [], []
    for g in (G1, G2):
        vals = [float(r["ctcf_mean_pair_enh"]) for r in perloc
                if r["dup_class4"] == g and r["ctcf_mean_pair_enh"] == r["ctcf_mean_pair_enh"]]
        data.append(vals); labels.append(f"{LAB[g]}\n(n={len(vals)})"); colors.append(COL[g])
    bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True,
                    medianprops=dict(color="black", lw=1.4))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.35)
    rng = np.random.default_rng(0)
    for i, (vals, c) in enumerate(zip(data, colors), start=1):
        ax.scatter(rng.normal(i, 0.06, size=len(vals)), vals, s=18, color=c,
                   edgecolor="black", lw=0.3, zorder=3, alpha=0.85)
    ax.set_xticks([1, 2]); ax.set_xticklabels(labels)
    ax.set_ylabel("mean CTCF signal (rE2G element feature)")
    ax.set_title(f"CTCF at pair-linked enhancers\n{cell}", fontsize=10)
    m = [c for c in comp if c["metric"] == "ctcf_mean_pair_enh"][0]
    pv = m["p_value"]
    ax.text(0.5, 0.98, f"MWU p={pv:.3g}" if pv == pv else "MWU p=NA",
            transform=ax.transAxes, ha="center", va="top", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(figd / "e2g_ctcf_pair_enhancers_by_group.pdf")
    fig.savefig(figd / "e2g_ctcf_pair_enhancers_by_group.png", dpi=200)
    plt.close(fig)

    logname = f"10_e2g_ctcf_pair_elements.{sub}.log" if sub else "10_e2g_ctcf_pair_elements.log"
    with open(ROOT / "log" / logname, "w") as f:
        f.write(f"E2G file: {args.e2g}\nlabel: {cell}\n")
        f.write(f"n links (primary chroms): {n_links}\nn unique elements: {n_unique}\n")
        f.write(f"n loci: {len(loci)} (ohnolog={sum(1 for r in loci if r['dup_class4']==G1)}, "
                f"SSD={sum(1 for r in loci if r['dup_class4']==G2)})\n")
        f.write("\nMann-Whitney U (clustered_cis_ohnolog vs clustered_SSD_tandem):\n")
        for c in comp:
            f.write(f"  {c['metric']:22s} med_ohno={c['median_ohnolog']:>9} med_SSD={c['median_SSD']:>9} "
                    f"n=({c['n_ohnolog']},{c['n_SSD']}) U={c['U']} p={c['p_value']:.4g} r_rb={c['rank_biserial']}\n")

    print(f"[{cell}] n unique elements={n_unique}")
    for c in comp:
        print(f"  {c['metric']:22s} med_o={c['median_ohnolog']:>9} med_s={c['median_SSD']:>9} "
              f"n=({c['n_ohnolog']},{c['n_SSD']}) p={c['p_value']:.4g} r={c['rank_biserial']}")
    print(f"wrote {out}/e2g_ctcf_pair_elements.tsv, e2g_ctcf_group_comparison.tsv")


if __name__ == "__main__":
    main()
