#!/usr/bin/env python3
"""
Co-expression vs. single-member dominance of CLUSTERED paralog clusters in human PBMC.

Question (from the dosage/buffering thread): within a genomic cluster of same-family paralogs,
how often are ALL members co-expressed in a given cell type, versus how often is only ONE
member expressed / dominant? Contrast clustered cis-OHNOLOGS (2R, co-regulated; e.g. ETS1/FLI1,
STAT cluster) against clustered SSD-TANDEM (e.g. chr2 SAND SP100/110/140/140L).

Per gene x PBMC major-lineage cell type we compute, from human_pbmc.h5ad:
  * pct  = fraction of cells in the cell type with detectable counts  (the "expressed" call)
  * cp10k= mean linear CP10K                                          (the "how much" / dominance)
A member is ON in a cell type if pct >= PCT_THR (default 0.10). For each (cluster, cell type)
that is ACTIVE (>=1 member ON) we record:
  n_on / n_members, all_on (all members ON), single_on (exactly 1 ON),
  top_share = max member cp10k / sum member cp10k  (dominance; 1/n .. 1),
  one_dominant = top_share >= DOM_THR (default 0.80).
Aggregated per category (cis-ohnolog vs SSD-tandem) over active (cluster, cell-type) instances,
with Fisher tests on all_on / single_on and Mann-Whitney on top_share; a size-2-only ("pair")
view is reported alongside the all-sizes view.

Universe = clustered paralog clusters with >=2 PBMC-detectable members (results/
TF_dup_2x2_classification.tsv). Cell types = major_lineage (as in script 07c).

Outputs:
  results/TF_cluster_coexpression.per_instance.tsv   one row per (cluster, cell type) active
  results/TF_cluster_coexpression.per_cluster.tsv     per-cluster summary across cell types
  results/TF_cluster_coexpression.summary.tsv         per-category stats + tests
  results/figures/TF_cluster_coexpression.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/16_clustered_ohnolog_coexpression.py [PCT_THR] [DOM_THR]
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
H5AD = RES / "intermediate" / "h5ad" / "human_pbmc.h5ad"

PCT_THR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.10
DOM_THR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.80
OUT_INST = RES / "TF_cluster_coexpression.per_instance.tsv"
OUT_CLUST = RES / "TF_cluster_coexpression.per_cluster.tsv"
OUT_SUM = RES / "TF_cluster_coexpression.summary.tsv"
OUT_FIG = RES / "figures" / "TF_cluster_coexpression"

CT_COL = "major_lineage"
CT_ORDER = ["progenitor", "monocyte", "cDC1", "cDC2", "DC_other", "pDC",
            "B", "plasma", "CD4 T", "CD8 T", "other T", "NK", "ILC",
            "platelet", "erythrocyte"]
CATS = {"clustered_cis_ohnolog": "cis-ohnolog", "clustered_SSD_tandem": "SSD-tandem"}
CCOL = {"clustered_cis_ohnolog": "#5aae61", "clustered_SSD_tandem": "#762a83"}


def log(*a):
    print("[16]", *a, flush=True)


def pseudobulk(genes):
    """Return (cts, pct, cp10k): pct & cp10k are (n_genes x n_celltypes)."""
    log(f"reading {H5AD} (backed) ...")
    A = ad.read_h5ad(H5AD, backed="r")
    vn = pd.Index(A.var_names)
    ct = A.obs[CT_COL].astype(str).values
    cts = [c for c in CT_ORDER if c in set(ct)] + [c for c in pd.unique(ct) if c not in CT_ORDER]
    idx = vn.get_indexer(genes)
    X = A[:, idx].to_memory().X
    X = X.tocsc() if sp.issparse(X) else sp.csr_matrix(X).tocsc()

    log("per-cell library size (chunked) ...")
    n = A.shape[0]
    lib = np.zeros(n)
    step = 20000
    for i in range(0, n, step):
        chunk = A[i:i + step].to_memory().X
        chunk = chunk if sp.issparse(chunk) else sp.csr_matrix(chunk)
        lib[i:i + step] = np.asarray(chunk.sum(1)).ravel()
    lib[lib == 0] = 1.0

    Xcsr = X.tocsr()
    cp = Xcsr.multiply(1e4 / lib[:, None]).tocsr()
    pct = np.zeros((len(genes), len(cts)))
    cp10k = np.zeros((len(genes), len(cts)))
    for j, c in enumerate(cts):
        m = ct == c
        if int(m.sum()) == 0:
            continue
        sub = Xcsr[m]
        pct[:, j] = np.asarray((sub > 0).mean(0)).ravel()
        cp10k[:, j] = np.asarray(cp[m].mean(0)).ravel()
    return cts, pct, cp10k


def main():
    c = pd.read_csv(CLASS, sep="\t")
    det = set(ad.read_h5ad(H5AD, backed="r").var_names)
    cl = c[(c.arrangement == "clustered") & (c.dup_class4.isin(CATS)) & (c.gene.isin(det))].copy()
    # clusters with >=2 detectable members
    grp = cl.groupby("cluster_id")
    clusters = {cid: (sub.dup_class4.iloc[0], list(sub.gene))
                for cid, sub in grp if len(sub) >= 2}
    log(f"clusters with >=2 detectable members: {len(clusters)} "
        f"({sum(v[0]=='clustered_cis_ohnolog' for v in clusters.values())} cis-ohnolog, "
        f"{sum(v[0]=='clustered_SSD_tandem' for v in clusters.values())} SSD-tandem)")

    genes = sorted({g for _, mem in clusters.values() for g in mem})
    gi = {g: k for k, g in enumerate(genes)}
    cts, pct, cp10k = pseudobulk(np.array(genes))

    rows = []
    for cid, (cat, mem) in clusters.items():
        mi = [gi[g] for g in mem]
        for j, ctn in enumerate(cts):
            p = pct[mi, j]
            e = cp10k[mi, j]
            on = p >= PCT_THR
            n_on = int(on.sum())
            if n_on == 0:
                continue   # cluster inactive in this cell type
            tot = e.sum()
            top_share = float(e.max() / tot) if tot > 0 else np.nan
            rows.append(dict(cluster_id=cid, dup_class4=cat, cell_type=ctn,
                             n_members=len(mem), n_on=n_on, frac_on=n_on / len(mem),
                             all_on=bool(n_on == len(mem)), single_on=bool(n_on == 1),
                             top_share=top_share, one_dominant=bool(top_share >= DOM_THR),
                             top_gene=mem[int(e.argmax())],
                             members="|".join(mem)))
    inst = pd.DataFrame(rows)
    inst.to_csv(OUT_INST, sep="\t", index=False)
    log(f"active (cluster x cell-type) instances: {len(inst)}")

    # per-cluster summary across active cell types
    cr = []
    for cid, (cat, mem) in clusters.items():
        s = inst[inst.cluster_id == cid]
        if len(s) == 0:
            continue
        cr.append(dict(cluster_id=cid, dup_class4=cat, n_members=len(mem), members="|".join(mem),
                       n_active_celltypes=len(s), frac_all_on=s.all_on.mean(),
                       frac_single_on=s.single_on.mean(), mean_top_share=s.top_share.mean()))
    pd.DataFrame(cr).sort_values(["dup_class4", "mean_top_share"]).to_csv(OUT_CLUST, sep="\t", index=False)

    # ---- per-category summary + tests ----
    def summarize(df, tag):
        out = []
        for cat in CATS:
            s = df[df.dup_class4 == cat]
            out.append(dict(view=tag, dup_class4=cat, n_instances=len(s),
                            pct_all_on=s.all_on.mean(), pct_single_on=s.single_on.mean(),
                            pct_one_dominant=s.one_dominant.mean(),
                            mean_frac_on=s.frac_on.mean(),
                            mean_top_share=s.top_share.mean(), median_top_share=s.top_share.median()))
        return pd.DataFrame(out)

    summ = pd.concat([summarize(inst, "all_sizes"),
                      summarize(inst[inst.n_members == 2], "pairs_only")], ignore_index=True)

    tests = []
    for tag, df in [("all_sizes", inst), ("pairs_only", inst[inst.n_members == 2])]:
        a = df[df.dup_class4 == "clustered_cis_ohnolog"]
        b = df[df.dup_class4 == "clustered_SSD_tandem"]
        # all_on Fisher
        ct_tab = [[a.all_on.sum(), len(a) - a.all_on.sum()], [b.all_on.sum(), len(b) - b.all_on.sum()]]
        _, p_all = stats.fisher_exact(ct_tab)
        ct_tab2 = [[a.single_on.sum(), len(a) - a.single_on.sum()],
                   [b.single_on.sum(), len(b) - b.single_on.sum()]]
        _, p_single = stats.fisher_exact(ct_tab2)
        u, p_ts = stats.mannwhitneyu(a.top_share.dropna(), b.top_share.dropna(), alternative="two-sided")
        tests.append(dict(view=tag, fisher_all_on_p=p_all, fisher_single_on_p=p_single,
                          mannwhitney_top_share_p=p_ts,
                          cis_pct_all_on=a.all_on.mean(), ssd_pct_all_on=b.all_on.mean(),
                          cis_median_top_share=a.top_share.median(),
                          ssd_median_top_share=b.top_share.median()))
    summ.to_csv(OUT_SUM, sep="\t", index=False)
    log("per-category summary:")
    print(summ.to_string(index=False))
    log("tests (cis-ohnolog vs SSD-tandem):")
    print(pd.DataFrame(tests).to_string(index=False))

    # notable clusters
    log("notable clusters (per-cluster co-expression):")
    pc = pd.read_csv(OUT_CLUST, sep="\t")
    for cid in pc.cluster_id:
        s = pc[pc.cluster_id == cid].iloc[0]
        if any(g in s.members for g in ("STAT1", "ETS1", "SP100", "STAT3", "IRF")):
            print(f"  {s.members} [{CATS[s.dup_class4]}]: active in {s.n_active_celltypes} cts, "
                  f"all_on {s.frac_all_on:.2f}, single_on {s.frac_single_on:.2f}, "
                  f"mean top_share {s.mean_top_share:.2f}")

    make_figure(inst, summ, pd.DataFrame(tests))
    log(f"wrote {OUT_INST} / {OUT_CLUST} / {OUT_SUM}")


def make_figure(inst, summ, tests):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.5, 5.0))
    s = summ[summ.view == "all_sizes"].set_index("dup_class4")
    t = tests[tests.view == "all_sizes"].iloc[0]
    cats = list(CATS)
    xlab = [CATS[c] for c in cats]

    # (A) % all_on vs % single_on
    x = np.arange(len(cats)); w = 0.38
    axA.bar(x - w / 2, [s.loc[c, "pct_all_on"] for c in cats], w, color=[CCOL[c] for c in cats],
            label="all members ON")
    axA.bar(x + w / 2, [s.loc[c, "pct_single_on"] for c in cats], w,
            color=[CCOL[c] for c in cats], alpha=0.45, hatch="//", label="only one ON")
    for i, c in enumerate(cats):
        axA.text(i - w / 2, s.loc[c, "pct_all_on"] + 0.01, f"{s.loc[c,'pct_all_on']:.2f}", ha="center", fontsize=8)
        axA.text(i + w / 2, s.loc[c, "pct_single_on"] + 0.01, f"{s.loc[c,'pct_single_on']:.2f}", ha="center", fontsize=8)
    axA.set_xticks(x); axA.set_xticklabels(xlab)
    axA.set_ylabel("fraction of active (cluster x cell type) instances")
    axA.set_title(f"(A) All members co-expressed vs. one only\nFisher all-ON p={t.fisher_all_on_p:.1e}",
                  fontsize=10)
    axA.legend(fontsize=8, frameon=False); axA.set_ylim(0, 1)
    for sp_ in ("top", "right"):
        axA.spines[sp_].set_visible(False)

    # (B) dominance (top_share) distribution
    data = [inst.loc[inst.dup_class4 == c, "top_share"].dropna().values for c in cats]
    bp = axB.boxplot(data, positions=x, widths=0.5, showfliers=False, patch_artist=True,
                     medianprops=dict(color="black", lw=1.4))
    for patch, c in zip(bp["boxes"], cats):
        patch.set_facecolor(CCOL[c]); patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data):
        axB.scatter(x[i] + rng.uniform(-0.15, 0.15, len(d)), d, s=5, color="black", alpha=0.2)
    axB.axhline(0.5, color="grey", ls=":", lw=1)
    axB.set_xticks(x); axB.set_xticklabels(xlab)
    axB.set_ylabel("dominance: top member's share of cluster expression")
    axB.set_title(f"(B) Single-member dominance\nMann-Whitney p={t.mannwhitney_top_share_p:.1e}", fontsize=10)
    for sp_ in ("top", "right"):
        axB.spines[sp_].set_visible(False)

    # (C) composition of active instances: all_on / partial / single
    axC2 = axC
    comp = []
    for c in cats:
        sub = inst[inst.dup_class4 == c]
        single = (sub.single_on).mean()
        allon = (sub.all_on).mean()
        partial = 1 - single - allon
        comp.append((allon, partial, single))
    comp = np.array(comp)
    labels = ["all ON", "partial", "only one ON"]
    colors = ["#1b7837", "#a6dba0", "#c2a5cf"]
    bottom = np.zeros(len(cats))
    for k in range(3):
        axC2.bar(x, comp[:, k], bottom=bottom, color=colors[k], label=labels[k], width=0.6,
                 edgecolor="white")
        for i in range(len(cats)):
            if comp[i, k] > 0.05:
                axC2.text(i, bottom[i] + comp[i, k] / 2, f"{comp[i,k]:.2f}", ha="center", va="center", fontsize=8)
        bottom += comp[:, k]
    axC2.set_xticks(x); axC2.set_xticklabels(xlab)
    axC2.set_ylabel("fraction of active instances"); axC2.set_ylim(0, 1)
    axC2.set_title("(C) Expression-pattern composition", fontsize=10)
    axC2.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
    for sp_ in ("top", "right"):
        axC2.spines[sp_].set_visible(False)

    fig.suptitle("Co-expression vs. single-member dominance of clustered paralog clusters in human PBMC "
                 f"(member ON if detected in >={PCT_THR:.0%} of a cell type; dominance = top member's "
                 "expression share)", fontsize=11, y=1.02)
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(OUT_FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(OUT_FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {OUT_FIG}.pdf/.png")


if __name__ == "__main__":
    main()
