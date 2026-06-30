#!/usr/bin/env python3
"""
Distribution of TF (duplication) AGE across the FOUR duplication categories (script 06).

Age axis = Ensembl Compara `youngest_paralog_node` carried in
results/TF_dup_2x2_classification.tsv (the most recent within-species paralog-duplication
node of each TF, from script 05). We bin the 21 fine Compara nodes into interpretable
evolutionary EPOCHS (youngest -> oldest), then show, per dup_class4 category:
  (A) 100%-stacked epoch composition  (the age DISTRIBUTION),
  (B) raw counts per epoch,
  (C) box + strip of the ordinal node-rank (older = deeper), with a Kruskal-Wallis test
      across the four categories.

Inputs : results/TF_dup_2x2_classification.tsv
Outputs: results/TF_age_distribution_by_category.tsv      (category x epoch counts + fractions)
         results/figures/TF_age_distribution_by_category.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/05c_TF_age_distribution_by_category.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
OUT_TSV = RES / "TF_age_distribution_by_category.tsv"
FIG = RES / "figures" / "TF_age_distribution_by_category"

# Compara nodes youngest -> oldest (= script 05 NODE_AGE); rank 0 = youngest.
NODE_AGE = [
    "Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini",
    "Simiiformes", "Primates", "Euarchontoglires", "Boreoeutheria", "Eutheria",
    "Theria", "Mammalia", "Amniota", "Tetrapoda", "Sarcopterygii", "Euteleostomi",
    "Vertebrata", "Chordata", "Bilateria", "Opisthokonta", "Eukaryota",
]
RANK = {n: i for i, n in enumerate(NODE_AGE)}

# Fine node -> evolutionary EPOCH (youngest -> oldest). 2R window kept as its own bin.
def epoch_of(node):
    if not isinstance(node, str) or node == "":
        return "no paralog node"
    if node in {"Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini",
                "Simiiformes", "Primates"}:
        return "primate"
    if node in {"Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammalia"}:
        return "mammal"
    if node in {"Amniota", "Tetrapoda", "Sarcopterygii"}:
        return "amniote–tetrapod"
    if node in {"Euteleostomi", "Vertebrata"}:
        return "2R vertebrate"
    if node == "Chordata":
        return "chordate (pre-2R)"
    if node == "Bilateria":
        return "bilaterian"
    if node in {"Opisthokonta", "Eukaryota"}:
        return "ancient (pre-metazoan)"
    return "other"

# epochs youngest -> oldest for plotting, plus the unknown bin last
EPOCHS = ["primate", "mammal", "amniote–tetrapod", "2R vertebrate",
          "chordate (pre-2R)", "bilaterian", "ancient (pre-metazoan)", "no paralog node"]
EPOCH_COLORS = {
    "primate": "#fee08b", "mammal": "#fdae61", "amniote–tetrapod": "#f46d43",
    "2R vertebrate": "#d73027", "chordate (pre-2R)": "#74add1",
    "bilaterian": "#4575b4", "ancient (pre-metazoan)": "#313695",
    "no paralog node": "#bdbdbd"}
GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
GLABEL = {"dispersed_ohnolog": "dispersed\nohnolog", "clustered_cis_ohnolog": "clustered\ncis-ohnolog",
          "clustered_SSD_tandem": "clustered\nSSD-tandem", "dispersed_SSD": "dispersed\nSSD"}


def log(*a):
    print("[05c]", *a, flush=True)


def main():
    m = pd.read_csv(CLASS, sep="\t")
    m = m[m.dup_class4.isin(GROUPS)].copy()
    m["epoch"] = m["youngest_paralog_node"].map(epoch_of)
    m["node_rank"] = m["youngest_paralog_node"].map(RANK)
    log(f"TFs: {len(m)} across {m.dup_class4.nunique()} categories")

    # category x epoch counts + within-category fractions
    counts = (m.groupby(["dup_class4", "epoch"]).size().unstack(fill_value=0)
              .reindex(index=GROUPS, columns=EPOCHS, fill_value=0))
    frac = counts.div(counts.sum(1), axis=0)
    out = counts.copy()
    out.columns = [f"n_{c}" for c in out.columns]
    for c in EPOCHS:
        out[f"frac_{c}"] = frac[c]
    out.insert(0, "n_total", counts.sum(1))
    out.to_csv(OUT_TSV, sep="\t")
    log(f"wrote {OUT_TSV}")

    # Kruskal-Wallis on node_rank across categories (drop genes with no paralog node)
    rk = [m.loc[(m.dup_class4 == g) & m.node_rank.notna(), "node_rank"].values for g in GROUPS]
    H, p = stats.kruskal(*rk)
    log(f"Kruskal-Wallis on youngest-paralog-node rank across 4 categories: H={H:.2f}, p={p:.2e}")

    make_figure(counts, frac, m, rk, H, p)


def make_figure(counts, frac, m, rk, H, p):
    fig = plt.figure(figsize=(15, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1.0], wspace=0.32)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))
    x = np.arange(len(GROUPS))
    xl = [GLABEL[g] for g in GROUPS]

    # ---- (A) 100% stacked epoch composition ----
    bottom = np.zeros(len(GROUPS))
    for e in EPOCHS:
        vals = frac[e].values
        axA.bar(x, vals, bottom=bottom, color=EPOCH_COLORS[e], label=e, width=0.74,
                edgecolor="white", linewidth=0.4)
        for i, v in enumerate(vals):
            if v > 0.045:
                axA.text(i, bottom[i] + v / 2, f"{v*100:.0f}", ha="center", va="center",
                         fontsize=7, color="white" if e in
                         {"2R vertebrate", "bilaterian", "ancient (pre-metazoan)"} else "black")
        bottom += vals
    axA.set_xticks(x); axA.set_xticklabels(xl, fontsize=8.5)
    axA.set_ylabel("fraction of category"); axA.set_ylim(0, 1)
    axA.set_title("(A) Age composition per category\n(youngest-paralog epoch, % of category)",
                  fontsize=10, loc="left")
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)

    # ---- (B) raw counts stacked ----
    bottom = np.zeros(len(GROUPS))
    for e in EPOCHS:
        vals = counts[e].values.astype(float)
        axB.bar(x, vals, bottom=bottom, color=EPOCH_COLORS[e], width=0.74,
                edgecolor="white", linewidth=0.4)
        bottom += vals
    for i, g in enumerate(GROUPS):
        axB.text(i, bottom[i] + max(bottom) * 0.01, f"n={int(counts.loc[g].sum())}",
                 ha="center", va="bottom", fontsize=8)
    axB.set_xticks(x); axB.set_xticklabels(xl, fontsize=8.5)
    axB.set_ylabel("number of TFs")
    axB.set_title("(B) Age composition per category\n(raw counts)", fontsize=10, loc="left")
    for s in ("top", "right"):
        axB.spines[s].set_visible(False)

    # ---- (C) node-rank distribution (box + strip); older = deeper ----
    bp = axC.boxplot(rk, positions=x, widths=0.5, showfliers=False, patch_artist=True,
                     medianprops=dict(color="black", lw=1.4))
    cat_col = {"dispersed_ohnolog": "#2166ac", "clustered_cis_ohnolog": "#5aae61",
               "clustered_SSD_tandem": "#762a83", "dispersed_SSD": "#bababa"}
    for patch, g in zip(bp["boxes"], GROUPS):
        patch.set_facecolor(cat_col[g]); patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate(rk):
        jit = rng.uniform(-0.16, 0.16, size=len(vals))
        axC.scatter(x[i] + jit, vals, s=5, color="black", alpha=0.25, linewidths=0)
    axC.set_xticks(x); axC.set_xticklabels(xl, fontsize=8.5)
    # label a few key ranks on the y-axis
    yticks = [RANK[n] for n in ["Homo sapiens", "Mammalia", "Vertebrata", "Chordata",
                                "Bilateria", "Opisthokonta"]]
    axC.set_yticks(yticks)
    axC.set_yticklabels(["Homo sapiens", "Mammalia", "Vertebrata (2R)", "Chordata",
                         "Bilateria", "Opisthokonta"], fontsize=7.5)
    axC.invert_yaxis()   # youngest at top, oldest at bottom
    axC.set_ylabel("youngest paralog node  (older →)")
    axC.set_title(f"(C) Paralog-age rank by category\nKruskal-Wallis H={H:.1f}, p={p:.1e}",
                  fontsize=10, loc="left")
    for s in ("top", "right"):
        axC.spines[s].set_visible(False)

    handles = [plt.Rectangle((0, 0), 1, 1, color=EPOCH_COLORS[e]) for e in EPOCHS]
    fig.legend(handles, EPOCHS, loc="lower center", ncol=len(EPOCHS), fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.04), title="evolutionary epoch (youngest → oldest)",
               title_fontsize=8.5)
    fig.suptitle("Distribution of TF duplication age across the four duplication categories "
                 "(youngest Ensembl-Compara paralog node; Lambert TFs, KRAB/HOX/readthrough/pseudo removed)",
                 fontsize=11.5, y=1.0)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {FIG}.pdf/.png")


if __name__ == "__main__":
    main()
