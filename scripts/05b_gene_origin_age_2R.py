#!/usr/bin/env python3
"""
GENE-ORIGIN (phylostratigraphic) age of selected paralogous TFs relative to the two
rounds of vertebrate whole-genome duplication (R1 / R2).

Motivation
----------
script 05's `age_category` is built from the YOUNGEST *paralog* duplication node -- it dates
the most recent in-genome duplication, NOT when the gene LINEAGE first appeared. For
singleton (dispersed_SSD) TFs that retained no 2R ohnolog (e.g. BCL6, IRF1, KLF2) we instead
want the gene's ORIGIN: the deepest taxonomic clade in which it still has an ortholog. We get
that from Ensembl Compara orthology across a species ladder that straddles the 2R window.

2R framework / decision rule (gene ORIGIN by deepest ortholog)
--------------------------------------------------------------
Modern consensus (Simakov 2020; Nakatani 2021): ONE WGD (call it R1) in the last common
ancestor of ALL vertebrates (shared by cyclostomes -- lamprey/hagfish -- and gnathostomes),
then a SECOND WGD (R2) in the gnathostome stem AFTER the cyclostome split. So presence in a
cyclostome but not invertebrates brackets a gene to the R1..R2 window; presence only in
gnathostomes brackets it to after R2; presence in any pre-vertebrate outgroup means it
predates R1 entirely.

  before_R1   : ortholog in a PRE-VERTEBRATE outgroup (yeast / worm / fly / tunicate)
                -> gene family is older than both WGDs (pre-2R ancient singleton).
  R1_to_R2    : NO invertebrate ortholog, but ortholog in a CYCLOSTOME (lamprey/hagfish)
                -> arose in the vertebrate ancestor, present before the gnathostome 2R.
  after_R2    : ortholog only in GNATHOSTOMES (elephant shark and younger), absent in
                cyclostomes + invertebrates -> gnathostome-specific, arose after R2.

CAVEAT: ortholog ABSENCE is weaker than presence (cyclostome genomes are GC-biased,
fragmentary, lossy). `before_R1` rests on positive invertebrate presence and is robust;
`R1_to_R2` vs `after_R2` hinges on a cyclostome call that can be a false negative, so we also
report the raw per-species presence matrix and an orthology-confidence flag.

Species ladder (deepest clade -> most recent), Ensembl vertebrates Compara:
  scerevisiae(Fungi) celegans,dmelanogaster(Protostomia) cintestinalis,csavignyi(Tunicata)
  | R1 |  eburgeri,pmarinus(Cyclostomata)  | R2 |  cmilii(Chondrichthyes) lchalumnae,drerio
  (bony fish) xtropicalis(Amphibia) ggallus(Aves) mmusculus(Mammalia)

Inputs : results/TF_dup_2x2_classification.tsv (for class context); genes via argv or default.
Output : results/TF_gene_origin_age_2R.tsv          (per-gene call + per-species presence)
         results/intermediate/gene_origin_orthologs.tsv   (raw ortholog rows, biomart cache)
         results/figures/TF_gene_origin_age_2R.{pdf,png}  (presence matrix + call)
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/05b_gene_origin_age_2R.py [GENE ...]
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
AGE = RES / "TF_duplication_age.tsv"
CACHE = RES / "intermediate" / "gene_origin_orthologs.tsv"
OUT = RES / "TF_gene_origin_age_2R.tsv"
FIG = RES / "figures" / "TF_gene_origin_age_2R"

DEFAULT_GENES = ["BCL6", "IRF1", "KLF2"]

# (ensembl species key, display label, clade, WGD-window class), deepest -> most recent.
LADDER = [
    ("scerevisiae",   "S. cerevisiae (yeast)",    "Fungi",          "pre_vertebrate"),
    ("celegans",      "C. elegans (nematode)",    "Protostomia",    "pre_vertebrate"),
    ("dmelanogaster", "D. melanogaster (fly)",    "Protostomia",    "pre_vertebrate"),
    ("cintestinalis", "C. intestinalis (tunicate)", "Tunicata",     "pre_vertebrate"),
    ("csavignyi",     "C. savignyi (tunicate)",   "Tunicata",       "pre_vertebrate"),
    ("eburgeri",      "E. burgeri (hagfish)",     "Cyclostomata",   "cyclostome"),
    ("pmarinus",      "P. marinus (lamprey)",     "Cyclostomata",   "cyclostome"),
    ("cmilii",        "C. milii (elephant shark)", "Chondrichthyes", "gnathostome"),
    ("lchalumnae",    "L. chalumnae (coelacanth)", "Sarcopterygii", "gnathostome"),
    ("drerio",        "D. rerio (zebrafish)",     "Actinopterygii", "gnathostome"),
    ("xtropicalis",   "X. tropicalis (frog)",     "Amphibia",       "gnathostome"),
    ("ggallus",       "G. gallus (chicken)",      "Aves",           "gnathostome"),
    ("mmusculus",     "M. musculus (mouse)",      "Mammalia",       "gnathostome"),
]
WINDOW_LABEL = {"before_R1": "before R1 (pre-2R ancient)",
                "R1_to_R2": "between R1 and R2 (vertebrate stem)",
                "after_R2": "after R2 (gnathostome-specific)",
                "human_only": "no ortholog found (human-only / unresolved)"}


def log(*a):
    print("[05b]", *a, flush=True)


def resolve_ids(genes):
    """symbol -> ENSG id, preferring the age table, else a biomart name lookup."""
    ids = {}
    if AGE.exists():
        a = pd.read_csv(AGE, sep="\t")
        m = dict(zip(a.gene, a.get("ensembl_id", pd.Series(dtype=str))))
        for g in genes:
            v = m.get(g)
            if isinstance(v, str) and v.startswith("ENSG"):
                ids[g] = v.split(".")[0]
    missing = [g for g in genes if g not in ids]
    if missing:
        from pybiomart import Server
        ds = Server(host="http://www.ensembl.org")["ENSEMBL_MART_ENSEMBL"]["hsapiens_gene_ensembl"]
        gm = ds.query(attributes=["ensembl_gene_id", "external_gene_name"])
        gm.columns = ["id", "sym"]
        s2i = dict(zip(gm.sym, gm.id))
        for g in missing:
            if g in s2i:
                ids[g] = s2i[g]
    return ids


def fetch_orthologs(id2sym):
    if CACHE.exists():
        log(f"using cached {CACHE}")
        return pd.read_csv(CACHE, sep="\t")
    from pybiomart import Server
    ds = Server(host="http://www.ensembl.org")["ENSEMBL_MART_ENSEMBL"]["hsapiens_gene_ensembl"]
    idlist = sorted(id2sym)
    rows = []
    for sp, _, clade, window in LADDER:
        log(f"querying orthologs in {sp} ...")
        df = ds.query(attributes=["ensembl_gene_id",
                                  f"{sp}_homolog_ensembl_gene",
                                  f"{sp}_homolog_orthology_type",
                                  f"{sp}_homolog_orthology_confidence"],
                      filters={"link_ensembl_gene_id": idlist})
        df.columns = ["id", "ortholog", "type", "confidence"]
        df = df[df["id"].isin(idlist)].copy()
        df["species"] = sp
        df["gene"] = df["id"].map(id2sym)
        rows.append(df)
    res = pd.concat(rows, ignore_index=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(CACHE, sep="\t", index=False)
    log(f"wrote {CACHE} ({len(res)} rows)")
    return res


def classify(present_windows):
    """present_windows: set of WGD-window classes with >=1 ortholog."""
    if "pre_vertebrate" in present_windows:
        return "before_R1"
    if "cyclostome" in present_windows:
        return "R1_to_R2"
    if "gnathostome" in present_windows:
        return "after_R2"
    return "human_only"


def main():
    genes = [g.upper() for g in sys.argv[1:]] or DEFAULT_GENES
    id2sym = {v: k for k, v in resolve_ids(genes).items()}
    log(f"resolved {len(id2sym)}/{len(genes)} genes: " +
        ", ".join(f"{s}={i}" for i, s in id2sym.items()))

    res = fetch_orthologs(id2sym)
    res = res[res["ortholog"].notna() & (res["ortholog"].astype(str) != "")]

    sp_window = {sp: win for sp, _, _, win in LADDER}
    sp_order = [sp for sp, *_ in LADDER]

    # presence matrix gene x species (1 if >=1 ortholog of any type)
    pres = pd.DataFrame(0, index=genes, columns=sp_order)
    hiconf = pd.DataFrame(0, index=genes, columns=sp_order)
    for (g, sp), sub in res.groupby(["gene", "species"]):
        if g in pres.index and sp in pres.columns:
            pres.loc[g, sp] = 1
            hiconf.loc[g, sp] = int((sub["confidence"].astype(str) == "1").any())

    out_rows = []
    for g in genes:
        wins = {sp_window[sp] for sp in sp_order if pres.loc[g, sp] == 1}
        call = classify(wins)
        deepest = next((sp for sp in sp_order if pres.loc[g, sp] == 1), "")
        row = dict(gene=g, origin_window=call, origin_window_label=WINDOW_LABEL[call],
                   deepest_ortholog_species=deepest,
                   deepest_clade=dict((s, c) for s, _, c, _ in LADDER).get(deepest, ""),
                   n_species_with_ortholog=int(pres.loc[g].sum()),
                   invertebrate_ortholog=bool({"pre_vertebrate"} & wins),
                   cyclostome_ortholog=bool({"cyclostome"} & wins),
                   gnathostome_ortholog=bool({"gnathostome"} & wins))
        for sp in sp_order:
            row[f"in_{sp}"] = int(pres.loc[g, sp])
        out_rows.append(row)
    res_df = pd.DataFrame(out_rows)

    # join the existing duplication context
    if CLASS.exists():
        cls = pd.read_csv(CLASS, sep="\t")[["gene", "dup_class4", "TF_subfamily", "cluster_id"]]
        res_df = res_df.merge(cls, on="gene", how="left")
    if AGE.exists():
        a = pd.read_csv(AGE, sep="\t")[["gene", "youngest_paralog_node", "age_category", "n_paralogs"]]
        res_df = res_df.merge(a, on="gene", how="left")
    res_df.to_csv(OUT, sep="\t", index=False)
    log(f"wrote {OUT}")

    cols = ["gene", "dup_class4", "origin_window_label", "deepest_ortholog_species",
            "invertebrate_ortholog", "cyclostome_ortholog", "youngest_paralog_node"]
    print(res_df[[c for c in cols if c in res_df.columns]].to_string(index=False))

    make_figure(genes, pres, hiconf, res_df, sp_order)


def make_figure(genes, pres, hiconf, res_df, sp_order):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    labels = [lab for _, lab, _, _ in LADDER]
    windows = [w for *_, w in LADDER]
    wcol = {"pre_vertebrate": "#bdbdbd", "cyclostome": "#fdae61", "gnathostome": "#74add1"}

    M = pres.loc[genes, sp_order].values.astype(float)
    fig, ax = plt.subplots(figsize=(10.5, 0.7 * len(genes) + 2.6))
    ax.imshow(np.ones_like(M), cmap="Greys", vmin=0, vmax=1, aspect="auto")
    for i in range(len(genes)):
        for j, sp in enumerate(sp_order):
            if M[i, j]:
                c = wcol[windows[j]]
                ax.scatter(j, i, s=230, marker="s", color=c,
                           edgecolor="black", linewidth=0.6, zorder=3)
                if hiconf.loc[genes[i], sp]:
                    ax.scatter(j, i, s=22, marker="o", color="black", zorder=4)

    # WGD boundary lines: R1 between tunicates & cyclostomes; R2 between cyclostomes & gnathostomes
    r1 = max(j for j, w in enumerate(windows) if w == "pre_vertebrate") + 0.5
    r2 = max(j for j, w in enumerate(windows) if w == "cyclostome") + 0.5
    for x, lab in [(r1, "R1"), (r2, "R2")]:
        ax.axvline(x, color="crimson", lw=2.2, ls="--", zorder=5)
        ax.text(x, -0.85, lab, color="crimson", ha="center", va="bottom",
                fontsize=12, fontweight="bold")

    ax.set_xticks(range(len(sp_order)))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ylab = [f"{g}\n[{res_df.set_index('gene').loc[g,'origin_window']}]" for g in genes]
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(ylab, fontsize=9)
    ax.set_xlim(-0.5, len(sp_order) - 0.5); ax.set_ylim(len(genes) - 0.5, -1.2)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title("Gene-origin (deepest ortholog) of singleton TFs relative to the 2R WGD rounds\n"
                 "filled square = ortholog present (Ensembl Compara); dot = high-confidence ortholog",
                 fontsize=11, loc="left")
    leg = [Patch(facecolor=wcol["pre_vertebrate"], edgecolor="k", label="pre-vertebrate (before R1)"),
           Patch(facecolor=wcol["cyclostome"], edgecolor="k", label="cyclostome (R1–R2 window)"),
           Patch(facecolor=wcol["gnathostome"], edgecolor="k", label="gnathostome (after R2)")]
    ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -1.15),
              fontsize=8.5, frameon=False, ncol=3)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {FIG}.pdf/.png")


if __name__ == "__main__":
    main()
