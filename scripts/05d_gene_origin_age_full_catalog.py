#!/usr/bin/env python3
"""
GENE-ORIGIN (phylostratigraphic) age of EVERY paralogous TF relative to the two rounds of
vertebrate whole-genome duplication (R1 / R2), and its distribution across the FOUR
duplication categories (script 06). Full-catalog generalization of script 05b.

Method (identical to 05b): for each TF, find the deepest taxonomic clade in which it still
has an Ensembl-Compara ortholog, using a species ladder that straddles the 2R window, then
bracket gene origin:
  before_R1 : ortholog in a PRE-VERTEBRATE outgroup (yeast/worm/fly/tunicate) -> predates both WGDs
  R1_to_R2  : no invertebrate ortholog, but ortholog in a CYCLOSTOME (lamprey/hagfish)
  after_R2  : ortholog only in GNATHOSTOMES (elephant shark and younger)
  human_only: no ortholog found anywhere on the ladder (unresolved / human-specific)
See 05b header for the 2R framework and the caveat that absence is weaker than presence.

Inputs : results/TF_dup_2x2_classification.tsv  (the 4-group catalog)
Output : results/TF_gene_origin_age_2R.full.tsv             (per-gene call + per-species presence)
         results/intermediate/gene_origin_orthologs.full.tsv (raw ortholog rows, biomart cache)
         results/TF_gene_origin_age_by_category.tsv          (category x origin-window counts/fracs)
         results/figures/TF_gene_origin_age_by_category.{pdf,png}
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/05d_gene_origin_age_full_catalog.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLASS = RES / "TF_dup_2x2_classification.tsv"
AGE = RES / "TF_duplication_age.tsv"
CACHE = RES / "intermediate" / "gene_origin_orthologs.full.tsv"
OUT_GENE = RES / "TF_gene_origin_age_2R.full.tsv"
OUT_CAT = RES / "TF_gene_origin_age_by_category.tsv"
FIG = RES / "figures" / "TF_gene_origin_age_by_category"

# (ensembl species key, label, clade, WGD-window), deepest -> most recent (= 05b LADDER).
LADDER = [
    ("scerevisiae",   "S. cerevisiae",   "Fungi",          "pre_vertebrate"),
    ("celegans",      "C. elegans",      "Protostomia",    "pre_vertebrate"),
    ("dmelanogaster", "D. melanogaster", "Protostomia",    "pre_vertebrate"),
    ("cintestinalis", "C. intestinalis", "Tunicata",       "pre_vertebrate"),
    ("csavignyi",     "C. savignyi",     "Tunicata",       "pre_vertebrate"),
    ("eburgeri",      "E. burgeri",      "Cyclostomata",   "cyclostome"),
    ("pmarinus",      "P. marinus",      "Cyclostomata",   "cyclostome"),
    ("cmilii",        "C. milii",        "Chondrichthyes", "gnathostome"),
    ("lchalumnae",    "L. chalumnae",    "Sarcopterygii",  "gnathostome"),
    ("drerio",        "D. rerio",        "Actinopterygii", "gnathostome"),
    ("xtropicalis",   "X. tropicalis",   "Amphibia",       "gnathostome"),
    ("ggallus",       "G. gallus",       "Aves",           "gnathostome"),
    ("mmusculus",     "M. musculus",     "Mammalia",       "gnathostome"),
]
WINDOWS = ["before_R1", "R1_to_R2", "after_R2", "human_only"]
WLABEL = {"before_R1": "before R1\n(pre-2R ancient)", "R1_to_R2": "R1–R2 window\n(vertebrate stem)",
          "after_R2": "after R2\n(gnathostome)", "human_only": "unresolved\n(no ortholog)"}
WCOLOR = {"before_R1": "#4575b4", "R1_to_R2": "#fdae61", "after_R2": "#d73027",
          "human_only": "#bdbdbd"}
GROUPS = ["dispersed_ohnolog", "clustered_cis_ohnolog", "clustered_SSD_tandem", "dispersed_SSD"]
GLABEL = {"dispersed_ohnolog": "dispersed\nohnolog", "clustered_cis_ohnolog": "clustered\ncis-ohnolog",
          "clustered_SSD_tandem": "clustered\nSSD-tandem", "dispersed_SSD": "dispersed\nSSD"}
CHUNK = 250


def log(*a):
    print("[05d]", *a, flush=True)


def resolve_ids(genes):
    """symbol -> ENSG id; prefer the age table's ensembl_id, fall back to a biomart name pull."""
    ids = {}
    if AGE.exists():
        a = pd.read_csv(AGE, sep="\t")
        if "ensembl_id" in a.columns:
            for g, v in zip(a.gene, a.ensembl_id):
                if g in genes and isinstance(v, str) and v.startswith("ENSG"):
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


def query_species(ds, sp, idlist):
    """presence query for one species, chunked + retried; returns DataFrame[id,ortholog,type,conf]."""
    frames = []
    for i in range(0, len(idlist), CHUNK):
        chunk = idlist[i:i + CHUNK]
        for attempt in range(3):
            try:
                df = ds.query(attributes=["ensembl_gene_id",
                                          f"{sp}_homolog_ensembl_gene",
                                          f"{sp}_homolog_orthology_type",
                                          f"{sp}_homolog_orthology_confidence"],
                              filters={"link_ensembl_gene_id": chunk})
                # biomart collapses to a single column when NO gene in the chunk has an
                # ortholog in this species -> that chunk simply contributes nothing.
                if df.shape[1] < 4:
                    break
                df = df.iloc[:, :4]
                df.columns = ["id", "ortholog", "type", "confidence"]
                frames.append(df[df["id"].isin(chunk)])
                break
            except Exception as e:
                log(f"  {sp} chunk {i}-{i+len(chunk)} attempt {attempt+1} failed: {e}")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["id", "ortholog", "type", "confidence"])


def fetch_orthologs(id2sym):
    if CACHE.exists():
        log(f"using cached {CACHE}")
        return pd.read_csv(CACHE, sep="\t")
    from pybiomart import Server
    ds = Server(host="http://www.ensembl.org")["ENSEMBL_MART_ENSEMBL"]["hsapiens_gene_ensembl"]
    idlist = sorted(id2sym)
    rows = []
    for sp, *_ in LADDER:
        log(f"querying orthologs in {sp} ({len(idlist)} ids, {((len(idlist)-1)//CHUNK)+1} chunks) ...")
        df = query_species(ds, sp, idlist)
        df = df[df["ortholog"].notna() & (df["ortholog"].astype(str) != "")].copy()
        df["species"] = sp
        df["gene"] = df["id"].map(id2sym)
        rows.append(df)
        log(f"  {sp}: {df['gene'].nunique()} TFs with an ortholog")
    res = pd.concat(rows, ignore_index=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(CACHE, sep="\t", index=False)
    log(f"wrote {CACHE} ({len(res)} rows)")
    return res


def classify(present_windows):
    if "pre_vertebrate" in present_windows:
        return "before_R1"
    if "cyclostome" in present_windows:
        return "R1_to_R2"
    if "gnathostome" in present_windows:
        return "after_R2"
    return "human_only"


def main():
    cls = pd.read_csv(CLASS, sep="\t")
    cls = cls[cls.dup_class4.isin(GROUPS)].reset_index(drop=True)
    genes = set(cls.gene)
    log(f"4-group catalog: {len(cls)} TFs")

    id2sym = {}
    for g, i in resolve_ids(genes).items():
        id2sym[i] = g                       # if two symbols map to one id, last wins (rare)
    log(f"resolved {len(set(id2sym.values()))}/{len(genes)} TFs to Ensembl IDs")

    res = fetch_orthologs(id2sym)
    res = res[res["ortholog"].notna() & (res["ortholog"].astype(str) != "")]

    sp_window = {sp: win for sp, _, _, win in LADDER}
    sp_clade = {sp: cl for sp, _, cl, _ in LADDER}
    sp_order = [sp for sp, *_ in LADDER]

    pres = {}      # gene -> set(species)
    for g, sub in res.groupby("gene"):
        pres[g] = set(sub["species"])

    rows = []
    for g in sorted(genes):
        sps = pres.get(g, set())
        wins = {sp_window[s] for s in sps}
        call = classify(wins)
        deepest = next((s for s in sp_order if s in sps), "")
        row = dict(gene=g, origin_window=call,
                   deepest_ortholog_species=deepest, deepest_clade=sp_clade.get(deepest, ""),
                   n_species_with_ortholog=len(sps),
                   invertebrate_ortholog="pre_vertebrate" in wins,
                   cyclostome_ortholog="cyclostome" in wins,
                   gnathostome_ortholog="gnathostome" in wins)
        for s in sp_order:
            row[f"in_{s}"] = int(s in sps)
        rows.append(row)
    g_df = pd.DataFrame(rows).merge(
        cls[["gene", "dup_class4", "provenance", "arrangement", "TF_subfamily", "cluster_id"]],
        on="gene", how="left")
    g_df.to_csv(OUT_GENE, sep="\t", index=False)
    log(f"wrote {OUT_GENE}")

    counts = (g_df.groupby(["dup_class4", "origin_window"]).size().unstack(fill_value=0)
              .reindex(index=GROUPS, columns=WINDOWS, fill_value=0))
    frac = counts.div(counts.sum(1), axis=0)
    out = counts.copy(); out.columns = [f"n_{c}" for c in out.columns]
    for c in WINDOWS:
        out[f"frac_{c}"] = frac[c]
    out.insert(0, "n_total", counts.sum(1))
    out.to_csv(OUT_CAT, sep="\t")
    log(f"wrote {OUT_CAT}")
    log("origin-window x category counts:")
    print(counts.to_string())

    # chi-square: does origin-window depend on category? (drop 'human_only' unresolved bin)
    tab = counts[["before_R1", "R1_to_R2", "after_R2"]]
    chi2, p, dof, _ = stats.chi2_contingency(tab.values)
    log(f"chi-square origin-window x category (3 resolved windows): chi2={chi2:.2f}, dof={dof}, p={p:.2e}")

    make_figure(counts, frac, chi2, p)


def make_figure(counts, frac, chi2, p):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 5.4))
    x = np.arange(len(GROUPS)); xl = [GLABEL[g] for g in GROUPS]

    # (A) 100% stacked origin-window composition
    bottom = np.zeros(len(GROUPS))
    for w in WINDOWS:
        vals = frac[w].values
        axA.bar(x, vals, bottom=bottom, color=WCOLOR[w], width=0.72,
                edgecolor="white", linewidth=0.5, label=WLABEL[w].replace("\n", " "))
        for i, v in enumerate(vals):
            if v > 0.03:
                axA.text(i, bottom[i] + v / 2, f"{v*100:.0f}", ha="center", va="center",
                         fontsize=8, color="white" if w in {"before_R1", "after_R2"} else "black")
        bottom += vals
    axA.set_xticks(x); axA.set_xticklabels(xl, fontsize=9)
    axA.set_ylabel("fraction of category"); axA.set_ylim(0, 1)
    axA.set_title("(A) Gene-origin window per category\n(deepest ortholog, % of category)",
                  fontsize=10, loc="left")
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)

    # (B) raw counts stacked
    bottom = np.zeros(len(GROUPS))
    for w in WINDOWS:
        vals = counts[w].values.astype(float)
        axB.bar(x, vals, bottom=bottom, color=WCOLOR[w], width=0.72,
                edgecolor="white", linewidth=0.5)
        bottom += vals
    for i, g in enumerate(GROUPS):
        axB.text(i, bottom[i] + max(bottom) * 0.01, f"n={int(counts.loc[g].sum())}",
                 ha="center", va="bottom", fontsize=8)
    axB.set_xticks(x); axB.set_xticklabels(xl, fontsize=9)
    axB.set_ylabel("number of TFs")
    axB.set_title(f"(B) Gene-origin window per category (counts)\n"
                  f"chi-square (3 resolved windows) p={p:.1e}", fontsize=10, loc="left")
    for s in ("top", "right"):
        axB.spines[s].set_visible(False)

    handles = [plt.Rectangle((0, 0), 1, 1, color=WCOLOR[w]) for w in WINDOWS]
    fig.legend(handles, [WLABEL[w].replace("\n", " ") for w in WINDOWS], loc="lower center",
               ncol=4, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Gene-origin age (ortholog phylostratigraphy) across the four duplication categories\n"
                 "deepest Ensembl-Compara ortholog on a yeast→mammal ladder straddling the 2R WGD rounds",
                 fontsize=11.5, y=1.01)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(FIG) + ".pdf", bbox_inches="tight")
    fig.savefig(str(FIG) + ".png", dpi=190, bbox_inches="tight")
    log(f"wrote {FIG}.pdf/.png")


if __name__ == "__main__":
    main()
