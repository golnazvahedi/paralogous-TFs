#!/usr/bin/env python3
"""
Annotate every paralogous TF (clustered + singleton, KRAB/HOX-free) with its
DUPLICATION MODE and DUPLICATION AGE, to test the 2R-WGD / adaptive-immunity hypothesis.

Two reference layers
--------------------
1) OHNOLOGS v2 (Singh & Isambert; ohnologs.curie.fr) human pair lists at three
   confidence levels (crit0=Strict, critA=Intermediate, critC=Relaxed). A TF is an
   **ohnolog** (2R-WGD-derived) if it has an ohnolog partner anywhere in the genome.
   OHNOLOGS dates ~99% of pairs to the "Vertebrata" node = the R1/R2 (2R) window.
   (files reused from ../iteration7_05.20.2026/inputs/ohnologs/)
2) Ensembl Compara within-species paralogy (pybiomart attribute
   "hsapiens_paralog_subtype" = taxonomic node of the duplication). For each TF we take
   the YOUNGEST paralog node (most recent duplication it participated in) as a finer age
   scale, since OHNOLOGS only resolves the 2R bin. Cached to results/intermediate.

Per-TF output columns
---------------------
  gene, class (clustered/singleton), is_ohnolog_{strict,interm,relaxed}, n_ohnolog_partners,
  ohnolog_dup_time, dup_mode (WGD_ohnolog / tandem / SSD_dispersed),
  youngest_paralog_node, youngest_node_rank, has_2R_window_paralog, n_paralogs, age_category

Output : results/TF_duplication_age.tsv
         results/intermediate/ensembl_paralog_nodes.tsv (biomart cache)
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/05_TF_duplication_age.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
CLUST = RES / "clustered_paralogs.noKRABHOX.tsv"
SING = RES / "singleton_paralogs.noKRABHOX.tsv"
OHNO_DIR = ROOT.parent / "iteration7_05.20.2026" / "inputs" / "ohnologs"
CACHE = RES / "intermediate" / "ensembl_paralog_nodes.tsv"
OUT = RES / "TF_duplication_age.tsv"

# Vertebrate lineage nodes, YOUNGEST -> OLDEST (Ensembl Compara taxon names)
NODE_AGE = [
    "Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini",
    "Simiiformes", "Primates", "Euarchontoglires", "Boreoeutheria", "Eutheria",
    "Theria", "Mammalia", "Amniota", "Tetrapoda", "Sarcopterygii", "Euteleostomi",
    "Vertebrata", "Chordata", "Bilateria", "Opisthokonta", "Eukaryota",
]
RANK = {n: i for i, n in enumerate(NODE_AGE)}      # 0 youngest
TWO_R_NODES = {"Euteleostomi", "Vertebrata"}        # canonical 2R window in Ensembl
PRE2R_OR_OLDER = {"Chordata", "Bilateria", "Opisthokonta", "Eukaryota"}


def log(*a):
    print("[24]", *a, flush=True)


def load_universe():
    cl = pd.read_csv(CLUST, sep="\t")
    rows = []
    for _, r in cl.iterrows():
        for g in str(r["members"]).split(","):
            rows.append(dict(gene=g, cls="clustered"))
    sg = pd.read_csv(SING, sep="\t")
    for _, r in sg.iterrows():
        rows.append(dict(gene=r["gene_symbol"], cls="singleton", ensembl_id=r["ensembl_id"]))
    u = pd.DataFrame(rows).drop_duplicates("gene").reset_index(drop=True)
    return u


def load_ohnologs():
    """Return dict gene_symbol(upper) -> (relaxed,interm,strict bools, n_partners, dup_time)."""
    levels = {"0": "strict", "A": "interm", "C": "relaxed"}
    sym_level = {}     # sym -> set of levels
    sym_partners = {}  # sym -> set of partner syms (relaxed)
    sym_time = {}      # sym -> duplication time (relaxed)
    for crit, lvl in levels.items():
        path = OHNO_DIR / f"hsapiens_pairs_crit{crit}.tsv"
        df = pd.read_csv(path, sep="\t")
        for _, r in df.iterrows():
            for s in (str(r["Symbol1"]).upper(), str(r["Symbol2"]).upper()):
                sym_level.setdefault(s, set()).add(lvl)
            if lvl == "relaxed":
                s1, s2 = str(r["Symbol1"]).upper(), str(r["Symbol2"]).upper()
                sym_partners.setdefault(s1, set()).add(s2)
                sym_partners.setdefault(s2, set()).add(s1)
                sym_time.setdefault(s1, r["Duplication time"])
                sym_time.setdefault(s2, r["Duplication time"])
    return sym_level, sym_partners, sym_time


def fetch_paralog_nodes(u):
    if CACHE.exists():
        log(f"using cached {CACHE}")
        return pd.read_csv(CACHE, sep="\t")
    from pybiomart import Server
    srv = Server(host="http://www.ensembl.org")
    ds = srv["ENSEMBL_MART_ENSEMBL"]["hsapiens_gene_ensembl"]

    # 1) symbol -> ensembl_gene_id map (one full pull), with singleton ids as fallback
    log("pulling gene-id map from biomart ...")
    gm = ds.query(attributes=["ensembl_gene_id", "external_gene_name"])
    gm.columns = ["ensembl_gene_id", "gene"]
    sym2id = dict(zip(gm.gene, gm.ensembl_gene_id))
    ids = {}
    for _, r in u.iterrows():
        gid = sym2id.get(r.gene)
        if not gid and "ensembl_id" in u.columns and pd.notna(r.get("ensembl_id")):
            gid = str(r["ensembl_id"]).split(".")[0]
        if gid:
            ids[gid] = r.gene
    id2sym = ids
    log(f"resolved {len(id2sym)}/{len(u)} TFs to Ensembl IDs")

    # 2) homology query filtered by link_ensembl_gene_id, chunked
    log("querying paralog duplication nodes ...")
    out = []
    idlist = sorted(id2sym)
    step = 200
    for i in range(0, len(idlist), step):
        chunk = idlist[i:i + step]
        df = ds.query(attributes=["ensembl_gene_id",
                                   "hsapiens_paralog_associated_gene_name",
                                   "hsapiens_paralog_subtype"],
                      filters={"link_ensembl_gene_id": chunk})
        out.append(df)
        log(f"  {min(i+step,len(idlist))}/{len(idlist)} genes")
    res = pd.concat(out, ignore_index=True)
    res.columns = ["ensembl_gene_id", "paralog", "node"]
    res["gene"] = res["ensembl_gene_id"].map(id2sym)
    res = res.dropna(subset=["gene"])[["gene", "paralog", "node"]]
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(CACHE, sep="\t", index=False)
    log(f"wrote {CACHE} ({len(res)} paralog rows)")
    return res


def main():
    u = load_universe()
    log(f"universe: {len(u)} paralogous TFs ({(u.cls=='clustered').sum()} clustered, "
        f"{(u.cls=='singleton').sum()} singleton)")

    sym_level, sym_partners, sym_time = load_ohnologs()
    up = u.gene.str.upper()
    u["is_ohnolog_strict"] = up.map(lambda s: "strict" in sym_level.get(s, set()))
    u["is_ohnolog_interm"] = up.map(lambda s: "interm" in sym_level.get(s, set()))
    u["is_ohnolog_relaxed"] = up.map(lambda s: "relaxed" in sym_level.get(s, set()))
    u["n_ohnolog_partners"] = up.map(lambda s: len(sym_partners.get(s, set())))
    u["ohnolog_dup_time"] = up.map(lambda s: sym_time.get(s, ""))

    # Ensembl Compara paralog nodes -> youngest node per gene
    nodes = fetch_paralog_nodes(u)
    nodes = nodes.dropna(subset=["node"]).copy()
    nodes["rank"] = nodes["node"].map(RANK)
    nodes = nodes.dropna(subset=["rank"])
    g = nodes.groupby("gene")
    youngest = g["rank"].min()
    npar = g["paralog"].nunique()
    has2r = g["node"].apply(lambda s: bool(set(s) & TWO_R_NODES))
    u["youngest_node_rank"] = u.gene.map(youngest)
    u["youngest_paralog_node"] = u["youngest_node_rank"].map(
        lambda r: NODE_AGE[int(r)] if pd.notna(r) else "")
    u["n_paralogs"] = u.gene.map(npar).fillna(0).astype(int)
    u["has_2R_window_paralog"] = u.gene.map(has2r).fillna(False)

    # duplication mode (WGD ohnolog takes precedence; tandem = our clustered class)
    def mode(row):
        if row.is_ohnolog_relaxed:
            return "WGD_ohnolog"
        return "tandem" if row.cls == "clustered" else "SSD_dispersed"
    u["dup_mode"] = u.apply(mode, axis=1)

    # age category from youngest paralog node
    def age_cat(node):
        if node == "":
            return "no_paralog_node"
        if node in TWO_R_NODES:
            return "2R_vertebrata"
        if node in PRE2R_OR_OLDER:
            return "pre2R_ancient"
        return "post2R_young"
    u["age_category"] = u["youngest_paralog_node"].map(age_cat)

    u.to_csv(OUT, sep="\t", index=False)
    log(f"wrote {OUT}")

    # --- console summary ---
    log("ohnolog status (relaxed) by class:")
    print(pd.crosstab(u.cls, u.is_ohnolog_relaxed).to_string())
    log("duplication mode by class:")
    print(pd.crosstab(u.cls, u.dup_mode).to_string())
    log("age category by class:")
    print(pd.crosstab(u.cls, u.age_category).to_string())
    log(f"ohnolog fraction: clustered={u.loc[u.cls=='clustered','is_ohnolog_relaxed'].mean():.2%}, "
        f"singleton={u.loc[u.cls=='singleton','is_ohnolog_relaxed'].mean():.2%}")


if __name__ == "__main__":
    main()
