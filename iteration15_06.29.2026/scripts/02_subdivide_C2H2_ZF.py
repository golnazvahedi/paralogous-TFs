#!/usr/bin/env python3
"""
Part I, step 2: Subdivide the large 'C2H2 ZF' DBD family more finely.

Lambert et al. lump all 747 paralogous C2H2 zinc-finger TFs under a single DBD
label 'C2H2 ZF'. That family is dominated by tandem-C2H2 effector subfamilies
defined by an accessory repression/oligomerization domain. We split them using
HGNC gene-group membership for the three canonical C2H2-ZF effector domains:

    KRAB  -> 'KRAB domain containing'   (Kruppel-associated box; KRAB-ZNF)
    SCAN  -> 'SCAN domain containing'
    BTB   -> 'BTB domain containing'     (a.k.a. POZ; e.g. ZBTB family)

A single fine subfamily label is assigned by priority KRAB > SCAN > BTB; a C2H2
ZF TF carrying none of these accessory domains stays 'C2H2 ZF (other)'. Genes
can carry more than one accessory domain (e.g. KRAB+SCAN), so per-gene boolean
flags has_KRAB / has_SCAN / has_BTB are kept alongside the single label.

HGNC is joined to Lambert by Ensembl gene ID first, then by current symbol, then
by previous/alias symbol.

Inputs
------
  results/paralogous_TFs.tsv          (from step 1)
  inputs/hgnc_complete_set.txt        (downloaded HGNC complete set)

Outputs
-------
  results/paralogous_TFs.tsv                       (updated, extra columns)
  results/intermediate/c2h2_zf_subfamily_assignment.tsv   (per-gene detail)
  results/C2H2_ZF_subfamily_summary.tsv            (counts per fine subfamily)
  log/02_subdivide_C2H2_ZF.log
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
INT = RES / "intermediate"
LOG = ROOT / "log" / "02_subdivide_C2H2_ZF.log"
LOG.parent.mkdir(parents=True, exist_ok=True)

_logf = open(LOG, "w")
def log(*a):
    msg = " ".join(str(x) for x in a)
    print(msg); _logf.write(msg + "\n"); _logf.flush()

# HGNC gene-group name -> short accessory-domain tag
GROUP2TAG = {
    "KRAB domain containing": "KRAB",
    "SCAN domain containing": "SCAN",
    "BTB domain containing": "BTB",
    "Zinc fingers C2H2-type": "C2H2",   # for sanity-checking the join
}
# single-label priority (most specific / dominant first)
PRIORITY = ["KRAB", "SCAN", "BTB"]

log("# Step 2: subdivide C2H2 ZF using HGNC accessory-domain gene groups")

# --- load step-1 paralog table ------------------------------------------------
par = pd.read_csv(RES / "paralogous_TFs.tsv", sep="\t", dtype=str)
par["family_size"] = par["family_size"].astype(int)
c2h2 = par[par["DBD_family"] == "C2H2 ZF"].copy()
log(f"paralogous TFs total: {len(par)};  C2H2 ZF: {len(c2h2)}")

# --- build HGNC membership lookups -------------------------------------------
hgnc = pd.read_csv(ROOT / "inputs" / "hgnc_complete_set.txt",
                   sep="\t", dtype=str, low_memory=False)

def tags_for_row(group_field):
    if not isinstance(group_field, str):
        return set()
    return {GROUP2TAG[g.strip()] for g in group_field.split("|")
            if g.strip() in GROUP2TAG}

# tag lookups carry accessory-domain memberships; the *_known sets record which
# keys exist in HGNC at all, so we can tell "in HGNC, no accessory domain" apart
# from "not in HGNC".
ensg2tags, sym2tags, alias2tags = {}, {}, {}
ensg_known, sym_known, alias_known = set(), set(), set()
for _, r in hgnc.iterrows():
    tags = tags_for_row(r.get("gene_group"))
    ensg = r.get("ensembl_gene_id")
    if isinstance(ensg, str) and ensg:
        ensg_known.add(ensg)
        if tags:
            ensg2tags.setdefault(ensg, set()).update(tags)
    sym = r.get("symbol")
    if isinstance(sym, str) and sym:
        sym_known.add(sym.upper())
        if tags:
            sym2tags.setdefault(sym.upper(), set()).update(tags)
    for col in ("prev_symbol", "alias_symbol"):
        v = r.get(col)
        if isinstance(v, str) and v:
            for a in v.split("|"):
                if a.strip():
                    alias_known.add(a.strip().upper())
                    if tags:
                        alias2tags.setdefault(a.strip().upper(), set()).update(tags)

log(f"HGNC genes with a relevant group: by ENSG={len(ensg2tags)}, "
    f"by symbol={len(sym2tags)}, by prev/alias={len(alias2tags)}")

# --- match each C2H2 ZF TF ----------------------------------------------------
# Returns (accessory-domain tags, how-matched). 'how' records the key that hit
# HGNC; tags may be empty if the gene is in HGNC but has no accessory domain.
def match(ensg, sym):
    u = sym.upper()
    if isinstance(ensg, str) and ensg in ensg_known:
        return ensg2tags.get(ensg, set()), "ensembl_id"
    if u in sym_known:
        return sym2tags.get(u, set()), "symbol"
    if u in alias_known:
        return alias2tags.get(u, set()), "prev/alias_symbol"
    return set(), "not_in_hgnc"

rows = []
for _, r in c2h2.iterrows():
    tags, how = match(r["ensembl_id"], r["gene_symbol"])
    has = {t: (t in tags) for t in ("KRAB", "SCAN", "BTB")}
    sub = next((t for t in PRIORITY if has[t]), None)
    fine = f"{sub}-ZNF" if sub else "C2H2 ZF (other)"
    rows.append({
        "gene_symbol": r["gene_symbol"],
        "ensembl_id": r["ensembl_id"],
        "match_method": how,
        "in_hgnc_C2H2": "C2H2" in tags,
        "has_KRAB": has["KRAB"], "has_SCAN": has["SCAN"], "has_BTB": has["BTB"],
        "C2H2_subfamily": fine,
    })
det = pd.DataFrame(rows)
det.to_csv(INT / "c2h2_zf_subfamily_assignment.tsv", sep="\t", index=False)
log("wrote", INT / "c2h2_zf_subfamily_assignment.tsv")

log("match-method breakdown:")
for k, v in det["match_method"].value_counts().items():
    log(f"    {k}: {v}")
log("C2H2-type membership confirmed in HGNC for "
    f"{int(det['in_hgnc_C2H2'].sum())}/{len(det)} matched genes")

# --- fold fine labels back into the main paralog table ------------------------
sym2fine = dict(zip(det["gene_symbol"], det["C2H2_subfamily"]))
sym2flags = {r["gene_symbol"]: (r["has_KRAB"], r["has_SCAN"], r["has_BTB"])
             for _, r in det.iterrows()}

def fine_family(row):
    if row["DBD_family"] == "C2H2 ZF":
        return sym2fine.get(row["gene_symbol"], "C2H2 ZF (other)")
    return row["DBD_family"]

par["DBD_family_fine"] = par.apply(fine_family, axis=1)
for i, col in enumerate(("has_KRAB", "has_SCAN", "has_BTB")):
    par[col] = par["gene_symbol"].map(
        lambda s: sym2flags[s][i] if s in sym2flags else pd.NA)

# recompute family size on the fine grouping
fine_size = par.groupby("DBD_family_fine")["gene_symbol"].transform("nunique")
par["family_size_fine"] = fine_size

par = par[["gene_symbol", "ensembl_id", "DBD_family", "family_size",
           "DBD_family_fine", "family_size_fine",
           "has_KRAB", "has_SCAN", "has_BTB"]]
par.to_csv(RES / "paralogous_TFs.tsv", sep="\t", index=False)
log("updated", RES / "paralogous_TFs.tsv", "with fine-grained C2H2 subfamilies")

# --- regenerate the per-family summary on the FINE grouping -------------------
# Replaces the single 'C2H2 ZF' row from step 1 with one row per C2H2 subfamily
# (KRAB-ZNF, SCAN-ZNF, BTB-ZNF, C2H2 ZF (other)); all other families unchanged.
fam = par.copy()
fam["parent_DBD_family"] = fam["DBD_family"]
fine_summary = (
    fam.groupby("DBD_family_fine")
    .agg(
        parent_DBD_family=("parent_DBD_family",
                           lambda s: ", ".join(sorted(s.unique()))),
        n_members=("gene_symbol", "nunique"),
        members=("gene_symbol", lambda s: ", ".join(sorted(s.unique()))),
    )
    .sort_values("n_members", ascending=False)
    .reset_index()
    .rename(columns={"DBD_family_fine": "family"})
)
fine_summary.to_csv(RES / "paralogous_TF_family_summary.tsv", sep="\t", index=False)
log("regenerated", RES / "paralogous_TF_family_summary.tsv",
    "on fine grouping:", len(fine_summary), "families")
log("C2H2-derived rows now present:",
    ", ".join(sorted(det['C2H2_subfamily'].unique())))

# --- subfamily summary --------------------------------------------------------
summ = (det.groupby("C2H2_subfamily")
        .agg(n=("gene_symbol", "nunique"))
        .sort_values("n", ascending=False).reset_index())
summ.to_csv(RES / "C2H2_ZF_subfamily_summary.tsv", sep="\t", index=False)
log("")
log("=== C2H2 ZF split into fine subfamilies ===")
for _, r in summ.iterrows():
    log(f"    {r['C2H2_subfamily']:20s} {r['n']}")
log(f"    (genes carrying >=2 accessory domains counted once, by priority "
    f"{' > '.join(PRIORITY)})")
n_multi = int(((det[['has_KRAB','has_SCAN','has_BTB']].sum(axis=1)) >= 2).sum())
log(f"    genes with >=2 accessory domains: {n_multi}")
_logf.close()
