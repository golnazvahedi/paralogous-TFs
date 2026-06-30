#!/usr/bin/env python3
"""
Part I, step 3: Remove CSD (cold-shock domain) genes from the paralog catalog.

The CSD 'DBD family' in Lambert (LIN28A/B, YBX1/2/3) comprises cold-shock-domain
proteins that bind RNA / single-stranded nucleic acid rather than acting as
bona fide sequence-specific DNA-binding TFs (LIN28A is the flagged example).
They are dropped from the paralogous-TF catalog and all downstream tables.

This step is idempotent: it removes any row whose DBD_family is in EXCLUDE_FAMILIES
and regenerates the per-family summary, so re-running 01 -> 02 -> 03 reproduces
the cleaned catalog.

Inputs / outputs (all under results/)
  paralogous_TFs.tsv                      (rewritten, CSD rows removed)
  intermediate/lambert_isTF_yes.tsv       (rewritten, CSD rows removed)
  paralogous_TF_family_summary.tsv        (regenerated on fine grouping)
  intermediate/removed_genes.tsv          (audit of what was dropped)
  log/03_remove_CSD.log
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
INT = RES / "intermediate"
LOG = ROOT / "log" / "03_remove_CSD.log"
LOG.parent.mkdir(parents=True, exist_ok=True)

EXCLUDE_FAMILIES = {"CSD"}   # cold-shock-domain, RNA-binding -> not DNA-binding TFs

_logf = open(LOG, "w")
def log(*a):
    msg = " ".join(str(x) for x in a)
    print(msg); _logf.write(msg + "\n"); _logf.flush()

log("# Step 3: drop non-TF DBD families:", ", ".join(sorted(EXCLUDE_FAMILIES)))

# --- main paralog table -------------------------------------------------------
par = pd.read_csv(RES / "paralogous_TFs.tsv", sep="\t", dtype=str)
drop = par[par["DBD_family"].isin(EXCLUDE_FAMILIES)].copy()
keep = par[~par["DBD_family"].isin(EXCLUDE_FAMILIES)].copy()
log(f"paralogous_TFs.tsv: {len(par)} -> {len(keep)} (removed {len(drop)})")
log("removed genes:", ", ".join(drop["gene_symbol"].tolist()))

# removing whole families does not change any *other* family's size, so the
# pre-computed family_size / family_size_fine columns remain correct.
keep.to_csv(RES / "paralogous_TFs.tsv", sep="\t", index=False)
drop.assign(reason="CSD cold-shock-domain RNA-binding, not a DNA-binding TF") \
    .to_csv(INT / "removed_genes.tsv", sep="\t", index=False)

# --- intermediate curated list ------------------------------------------------
lam = pd.read_csv(INT / "lambert_isTF_yes.tsv", sep="\t", dtype=str)
# this file uses the original DBD column name
dbd_col = "DBD" if "DBD" in lam.columns else "DBD_family"
n0 = len(lam)
lam = lam[~lam[dbd_col].isin(EXCLUDE_FAMILIES)].copy()
lam.to_csv(INT / "lambert_isTF_yes.tsv", sep="\t", index=False)
log(f"lambert_isTF_yes.tsv: {n0} -> {len(lam)}")

# --- regenerate per-family summary on the fine grouping -----------------------
fam = keep.copy()
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
log(f"paralogous_TF_family_summary.tsv: {len(fine_summary)} family rows")
log("")
log("=== cleaned catalog ===")
log("paralogous TFs remaining:", len(keep))
log("DBD families (coarse):", keep["DBD_family"].nunique())
log("families (fine grouping):", len(fine_summary))
log("confirm CSD gone:", "CSD" not in set(keep["DBD_family"]))
_logf.close()
