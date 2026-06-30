#!/usr/bin/env python3
"""
Part I, step 1: Curate paralogous transcription factors from Lambert et al. 2018.

Source: inputs/1-s2.0-S0092867418301065-mmc2.xlsx
        tab 'Table S1. Related to Figure 1B' (Lambert et al., Cell 2018, human TF census)

Logic
-----
1. Read Table S1. The real column header is on the 2nd row of the sheet
   (row index 1); data starts on row index 2. Columns used:
        ID       -> Ensembl gene ID
        Name     -> gene symbol
        DBD      -> DNA-binding domain family
        Is TF?   -> curators' final TF assessment (Yes / No / ...)
2. Keep only entries with Is TF? == 'Yes'.
3. Define a "family" as all TFs sharing the same DBD value. Two or more TFs
   in the same DBD family are considered paralogous TFs.
4. Emit:
   - results/intermediate/lambert_isTF_yes.tsv : all curated TFs (Is TF? == Yes)
   - results/paralogous_TFs.tsv                : TFs that belong to a family of
                                                 size >= 2 (i.e. paralogous), with
                                                 their DBD family + family size.
   - results/paralogous_TF_family_summary.tsv  : per-family member counts.

Logs go to log/01_curate_paralogous_TFs.log
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "inputs" / "1-s2.0-S0092867418301065-mmc2.xlsx"
SHEET = "Table S1. Related to Figure 1B"
RES = ROOT / "results"
INT = RES / "intermediate"
LOG = ROOT / "log" / "01_curate_paralogous_TFs.log"

INT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)

_logf = open(LOG, "w")
def log(*a):
    msg = " ".join(str(x) for x in a)
    print(msg)
    _logf.write(msg + "\n")
    _logf.flush()

log("# Curate paralogous TFs from Lambert et al. 2018 Table S1")
log("source:", XLSX)

# --- read raw, take the second physical row (index 1) as the header ------------
raw = pd.read_excel(XLSX, sheet_name=SHEET, header=None, dtype=str)
log("raw sheet shape:", raw.shape)

# The sheet has a two-row header: row 0 holds group labels (e.g. 'Is TF?'),
# row 1 holds the per-column names (e.g. 'ID','Name','DBD'). Merge them so a
# label is matched whether it sits on row 0 or row 1.
header0 = raw.iloc[0].tolist()
header1 = raw.iloc[1].tolist()
def col_idx(label):
    for i in range(len(header1)):
        for h in (header1[i], header0[i]):
            if isinstance(h, str) and h.strip().lower() == label.lower():
                return i
    raise KeyError(f"header {label!r} not found")

i_id   = col_idx("ID")
i_name = col_idx("Name")
i_dbd  = col_idx("DBD")
i_istf = col_idx("Is TF?")
log(f"column indices -> ID={i_id} Name={i_name} DBD={i_dbd} Is TF?={i_istf}")

df = raw.iloc[2:, [i_id, i_name, i_dbd, i_istf]].copy()
df.columns = ["ensembl_id", "gene_symbol", "DBD", "is_TF"]
df = df.dropna(subset=["gene_symbol"])
for c in df.columns:
    df[c] = df[c].astype(str).str.strip()
log("total Lambert entries (rows with a gene symbol):", len(df))
log("Is TF? value counts:")
for k, v in df["is_TF"].value_counts(dropna=False).items():
    log(f"    {k!r}: {v}")

# --- step 2: keep only curated TFs --------------------------------------------
tfs = df[df["is_TF"].str.lower() == "yes"].copy()
log("entries with Is TF? == 'Yes':", len(tfs))

# Drop entries with no usable DBD family label.
bad_dbd = {"", "nan", "none", "unknown"}
tfs["DBD"] = tfs["DBD"].fillna("")
n_before = len(tfs)
tfs_named = tfs[~tfs["DBD"].str.lower().isin(bad_dbd)].copy()
log(f"of those, with a usable DBD family label: {len(tfs_named)} "
    f"(dropped {n_before - len(tfs_named)} with missing/Unknown DBD)")

tfs_named = tfs_named.drop_duplicates(subset=["gene_symbol"]).sort_values(
    ["DBD", "gene_symbol"]
)
tfs_named.to_csv(INT / "lambert_isTF_yes.tsv", sep="\t", index=False)
log("wrote", INT / "lambert_isTF_yes.tsv")

# --- step 3: family sizes; paralogous = family size >= 2 ----------------------
fam_size = tfs_named.groupby("DBD")["gene_symbol"].nunique().rename("family_size")
tfs_named = tfs_named.merge(fam_size, left_on="DBD", right_index=True)

paralogous = tfs_named[tfs_named["family_size"] >= 2].copy()
paralogous = paralogous.rename(columns={"DBD": "DBD_family"})
paralogous = paralogous[
    ["gene_symbol", "ensembl_id", "DBD_family", "family_size"]
].sort_values(["DBD_family", "gene_symbol"])
paralogous.to_csv(RES / "paralogous_TFs.tsv", sep="\t", index=False)
log("wrote", RES / "paralogous_TFs.tsv")

# --- family summary -----------------------------------------------------------
summary = (
    paralogous.groupby("DBD_family")
    .agg(
        n_members=("gene_symbol", "nunique"),
        members=("gene_symbol", lambda s: ", ".join(sorted(s.unique()))),
    )
    .sort_values("n_members", ascending=False)
    .reset_index()
)
summary.to_csv(RES / "paralogous_TF_family_summary.tsv", sep="\t", index=False)
log("wrote", RES / "paralogous_TF_family_summary.tsv")

# --- headline numbers ---------------------------------------------------------
n_singleton = int((tfs_named["family_size"] == 1).sum())
log("")
log("=== summary ===")
log("curated TFs (Is TF?=Yes, DBD known):", len(tfs_named))
log("paralogous TFs (family size >= 2):  ", len(paralogous))
log("singleton TFs (family size == 1):   ", n_singleton)
log("number of multi-member DBD families:", len(summary))
log("largest families:")
for _, r in summary.head(10).iterrows():
    log(f"    {r['DBD_family']}: {r['n_members']}")
_logf.close()
