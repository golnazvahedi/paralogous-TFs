#!/usr/bin/env python3
"""
Part I, Part 3: Separate clustered paralogs from singleton paralogs (hg38).

For each TF family (the FINE grouping DBD_family_fine, after CSD removal), find
genomic clusters of same-family members and split members into:
    results/clustered_paralogs.tsv   (clusters of >= 2 same-family members)
    results/singleton_paralogs.tsv   (everything else: mapped members not in a cluster)

Definitions (from CLAUDE.md Part 3)
-----------------------------------
1) Candidate genes come from the same family (paralogous_TF_family_summary.tsv /
   the DBD_family_fine column of paralogous_TFs.tsv).
2) Two same-family members A,B (A before B on the chromosome) are *adjacent* iff
   NO other protein-coding gene has its body strictly contained in the inter-gene
   gap (A.end, B.start). Other same-family members in the gap do NOT block;
   genes that merely overlap A or B (not strictly between) do NOT block. A cluster
   is a connected component of the per-family adjacency graph on one chromosome.
   (Testing consecutive members by genomic position is sufficient: if a blocker
   sits strictly between consecutive members it also sits strictly between any
   wider pair spanning them, so a broken consecutive edge truly splits clusters.)
3) Non-coding genes between members are tolerated; lncRNA genes whose body
   intersects any inter-member gap are listed in `lncRNAs_between`.
4) Promoter = TSS = strand-aware gene start from GENCODE v45 basic, hg38 primary
   assembly; recorded per member in `member_tss`.
5) `member_strands` records each member's strand, ordered by genomic position.

Coordinates come from GENCODE v45 basic (inputs/gencode.v45.basic.annotation.gtf.gz);
TF members are matched to GENCODE by versionless Ensembl gene ID, then by symbol.
"""
import gzip
import re
from bisect import bisect_right, bisect_left
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
INT = RES / "intermediate"
GTF = ROOT / "inputs" / "gencode.v45.basic.annotation.gtf.gz"
LOG = ROOT / "log" / "04_cluster_paralogs.log"
LOG.parent.mkdir(parents=True, exist_ok=True)

CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}

_logf = open(LOG, "w")
def log(*a):
    msg = " ".join(str(x) for x in a)
    print(msg); _logf.write(msg + "\n"); _logf.flush()

# ---------------------------------------------------------------------------
# 1. Parse GENCODE genes (cache to an intermediate TSV)
# ---------------------------------------------------------------------------
_re_id   = re.compile(r'gene_id "([^"]+)"')
_re_name = re.compile(r'gene_name "([^"]+)"')
_re_type = re.compile(r'gene_type "([^"]+)"')

def parse_gencode():
    genes = []
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene":
                continue
            chrom = f[0]
            if chrom not in CHROMS:
                continue
            attr = f[8]
            gid = _re_id.search(attr).group(1)
            genes.append({
                "ensg": gid.split(".")[0],
                "gene_id": gid,
                "name": _re_name.search(attr).group(1),
                "gtype": _re_type.search(attr).group(1),
                "chrom": chrom,
                "start": int(f[3]),           # GTF 1-based, start<=end
                "end": int(f[4]),
                "strand": f[6],
            })
    return pd.DataFrame(genes)

g = parse_gencode()
g.to_csv(INT / "gencode_v45_genes.tsv", sep="\t", index=False)
log(f"GENCODE genes parsed (main chroms): {len(g)}")
log("  protein_coding:", int((g.gtype == "protein_coding").sum()),
    "| lncRNA:", int((g.gtype == "lncRNA").sum()))

def tss(row):
    return row["start"] if row["strand"] == "+" else row["end"]

# per-chromosome protein-coding index for fast strict-containment queries
pc = g[g.gtype == "protein_coding"].sort_values(["chrom", "start"])
pc_by_chrom = {}
for chrom, sub in pc.groupby("chrom"):
    pc_by_chrom[chrom] = {
        "starts": sub["start"].tolist(),
        "ends": sub["end"].tolist(),
        "ensg": sub["ensg"].tolist(),
    }

# per-chromosome lncRNA index
lnc = g[g.gtype == "lncRNA"].sort_values(["chrom", "start"])
lnc_by_chrom = {}
for chrom, sub in lnc.groupby("chrom"):
    lnc_by_chrom[chrom] = {
        "starts": sub["start"].tolist(),
        "ends": sub["end"].tolist(),
        "name": sub["name"].tolist(),
    }

def coding_blocker(chrom, lo, hi, exclude_ensg):
    """True if some protein-coding gene (not in exclude_ensg) is strictly inside (lo,hi)."""
    if hi - lo <= 1:
        return False
    idx = pc_by_chrom.get(chrom)
    if idx is None:
        return False
    starts, ends, ensg = idx["starts"], idx["ends"], idx["ensg"]
    i0 = bisect_right(starts, lo)        # first gene with start > lo
    i1 = bisect_left(starts, hi)         # genes with start < hi
    for i in range(i0, i1):
        if ends[i] < hi and ensg[i] not in exclude_ensg:
            return True
    return False

def lncrnas_in_gap(chrom, lo, hi):
    """lncRNA gene names whose body intersects (overlaps) the open gap (lo,hi)."""
    if hi - lo <= 1:
        return []
    idx = lnc_by_chrom.get(chrom)
    if idx is None:
        return []
    starts, ends, names = idx["starts"], idx["ends"], idx["name"]
    out = []
    # overlap: g.start < hi and g.end > lo ; scan candidates with start < hi
    i1 = bisect_left(starts, hi)
    for i in range(i1):
        if ends[i] > lo:
            out.append(names[i])
    return out

# ---------------------------------------------------------------------------
# 2. Map paralogous TFs to GENCODE coordinates
# ---------------------------------------------------------------------------
par = pd.read_csv(RES / "paralogous_TFs.tsv", sep="\t", dtype=str)
log(f"paralogous TFs to place: {len(par)}  | families (fine): "
    f"{par['DBD_family_fine'].nunique()}")

by_ensg = {r.ensg: r for r in g.itertuples(index=False)}
by_name = {}
for r in g.itertuples(index=False):
    by_name.setdefault(r.name.upper(), r)   # first wins; protein_coding sorted first? not guaranteed

# prefer protein_coding when matching by name
by_name_pc = {}
for r in pc.itertuples(index=False):
    by_name_pc.setdefault(r.name.upper(), r)

members, unmapped = [], []
for _, row in par.iterrows():
    hit = by_ensg.get(row["ensembl_id"])
    how = "ensembl_id"
    if hit is None:
        hit = by_name_pc.get(row["gene_symbol"].upper()) or by_name.get(row["gene_symbol"].upper())
        how = "symbol" if hit is not None else None
    if hit is None:
        unmapped.append({"gene_symbol": row["gene_symbol"],
                         "ensembl_id": row["ensembl_id"],
                         "DBD_family_fine": row["DBD_family_fine"]})
        continue
    members.append({
        "gene_symbol": row["gene_symbol"],
        "ensembl_id": row["ensembl_id"],
        "DBD_family": row["DBD_family"],
        "family": row["DBD_family_fine"],
        "chrom": hit.chrom, "start": hit.start, "end": hit.end,
        "strand": hit.strand, "tss": hit.start if hit.strand == "+" else hit.end,
        "match": how,
    })
M = pd.DataFrame(members)
log(f"mapped to GENCODE: {len(M)}  (by ensembl_id "
    f"{int((M['match']=='ensembl_id').sum())}, by symbol "
    f"{int((M['match']=='symbol').sum())}) | unmapped: {len(unmapped)}")
if unmapped:
    pd.DataFrame(unmapped).to_csv(INT / "unmapped_members.tsv", sep="\t", index=False)
    log("  unmapped written to intermediate/unmapped_members.tsv:",
        ", ".join(u["gene_symbol"] for u in unmapped[:20]),
        "..." if len(unmapped) > 20 else "")

# ---------------------------------------------------------------------------
# 3. Cluster per family per chromosome
# ---------------------------------------------------------------------------
clusters = []          # one row per cluster (>=2 members)
singleton_rows = []    # mapped members not in a multi-member cluster
cluster_id = 0

for family, fam in M.groupby("family"):
    for chrom, ch in fam.groupby("chrom"):
        ch = ch.sort_values(["start", "end"]).reset_index(drop=True)
        if len(ch) == 1:
            singleton_rows.append(ch.iloc[0])
            continue
        # exclude-set: GENCODE ensg of this family's members on this chromosome
        excl = set()
        for _, m in ch.iterrows():
            hh = by_ensg.get(m["ensembl_id"])
            excl.add(hh.ensg if hh else m["ensembl_id"])
        # build consecutive adjacency, split into runs
        run = [0]
        runs = []
        for i in range(len(ch) - 1):
            A, B = ch.iloc[i], ch.iloc[i + 1]
            lo, hi = A["end"], B["start"]
            blocked = coding_blocker(chrom, lo, hi, excl)
            if blocked:
                runs.append(run); run = [i + 1]
            else:
                run.append(i + 1)
        runs.append(run)
        for idxs in runs:
            sub = ch.iloc[idxs].sort_values("start").reset_index(drop=True)
            if len(sub) == 1:
                singleton_rows.append(sub.iloc[0]); continue
            cluster_id += 1
            # lncRNAs across inter-member gaps
            lncs = []
            for i in range(len(sub) - 1):
                lo, hi = sub.iloc[i]["end"], sub.iloc[i + 1]["start"]
                lncs += lncrnas_in_gap(chrom, lo, hi)
            lncs = sorted(dict.fromkeys(lncs))   # dedup, keep order
            clusters.append({
                "cluster_id": f"clust{cluster_id:04d}",
                "family": family,
                "parent_DBD_family": ", ".join(sorted(sub["DBD_family"].unique())),
                "chrom": chrom,
                "n_members": len(sub),
                "members": ",".join(sub["gene_symbol"]),
                "member_strands": ",".join(sub["strand"]),
                "member_tss": ",".join(str(int(t)) for t in sub["tss"]),
                "member_starts": ",".join(str(int(s)) for s in sub["start"]),
                "member_ends": ",".join(str(int(e)) for e in sub["end"]),
                "cluster_start": int(sub["start"].min()),
                "cluster_end": int(sub["end"].max()),
                "span_bp": int(sub["end"].max() - sub["start"].min()),
                "n_lncRNAs_between": len(lncs),
                "lncRNAs_between": ",".join(lncs),
            })

# ---------------------------------------------------------------------------
# 4. Write outputs
# ---------------------------------------------------------------------------
cl = pd.DataFrame(clusters).sort_values(["family", "chrom", "cluster_start"])
cl.to_csv(RES / "clustered_paralogs.tsv", sep="\t", index=False)

sg = pd.DataFrame(singleton_rows)[
    ["gene_symbol", "ensembl_id", "DBD_family", "family",
     "chrom", "start", "end", "strand", "tss"]
].sort_values(["family", "chrom", "start"])
sg.to_csv(RES / "singleton_paralogs.tsv", sep="\t", index=False)

n_clustered_members = int(cl["n_members"].sum()) if len(cl) else 0
log("")
log("=== Part 3 result ===")
log("clusters (>=2 members):", len(cl))
log("clustered members:", n_clustered_members)
log("singleton members:", len(sg))
log(f"check: clustered + singleton = {n_clustered_members + len(sg)} "
    f"(mapped members = {len(M)})")
log("clusters with >=1 lncRNA between members:",
    int((cl["n_lncRNAs_between"] > 0).sum()) if len(cl) else 0)
log("")
log("families contributing the most clusters:")
if len(cl):
    for fam, n in cl["family"].value_counts().head(10).items():
        log(f"    {fam:20s} {n} clusters")
_logf.close()
