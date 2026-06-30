#!/usr/bin/env python3
"""
Curate paralogous TFs into the clean 2x2 of duplication PROVENANCE x genomic ARRANGEMENT,
decoupling the two axes that the per-gene `dup_mode` (script 05) wrongly collapsed.

  AXIS 1  provenance : ohnolog (2R-WGD, R1/R2)  vs  SSD (small-scale duplication)
  AXIS 2  arrangement: clustered (genomic array) vs dispersed (singleton)

  4 groups:
     dispersed_ohnolog      = ohnolog & dispersed   (the canonical 2R case)
     clustered_cis_ohnolog  = ohnolog & clustered   (ETS1/FLI1/ETS2/ERG: tandem within a
                                                      2R paralogon; each gene's WGD partner
                                                      is the TRANS copy elsewhere)
     clustered_SSD_tandem   = SSD     & clustered    (true local tandem, no WGD partner)
     dispersed_SSD          = SSD     & dispersed    (small-scale dispersed duplicate)

Provenance source: OHNOLOGS v2 relaxed PAIRS (script 05 `is_ohnolog_relaxed`) plus a
GENERAL, data-driven CIS-OHNOLOG RESCUE rule (replaces the earlier manual ETS list):

  RULE: a genomic cluster that contains >=1 OHNOLOGS ohnolog member is an "ohnolog-bearing"
  cluster; ALL of its members are then assigned ohnolog provenance. Rationale: a tandem
  block caught in a 2R whole-genome duplication is copied as a unit, so the silent tandem
  partners of a retained 2R ohnolog are themselves cis-ohnologs from the same event.

This recovers OHNOLOGS misses across ALL families with no hand-curation — e.g. FLI1/ERG
(their clusters hold ETS1/ETS2, which ARE OHNOLOGS ohnologs), exactly the case OHNOLOGS v2
drops. It is conservative where it should be: a true local-SSD array with no ohnolog
member (e.g. the chr2 SAND cluster SP100/110/140/140L) is NOT rescued and stays SSD-tandem.
LIMITATION: the rule rescues only CLUSTERED (cis) misses; a DISPERSED singleton ohnolog
that OHNOLOGS misses (e.g. FEV, the third ETS paralogon copy) has no cluster to rescue it
and remains SSD -- accepted, since avoiding manual curation is the explicit requirement.

Arrangement = our genomic clustered/singleton class (results catalogs, KRAB/HOX-free).

NOTE: this script ONLY curates and writes the lists; no enrichment/OR is run here.

Inputs : results/TF_duplication_age.tsv          (script 05; gene, cls, is_ohnolog_relaxed)
         results/clustered_paralogs.noKRABHOX.tsv (cluster_id <-> member map)
Outputs:
  results/TF_dup_2x2_classification.tsv         (per-TF: provenance, arrangement, dup_class4)
  results/TF_dup_2x2_classification.lists.tsv   (one row per group, comma-joined gene list)
  results/TF_dup_2x2_cis_ohnolog_rescued.tsv    (members rescued by the general rule, audit)
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/06_TF_dup_2x2_classification.py
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
DUP = RES / "TF_duplication_age.tsv"
CLUST = RES / "clustered_paralogs.noKRABHOX.tsv"
HGNC = ROOT / "inputs" / "hgnc_complete_set.txt"
OUT_PER = RES / "TF_dup_2x2_classification.tsv"
OUT_LISTS = RES / "TF_dup_2x2_classification.lists.tsv"
OUT_RESC = RES / "TF_dup_2x2_cis_ohnolog_rescued.tsv"
OUT_REMOVED = RES / "TF_dup_2x2_removed.tsv"
import re as _re
CLONE_RE = _re.compile(r"^(AC|AL|AP|AF|RP\d|CTD|CTC|CTA|LINC)\d|\.\d+$")

GROUP_ORDER = ["dispersed_ohnolog", "clustered_cis_ohnolog",
               "clustered_SSD_tandem", "dispersed_SSD"]


def log(*a):
    print("[27]", *a, flush=True)


def gene_to_cluster():
    cl = pd.read_csv(CLUST, sep="\t")
    m = {}
    for _, r in cl.iterrows():
        for g in str(r["members"]).split(","):
            m[g] = r["cluster_id"]
    return m


def lambert_map():
    """gene_symbol -> Lambert DBD (TF subfamily), and the set of Lambert Is-TF?==Yes symbols."""
    lam = pd.read_csv(RES / "intermediate" / "lambert_isTF_yes.tsv", sep="\t", dtype=str)
    return dict(zip(lam.gene_symbol, lam.DBD)), set(lam.gene_symbol)


def hgnc_locus_type():
    """symbol (incl prev/alias) -> HGNC locus_type, to flag readthrough/pseudogene."""
    h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
    m = dict(zip(h.symbol, h.locus_type))
    for col in ("prev_symbol", "alias_symbol"):
        for s_, lt_ in zip(h[col].fillna(""), h.locus_type):
            for s in s_.split("|"):
                if s:
                    m.setdefault(s, lt_)
    return m


def main():
    u = pd.read_csv(DUP, sep="\t")
    u["is_ohnolog_relaxed"] = u["is_ohnolog_relaxed"].astype(str).str.lower().isin(["true", "1"])

    # ---- AXIS 1: provenance with GENERAL cis-ohnolog rescue ----
    g2c = gene_to_cluster()
    u["cluster_id"] = u["gene"].map(g2c)
    # clusters that contain >=1 OHNOLOGS ohnolog member -> "ohnolog-bearing"
    ohno_clusters = set(u.loc[u.is_ohnolog_relaxed & u.cluster_id.notna(), "cluster_id"])
    in_ohno_cluster = u.cluster_id.isin(ohno_clusters) & u.cluster_id.notna()
    u["ohnolog_rescue_via_cluster"] = in_ohno_cluster & ~u.is_ohnolog_relaxed
    u["is_ohnolog"] = u.is_ohnolog_relaxed | in_ohno_cluster
    u["provenance"] = u["is_ohnolog"].map({True: "ohnolog", False: "SSD"})

    # ---- AXIS 2: arrangement (clustered vs dispersed) ----
    u["arrangement"] = u["cls"].map({"clustered": "clustered", "singleton": "dispersed"})

    # ---- 4-group label ----
    def grp(r):
        if r.provenance == "ohnolog":
            return "clustered_cis_ohnolog" if r.arrangement == "clustered" else "dispersed_ohnolog"
        return "clustered_SSD_tandem" if r.arrangement == "clustered" else "dispersed_SSD"
    u["dup_class4"] = u.apply(grp, axis=1)

    u = u.rename(columns={"cls": "genomic_class"})

    # ---- cross-reference Lambert 2018 (membership + TF subfamily = DBD) ----
    sym2dbd, lam_syms = lambert_map()
    u["in_lambert"] = u.gene.isin(lam_syms)
    u["TF_subfamily"] = u.gene.map(sym2dbd)

    # ---- flag readthroughs / pseudogenes / uncharacterized clone loci (HGNC) ----
    lt = hgnc_locus_type()
    u["hgnc_locus_type"] = u.gene.map(lt).fillna("NOT_IN_HGNC")

    def removal_reason(r):
        if not r.in_lambert:
            return "not_in_Lambert"
        if r.hgnc_locus_type == "readthrough":
            return "readthrough"
        if r.hgnc_locus_type == "pseudogene":
            return "pseudogene"
        if r.hgnc_locus_type == "NOT_IN_HGNC" and CLONE_RE.search(r.gene):
            return "uncharacterized_clone"
        return ""
    u["removal_reason"] = u.apply(removal_reason, axis=1)

    col_order = ["gene", "genomic_class", "arrangement", "provenance", "dup_class4",
                 "TF_subfamily", "in_lambert", "hgnc_locus_type", "is_ohnolog",
                 "is_ohnolog_relaxed", "ohnolog_rescue_via_cluster", "cluster_id",
                 "n_ohnolog_partners", "ohnolog_dup_time", "youngest_paralog_node",
                 "age_category"]
    col_order = [c for c in col_order if c in u.columns]

    removed = u[u.removal_reason != ""][["gene", "dup_class4", "TF_subfamily",
                                         "hgnc_locus_type", "removal_reason"]] \
        .sort_values("removal_reason").reset_index(drop=True)
    kept = u[u.removal_reason == ""][col_order] \
        .sort_values(["dup_class4", "TF_subfamily", "gene"]).reset_index(drop=True)

    kept.to_csv(OUT_PER, sep="\t", index=False)
    removed.to_csv(OUT_REMOVED, sep="\t", index=False)

    # wide lists (kept only)
    rows = []
    for g in GROUP_ORDER:
        genes = sorted(kept.loc[kept.dup_class4 == g, "gene"])
        rows.append(dict(dup_class4=g, n_TFs=len(genes), genes=",".join(genes)))
    pd.DataFrame(rows).to_csv(OUT_LISTS, sep="\t", index=False)

    # audit: members rescued to cis-ohnolog (kept only)
    resc = kept[kept.ohnolog_rescue_via_cluster][["gene", "cluster_id", "TF_subfamily", "dup_class4"]]
    resc.to_csv(OUT_RESC, sep="\t", index=False)

    log(f"universe before cleaning: {len(u)} TFs")
    log(f"removed {len(removed)} (readthrough/pseudogene/clone): "
        f"{', '.join(removed.gene.tolist())}")
    print(removed.to_string(index=False))
    log(f"all kept are Lambert 2018 Is-TF?==Yes: {bool(kept.in_lambert.all())}")
    log(f"cleaned universe: {len(kept)} paralogous TFs")
    log("2x2 group counts (cleaned):")
    print(kept.dup_class4.value_counts().reindex(GROUP_ORDER).to_string())
    log(f"cis-ohnolog rescued by general cluster rule: {len(resc)} members")
    log("ETS paralogon check:")
    print(kept[kept.gene.isin(['ETS1', 'ETS2', 'FLI1', 'ERG'])]
          [["gene", "TF_subfamily", "provenance", "dup_class4"]].to_string(index=False))
    log(f"wrote {OUT_PER} / {OUT_LISTS} / {OUT_RESC} / {OUT_REMOVED}")


if __name__ == "__main__":
    main()
