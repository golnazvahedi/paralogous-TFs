#!/usr/bin/env python3
"""
Curate GWAS-based genes for COMPLEX (polygenic) immune diseases -- the autoimmunity/allergy
counterpart to the monogenic-IEI curation (script 13). This is the polygenic layer the IEI
panel cannot capture (e.g. the ETS1/FLI1 lupus locus).

Source = NHGRI-EBI GWAS Catalog, full ontology-annotated associations, read IN PLACE from
iteration7 (project convention, like OHNOLOGS -- not copied):
  ../iteration7_05.20.2026/inputs/gwas_catalog_full.tsv  (1.12M associations, v current).

Five disease groups, matched case-insensitively against the EFO-normalized MAPPED_TRAIT plus
the author DISEASE/TRAIT string:
  SLE_lupus | rheumatoid_arthritis | allergy | asthma | type_1_diabetes
Associations are kept only at GENOME-WIDE SIGNIFICANCE (P <= 5e-8, i.e. PVALUE_MLOG >= 7.3).
Genes are taken from the catalog MAPPED_GENE field (nearest/overlapping gene annotation;
intergenic SNPs contribute both flanking genes), split on the catalog's separators. HLA/MHC
genes are flagged (not dropped) so the MHC can be excluded downstream.

Outputs:
  results/gwas_immune_disease_genes.long.tsv     one row per (gene, disease): n_assoc, n_studies,
                                                 n_snps, max_mlog, min_p, is_MHC
  results/gwas_immune_disease_genes.by_disease.tsv  per-disease summary (n_assoc, n_genes)
  inputs/gwas/gwas_immune_disease_genes.tsv      light curated copy of the long table (reusable)
Run: /mnt/alvand/apps/anaconda2/envs/py3/bin/python3 scripts/14_gwas_curate_immune_disease_genes.py
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
HGNC = ROOT / "inputs" / "hgnc_complete_set.txt"
CATALOG = ROOT.parent / "iteration7_05.20.2026" / "inputs" / "gwas_catalog_full.tsv"
OUT_LONG = RES / "gwas_immune_disease_genes.long.tsv"
OUT_SUM = RES / "gwas_immune_disease_genes.by_disease.tsv"
OUT_INPUTS = ROOT / "inputs" / "gwas" / "gwas_immune_disease_genes.tsv"

MLOG_THR = 7.30103          # -log10(5e-8), genome-wide significance

# disease -> (include regex, exclude regex or None). Matched on lowercased trait text.
DISEASES = {
    "SLE_lupus":           (r"lupus", None),
    "rheumatoid_arthritis": (r"rheumatoid arthritis", None),
    "allergy":             (r"allerg|atopic|eczema|hay fever", None),
    "asthma":              (r"asthma", None),
    "type_1_diabetes":     (r"type 1 diabetes|type i diabetes", r"type 2|type ii"),
}
USECOLS = ["DISEASE/TRAIT", "MAPPED_TRAIT", "MAPPED_GENE", "SNPS", "P-VALUE",
           "PVALUE_MLOG", "STUDY ACCESSION"]
GENE_SPLIT = re.compile(r"\s+-\s+|[,;/]\s*")          # ' - ' intergenic, ',' ';' '/' multi
SKIP_GENES = {"", "NR", "INTERGENIC", "NA"}


def log(*a):
    print("[14]", *a, flush=True)


def protein_coding_symbols():
    """UPPER(symbol|alias|prev) -> True for HGNC protein-coding genes."""
    h = pd.read_csv(HGNC, sep="\t", dtype=str, low_memory=False)
    h = h[h["locus_group"] == "protein-coding gene"]
    pc = set()
    for _, r in h.iterrows():
        for col in ("symbol", "alias_symbol", "prev_symbol"):
            v = r.get(col)
            if isinstance(v, str) and v:
                pc.update(s.upper() for s in v.split("|"))
    return pc


def parse_genes(s):
    if not isinstance(s, str) or not s.strip():
        return []
    out = []
    for g in GENE_SPLIT.split(s):
        g = g.strip()
        if g.upper() not in SKIP_GENES and not g.startswith("LOC"):
            out.append(g)
    return out


def main():
    log(f"reading {CATALOG} (selected columns) ...")
    df = pd.read_csv(CATALOG, sep="\t", usecols=USECOLS, dtype=str,
                     engine="c", on_bad_lines="skip", low_memory=False)
    df["mlog"] = pd.to_numeric(df["PVALUE_MLOG"], errors="coerce")
    df = df[df["mlog"] >= MLOG_THR].copy()
    log(f"genome-wide-significant associations (P<=5e-8): {len(df)}")

    trait = (df["MAPPED_TRAIT"].fillna("") + " || " + df["DISEASE/TRAIT"].fillna("")).str.lower()

    rows = []
    dsum = []
    for dis, (inc, exc) in DISEASES.items():
        m = trait.str.contains(inc, regex=True, na=False)
        if exc:
            m &= ~trait.str.contains(exc, regex=True, na=False)
        sub = df[m]
        dsum.append(dict(disease=dis, n_associations=int(m.sum()),
                         n_studies=sub["STUDY ACCESSION"].nunique()))
        log(f"  {dis}: {int(m.sum())} assoc, {sub['STUDY ACCESSION'].nunique()} studies")
        # explode to genes
        per = {}
        for _, r in sub.iterrows():
            for g in parse_genes(r["MAPPED_GENE"]):
                d = per.setdefault(g, dict(n_assoc=0, snps=set(), studies=set(), mlog=[]))
                d["n_assoc"] += 1
                d["snps"].add(r["SNPS"]); d["studies"].add(r["STUDY ACCESSION"])
                d["mlog"].append(r["mlog"])
        for g, d in per.items():
            rows.append(dict(gene=g, disease=dis, n_assoc=d["n_assoc"],
                             n_snps=len(d["snps"]), n_studies=len(d["studies"]),
                             max_mlog=max(d["mlog"]), min_p=10 ** (-max(d["mlog"]))))

    long = pd.DataFrame(rows).sort_values(["disease", "n_assoc"], ascending=[True, False])
    long["is_MHC"] = long.gene.str.upper().str.startswith("HLA-") | long.gene.str.upper().isin(
        {"MICA", "MICB", "C4A", "C4B", "TNF", "NOTCH4", "AGER", "HLA-DRA"})
    pc = protein_coding_symbols()
    long["is_protein_coding"] = long.gene.str.upper().isin(pc)
    long.to_csv(OUT_LONG, sep="\t", index=False)
    OUT_INPUTS.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT_INPUTS, sep="\t", index=False)

    sumdf = pd.DataFrame(dsum)
    sumdf["n_genes"] = [long[long.disease == d].gene.nunique() for d in sumdf.disease]
    sumdf["n_genes_proteincoding_noMHC"] = [
        long[(long.disease == d) & long.is_protein_coding & (~long.is_MHC)].gene.nunique()
        for d in sumdf.disease]
    sumdf.to_csv(OUT_SUM, sep="\t", index=False)

    clean = long[long.is_protein_coding & (~long.is_MHC)]
    log("per-disease summary:")
    print(sumdf.to_string(index=False))
    log("union of curated GWAS genes: %d total | %d protein-coding non-MHC" %
        (long.gene.nunique(), clean.gene.nunique()))
    log("top protein-coding non-MHC genes per disease (by n_assoc):")
    for d in DISEASES:
        top = clean[clean.disease == d].nlargest(15, "n_assoc")["gene"].tolist()
        print(f"  {d}: {', '.join(top)}")
    log(f"wrote {OUT_LONG} / {OUT_SUM} / {OUT_INPUTS}")


if __name__ == "__main__":
    main()
