#!/usr/bin/env python3
"""
Part III -- HUMAN cytokine dictionary (Oesinghaus et al. 2025, bioRxiv;
theislab/HumanCytokineDict). Standardize the published pseudobulk DEG table into the
schema used by the OR script (08), so singleton-vs-clustered can be tested DIRECTLY in
human with NO mouse->human ortholog assignment.

iteration14 clean human-only run. Provenance: iteration13 script 31 (BASE made
relative to the repo root).

Source: inputs/human_cytokine_dict/human_cytokine_dict_mini.csv
  (per celltype x cytokine x gene pseudobulk DE: log_fc, logCPM, F, adj_p_value; HUMAN
   symbols; 12 donors; well_biased already filtered out in the 'mini' table.)
Gate: |log_fc| >= 1 & adj_p_value < 0.05, well_biased=False.

Output: results/intermediate/human_cytokine_dict_DEGs.tsv
  cols: sheet, Celltype, Cytokine, Gene, human_symbol, map_method, Avg_log2FC, FDR, direction
Run: python3 scripts/07_human_cytokine_dict_standardize.py
"""
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CSV  = ROOT / "inputs" / "human_cytokine_dict" / "human_cytokine_dict_mini.csv"
OUT  = ROOT / "results" / "intermediate" / "human_cytokine_dict_DEGs.tsv"
LFC, PADJ = 1.0, 0.05


def log(*a): print("[09]", *a, flush=True)


def main():
    d = pd.read_csv(CSV, index_col=0)
    log(f"rows={len(d)}, celltypes={d.celltype.nunique()}, cytokines={d.cytokine.nunique()}")
    d = d[(~d.well_biased) & (d.adj_p_value < PADJ) & (d.log_fc.abs() >= LFC)].copy()
    out = pd.DataFrame({
        "sheet": d.celltype.astype(str),
        "Celltype": d.celltype.astype(str),
        "Cytokine": d.cytokine.astype(str),
        "Gene": d.gene.astype(str),          # already human symbol
        "human_symbol": d.gene.astype(str),
        "map_method": "human_native",
        "Avg_log2FC": d.log_fc,
        "FDR": d.adj_p_value,
        "direction": np.where(d.log_fc > 0, "up", "down"),
    })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    log(f"wrote {OUT}: {len(out)} gated DEGs "
        f"({(out.direction=='up').sum()} up / {(out.direction=='down').sum()} down)")
    log("DEGs per cell type:")
    print(out.groupby("sheet").size().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
