# inputs/ — data sources (not version-controlled)

Large reference inputs and soft-linked datasets are excluded from git (see `../.gitignore`).
To reproduce the pipeline, place the following here:

| File / dir | Source |
|---|---|
| `1-s2.0-S0092867418301065-mmc2.xlsx` | Lambert et al., *Cell* 2018, Table S1 (human TF catalog) |
| `gencode.v45.basic.annotation.gtf.gz` | GENCODE release 45, basic annotation (hg38) |
| `hgnc_complete_set.txt` | HGNC complete set (gene-group / accessory-domain annotation) |
| `human_cytokine_dict/human_cytokine_dict_mini.csv` | Oesinghaus et al. 2025 human cytokine dictionary (pseudobulk) |
| `pbmc/` | Human PBMC raw scRNA panel (only needed to rebuild the h5ad) |

Also required (read-only, not in repo):
- `results/intermediate/h5ad/human_pbmc.h5ad`, `human_spleen.h5ad` — only `var_names` are read.
- OHNOLOGS v2 human pair lists (`hsapiens_pairs_crit{0,A,C}.tsv`).

Small precomputed intermediates needed to run offline (DEG marker tables, label maps,
`ensembl_paralog_nodes.tsv`) ARE committed under `results/intermediate/`.
