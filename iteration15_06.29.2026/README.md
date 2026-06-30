# iteration15 — Paralogous-TF duplication framework (clean, human-only)

A minimal, human-only distillation of the paralogous transcription-factor (TF)
duplication analysis. This iteration keeps only the load-bearing pieces of the
duplication story and drops the singleton-vs-clustered 2-group OR, the z-score DEG
thread, the RNA/cytokine co-expression and PBMC dosage threads, the GWAS/Mendelian
checks, the macaque cross-species checks, and all narrative reports.

Part of the broader **ETS P01 (PAR-27-082, NIAID)** program-project effort, which
develops the thesis that ETS-family TFs form a dosage-sensitive, paralog-buffered
regulatory network whose perturbation spans the immune-disease spectrum
(immunodeficiency ↔ autoinflammation ↔ autoimmunity).

---

## Central idea

Every paralogous TF is classified on a **2×2 of duplication PROVENANCE × ARRANGEMENT**:

|                | clustered (cis)         | dispersed              |
|----------------|-------------------------|------------------------|
| **2R-WGD ohnolog** | `clustered_cis_ohnolog` | `dispersed_ohnolog`    |
| **small-scale dup (SSD)** | `clustered_SSD_tandem`  | `dispersed_SSD`        |

This `dup_class4` label (in `results/TF_dup_2x2_classification.tsv`) is the unit of
all downstream analysis. We then ask which category is enriched for cell-type
**identity** DEGs vs **cytokine-response** DEGs, and which categories show dosage
**buffering** vs **amplification**.

Three design decisions define this iteration:

1. **Four-category framework** (the 2×2 above) supersedes the old per-gene 3-way
   `dup_mode`.
2. **KRAB-ZNF + HOX removed from both groups, up front** (inside `04`, before
   clustering). In iteration13 these silent large arrays were shown to drive most of
   the raw singleton > clustered enrichment, so every result here is the
   KRAB/HOX-controlled effect. Canonical catalogs carry the `.noKRABHOX` suffix.
3. **Human only** — human PBMC + spleen identity DEGs and the human cytokine
   dictionary. No mouse / macaque / lamprey.

---

## The four analyses (deliverables)

- **Part I** — generate the paralogous-TF catalogs (clustered vs singleton).
- **Part II** — classify every TF into the 4 duplication categories (`dup_class4`).
- **Part III** — DEG-enrichment odds-ratio (OR) across the 4 categories for: human
  PBMC identity markers, human spleen identity markers, and cytokine-response DEGs.
- **Part IV** — dosage-buffering analysis (by category) on cytokine-dictionary data.

---

## Headline findings (expected; reproduced from iteration14)

- **Identity (PBMC + spleen) enrichment is a property of DISPERSED ohnologs.**
  `dispersed_ohnolog` is the only FDR-significant enriched group (PBMC rate ~28.7%,
  OR ~1.79; spleen ~16.3%, OR ~2.87). **Clustering erases it** — within-ohnolog
  clustered-vs-dispersed OR ~0.47 (PBMC) / ~0.33 (spleen).
- **Cytokine-response enrichment follows ohnolog PROVENANCE and survives clustering.**
  `dispersed_ohnolog` most enriched (~37.5%, OR ~1.67, FDR-sig) but
  `clustered_cis_ohnolog` stays high (~29.2%, ns vs dispersed). Interpretation:
  clustering converts a 2R ohnolog from a candidate identity factor into a dedicated
  cytokine **effector** — the clustered ETS/STAT cis-ohnologs lose identity signal
  but keep cytokine response.
- **Dosage buffering is not a category-level trait** (every category median B ≈ 1.0),
  but among *neighbors* the two clustered categories diverge: SSD-tandem neighbors
  co-**amplify** (the chr2 SAND cluster SP100/110/140/140L, B 1.4–1.8) while
  cis-ohnologs harbor the cleanest **buffer**, STAT1/STAT4 (B=0.78, partial r=−0.41).

---

## Layout

```
inputs/        reference inputs (Lambert TF xlsx, GENCODE v45 gtf, HGNC set) +
               soft-linked large data (cytokine dict, PBMC raw panel)
scripts/       the 11 numbered pipeline scripts (01–11)
results/       output TSVs + figures/ ; intermediate/ holds caches & DEG markers
log/           run logs
reports/       LaTeX reports (none in this iteration)
```

> **Note:** `results/` currently holds only the precomputed *inputs*
> (`intermediate/deg`, `label_maps`, `ensembl_paralog_nodes.tsv`, and the `h5ad`
> symlink). Final TSVs and figures are produced by running the pipeline below.

---

## Pipeline

| Script | Part | Purpose | Key output |
|--------|------|---------|------------|
| `01_curate_paralogous_TFs.py`        | I  | Curate human TFs (Lambert `Is TF? == Yes`); family = shared DBD | `results/paralogous_TFs.tsv` |
| `02_subdivide_C2H2_ZF.py`            | I  | Subdivide C2H2-ZF by HGNC accessory domains (KRAB/SCAN/BTB) | — |
| `03_remove_CSD.py`                   | I  | Drop CSD RNA-binders | — |
| `04_cluster_paralogs.py`             | I  | Remove KRAB/HOX, split into clusters vs singletons | `results/clustered_paralogs.noKRABHOX.tsv`, `singleton_paralogs.noKRABHOX.tsv` |
| `05_TF_duplication_age.py`           | II | Annotate duplication MODE + AGE (OHNOLOGS v2 + Ensembl Compara cache; OFFLINE) | `results/TF_duplication_age.tsv` |
| `06_TF_dup_2x2_classification.py`    | II | Build 4-group `dup_class4` (provenance × arrangement; cis-ohnolog rescue) | `results/TF_dup_2x2_classification.tsv` |
| `07_human_pbmc_dup4_DEG_OR.py`       | III| PBMC identity-marker OR across 4 groups | `results/TF_dup4_DEG_OR.*` + figure |
| `08_human_spleen_dup4_DEG_OR.py`     | III| Spleen identity-marker OR across 4 groups | `results/TF_dup4_spleen_DEG_OR.*` |
| `09_human_cytokine_dict_standardize.py` | III | Standardize Oesinghaus 2025 cytokine pseudobulk (gate \|log_fc\|≥1 & adj_p<0.05) | `results/intermediate/human_cytokine_dict_DEGs.tsv` |
| `10_human_cytokine_dup4_DEG_OR.py`   | III| Cytokine-response DEG OR across 4 groups | `results/TF_dup4_cytokine_DEG_OR.*` + figure |
| `11_cytokine_dosage_buffering_by_category.py` | IV | Per-pair dosage buffering (B index + partial r), by category | `results/cytokine_dosage_buffering_by_category.*` + figure |

### How to run (from the iteration root)

```bash
PY=/mnt/alvand/apps/anaconda2/envs/py3/bin/python3
$PY scripts/01_curate_paralogous_TFs.py
$PY scripts/02_subdivide_C2H2_ZF.py
$PY scripts/03_remove_CSD.py
$PY scripts/04_cluster_paralogs.py
$PY scripts/05_TF_duplication_age.py            # OFFLINE via cached ensembl_paralog_nodes.tsv + OHNOLOGS
$PY scripts/06_TF_dup_2x2_classification.py     # 4-group catalog
$PY scripts/07_human_pbmc_dup4_DEG_OR.py
$PY scripts/08_human_spleen_dup4_DEG_OR.py
$PY scripts/09_human_cytokine_dict_standardize.py
$PY scripts/10_human_cytokine_dup4_DEG_OR.py
$PY scripts/11_cytokine_dosage_buffering_by_category.py
```

**Environment:** use `/mnt/alvand/apps/anaconda2/envs/py3/bin/python3` (scanpy 1.9.1,
anndata 0.8.0). Do **not** use the default `snakemake` env (broken sklearn). All
scripts resolve paths relative to the repo root via `Path(__file__).parents[1]`, so
run them from the iteration root.

---

## Cleaned catalog counts (n ≈ 1153)

| `dup_class4`          | n   |
|-----------------------|-----|
| dispersed_ohnolog     | 524 |
| dispersed_SSD         | 450 |
| clustered_cis_ohnolog | 106 |
| clustered_SSD_tandem  | 73  |

Reference cases: the ETS quartet (ETS1/FLI1/ETS2/ERG) and STAT1/3/4/5A/5B are
`clustered_cis_ohnolog`; the chr2 SAND cluster (SP100/110/140/140L) is
`clustered_SSD_tandem`.

---

## Provenance & data dependencies

Created 2026-06-29 from `../iteration14_06.27.26/`. **Code and small text outputs are
real copies; large data are soft links** resolved to their real targets on
`/mnt/alvand` (iteration11/13/7). The iteration survives iteration14 being moved but
still depends on iteration11/13/7.

Key soft links:
- `inputs/human_cytokine_dict/human_cytokine_dict_mini.csv` → iteration13 (Oesinghaus 2025).
- `inputs/pbmc/` → iteration11 raw scRNA panel (only needed to rebuild the h5ad).
- `results/intermediate/h5ad/` → **iteration11's REAL h5ad dir** (read-only; the OR
  scripts read only `var_names` for the detectable-gene universe).
- OHNOLOGS v2 human pairs read in place from `../iteration7_05.20.2026/inputs/ohnologs/`.

> ⚠️ **HAZARD:** `results/intermediate/h5ad/` is a *directory* symlink into a SHARED
> real dir. Never write files inside it — that lands in iteration11's real data and can
> clobber shared files. Reference its files read-only.

---

## Caveats

1. The OR denominator is the **detectable** paralogous-TF universe (catalog ∩ h5ad
   `var_names`), not an expression-matched background. The robust, FDR-backed results
   are the `dispersed_ohnolog` identity enrichment and the provenance-is-the-driver
   factorial; within-clustered contrasts and the interaction are **trends** (small n).
2. Spleen cell-type labels come from a blood-PBMC Azimuth transfer, so its cell-type
   universe matches PBMC (good for a reproducibility contrast, not added diversity).
   The cytokine dictionary is also PBMC — no tissue-resident subsets anywhere here.
3. Provenance depends on OHNOLOGS v2, which is conservative; the cis-ohnolog rescue in
   `06` fixes only *clustered* misses. A dispersed singleton OHNOLOGS miss (e.g. FEV)
   stays SSD. **Use `dup_class4`, never the per-gene `dup_mode`.**
