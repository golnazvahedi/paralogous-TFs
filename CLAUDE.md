iteration15 (2026-06-29) is the **clean, minimal, human-only** distillation of iteration14.
It keeps ONLY the load-bearing pieces of the paralogous-TF duplication story and drops
everything else (the singleton-vs-clustered 2-group OR, the z-score DEG thread, the
RNA/cytokine co-expression and PBMC dosage threads, the GWAS/Mendelian and macaque
cross-species checks, all reports).

THREE design decisions are carried verbatim from iteration14 and define this iteration:

  1. **Four-category duplication framework.** Every paralogous TF is classified on the 2x2
     of PROVENANCE (2R-WGD ohnolog vs SSD) x ARRANGEMENT (clustered vs dispersed):
        dispersed_ohnolog | clustered_cis_ohnolog | clustered_SSD_tandem | dispersed_SSD
     This is the `dup_class4` label in `results/TF_dup_2x2_classification.tsv`; it supersedes
     the old per-gene 3-way `dup_mode`.
  2. **KRAB-ZNF + HOX removed from BOTH groups, up front.** The exclusion happens inside
     `04_cluster_paralogs.py` BEFORE clustering (fine_family contains "KRAB"; HOX[ABCD]\d),
     so the canonical catalogs are `.noKRABHOX.tsv`. In iteration13 those silent large
     arrays were shown to drive most of the raw singleton>clustered enrichment, so every
     result here is the KRAB/HOX-controlled effect.
  3. **Human only.** PBMC + spleen identity DEGs and the human cytokine dictionary. No
     mouse / macaque / lamprey, no Cui mouse Immune-Dictionary.

The deliverables are exactly four analyses, all framed on `dup_class4`:
  Part I   — generate the paralogous-TF catalogs (clustered vs singleton).
  Part II  — classify every TF into the 4 duplication categories.
  Part III — DEG-enrichment odds-ratio (OR) across the 4 categories for: human PBMC
             identity markers, human spleen identity markers, human cytokine-response DEGs.
  Part IV  — dosage-buffering analysis (by category) on the cytokine-dictionary data.

############################ Provenance: set up 2026-06-29 from iteration14 ############################

Created from `../iteration14_06.27.26/`. **Code and small text outputs are real copies;
large data are soft links** back into the chain (iteration14 -> iteration13 -> iteration11,
all on `/mnt/alvand`). Links are resolved to their REAL targets here (not chained through
iteration14), so iteration15 survives even if iteration14 is moved — but it still depends on
iteration11/iteration13/iteration7 existing. If those are moved/deleted, re-point the links.

  - **Real copies (light):**
    - `scripts/` — the 11 scripts for the four parts (renumbered 01–11; mapping below).
    - reference inputs: Lambert TF xlsx (`1-s2.0-S0092867418301065-mmc2.xlsx`),
      GENCODE v45 basic gtf.gz, `hgnc_complete_set.txt`.
    - precomputed one-vs-rest DEG marker tables: human PBMC + human spleen cell-type
      markers under `results/intermediate/deg/*.top.tsv` (inputs to the OR scripts;
      produced in iteration13 by one-vs-rest Wilcoxon on the annotated h5ad — gate
      adj-p<0.05 & log2FC>=1 & pct_in>=0.25 & pct_out<=0.60).
    - label maps under `results/intermediate/label_maps/`.
    - `results/intermediate/ensembl_paralog_nodes.tsv` (Ensembl-Compara biomart cache) +
      `.pybiomart.sqlite`, so `05` runs OFFLINE (no biomart call needed).
  - **Soft links (large data), resolved to real targets:**
    - `inputs/human_cytokine_dict/human_cytokine_dict_mini.csv` -> iteration13 real
      (78 MB, Oesinghaus 2025).
    - `inputs/pbmc/` -> iteration11 real (raw scRNA panel; only needed to *rebuild* the
      h5ad, which we never do here).
    - `results/intermediate/h5ad/` -> **iteration11's REAL h5ad directory** (human_pbmc.h5ad
      5.7 GB + human_spleen.h5ad 6.4 GB; the OR scripts read only `var_names` for the
      detectable-gene universe, backed read).
      **HAZARD: this is a DIRECTORY symlink into a SHARED dir.** Writing a file *inside* it
      (e.g. `ln -sf X results/intermediate/h5ad/foo`) lands in iteration11's REAL directory
      and can clobber shared data. Do NOT create files under that dir; reference the
      existing files through it read-only.
    - OHNOLOGS v2 human pair lists are read in place from
      `../iteration7_05.20.2026/inputs/ohnologs/hsapiens_pairs_crit{0,A,C}.tsv` (no copy).

Script renumber map (iteration14 -> iteration15):
  04→04 nochange path, 24→05, 27→06, 28→07, 30→08, 07→09, 29→10, 32→11.

### Inputs and Outputs
- inputs -> `/inputs` | code -> `/scripts` | figures & spreadsheets -> `/results`
  (intermediates under `/results/intermediate`) | logs -> `/log` | reports (if any) in
  LaTeX compiled with `tectonic` -> `/reports`.

### Environment
- Python: `/mnt/alvand/apps/anaconda2/envs/py3/bin/python3` (scanpy 1.9.1, anndata 0.8.0;
  the default `snakemake` env has a broken sklearn — do not use it). Run scripts from the
  iteration root (`python3 scripts/NN_*.py`); all scripts resolve paths relative to the repo
  root via `Path(__file__).parents[1]`.

########################################## Part I — Paralogous-TF catalogs (clustered vs singleton) ##########################################

Scripts `01`–`04` (verbatim from iteration14). From Lambert et al. (Cell 2018, Table S1)
curate human TFs where `Is TF? == Yes`, define **family = shared DBD**, subdivide C2H2-ZF by
HGNC accessory domains (KRAB/SCAN/BTB), remove the CSD family. Then `04_cluster_paralogs.py`
places members on hg38 (GENCODE v45 basic), **removes KRAB-ZNF + HOX up front**, and splits
the survivors into **clusters vs singletons**: two same-family members are *adjacent* iff no
*other* protein-coding gene body lies strictly in the inter-gene gap (other same-family
members and lncRNAs do not block); clusters are connected components.

  - `01_curate_paralogous_TFs.py` -> `results/paralogous_TFs.tsv`,
    `results/paralogous_TF_family_summary.tsv`, `results/intermediate/lambert_isTF_yes.tsv`
    (the last is consumed by `06`).
  - `02_subdivide_C2H2_ZF.py`     -> C2H2-ZF subfamily assignment.
  - `03_remove_CSD.py`            -> drops CSD RNA-binders from the universe.
  - `04_cluster_paralogs.py`      -> the canonical catalogs:
      * `results/clustered_paralogs.noKRABHOX.tsv` (78 clusters / 184 members) +
        `results/singleton_paralogs.noKRABHOX.tsv` (979 singletons) — KRAB/HOX-free,
        inputs to Parts II–III.
      * `results/clustered_paralogs.full.tsv` + `results/singleton_paralogs.full.tsv` —
        full un-excluded catalog, audit only.

########################################## Part II — Four-category duplication classification ##########################################

  - `05_TF_duplication_age.py` (= iteration14 `24`). Annotates every paralogous TF with
    duplication MODE + AGE. Ohnolog (2R-WGD) calls from OHNOLOGS v2 human pair lists (reused
    from `../iteration7_05.20.2026/inputs/ohnologs/`, crit0/A/C). Finer age = Ensembl Compara
    youngest within-species paralog node (read from the cached
    `results/intermediate/ensembl_paralog_nodes.tsv` — runs OFFLINE; only hits biomart if the
    cache is missing). Output `results/TF_duplication_age.tsv` (per-gene `dup_mode` =
    WGD_ohnolog / tandem / SSD_dispersed).
    CAVEAT: OHNOLOGS v2 is conservative and has real misses — notably **FLI1/ERG/FEV are
    absent from all 3 crit levels** though FLI1/ERG is a textbook 2R pair. The per-gene
    `dup_mode` is therefore SUPERSEDED by the 2x2 `dup_class4` from `06`, which fixes the
    clustered misses via the general cis-ohnolog rescue rule. Prefer `dup_class4`.
  - `06_TF_dup_2x2_classification.py` (= iteration14 `27`). Builds the 4-group catalog by
    decoupling PROVENANCE x ARRANGEMENT. Provenance = OHNOLOGS relaxed pairs PLUS a GENERAL
    cis-ohnolog RESCUE rule: *a genomic cluster containing >=1 OHNOLOGS ohnolog is
    "ohnolog-bearing"; ALL its members get ohnolog provenance* (a tandem block copied by 2R
    carries silent partners through). Recovers FLI1/ERG with NO hand curation; leaves the
    chr2 SAND cluster as SSD-tandem (no ohnolog member). Cross-references Lambert (adds
    `TF_subfamily` = Lambert DBD) and removes 10 non-protein-coding entries by HGNC
    locus_type (pseudogene/readthrough/clone). Inputs: `TF_duplication_age.tsv` (from `05`),
    `clustered_paralogs.noKRABHOX.tsv` (from `04`), `hgnc_complete_set.txt`,
    `results/intermediate/lambert_isTF_yes.tsv` (from `01`). Outputs:
    `results/TF_dup_2x2_classification.tsv` (per-TF: provenance/arrangement/dup_class4/
    TF_subfamily/cluster_id/...), `.lists.tsv`, `.cis_ohnolog_rescued.tsv`, `_removed.tsv`.
    CLEANED COUNTS (n≈1153): dispersed_ohnolog 524, clustered_cis_ohnolog 106,
    clustered_SSD_tandem 73, dispersed_SSD 450. ETS quartet (ETS1/FLI1/ETS2/ERG) +
    STAT1/3/4/5A/5B = clustered_cis_ohnolog; chr2 SAND (SP100/110/140/140L) =
    clustered_SSD_tandem. KRAB/HOX already excluded upstream.

########################################## Part III — DEG-enrichment OR across the 4 categories ##########################################

Question: among paralogous TFs, which duplication category is enriched for cell-type
identity DEGs (PBMC, spleen) and for cytokine-response DEGs? For each category build a
2x2 (category vs rest) x (DEG / non-DEG) table -> one-vs-rest Fisher OR (+95% CI, p, BH-FDR);
plus the key PAIRWISE contrasts (provenance within arrangement, arrangement within
provenance, and the two marginals) and a provenance x arrangement FACTORIAL logit.

**Detectable-gene universe (UNIFIED in iteration15):** all three OR scripts intersect the
4-group catalog with the genes present in the relevant h5ad `var_names` (backed read).
This is the one analysis change vs iteration14: scripts `07`/`10` previously took the
detectable universe from the out-of-scope z-score-DEG table (`zscore_DEG_TF.per_TF.tsv`);
here they read `human_pbmc.h5ad var_names` directly, matching the spleen script `08`, so all
three share ONE universe rule and the z-score thread is fully dropped. Numbers may shift by
a hair vs iteration14's `28`/`29`, but the design and conclusions are unchanged.

  - `07_human_pbmc_dup4_DEG_OR.py` (= iteration14 `28`). Marker = canonical one-vs-rest PBMC
    lineage marker (`results/intermediate/deg/human_markers.major_lineage.top.tsv`).
    Universe = catalog ∩ `human_pbmc.h5ad` var_names. Outputs `results/TF_dup4_DEG_OR.{summary,
    by_celltype,per_TF}.tsv` + `results/figures/TF_dup4_DEG_OR.{pdf,png}`.
  - `08_human_spleen_dup4_DEG_OR.py` (= iteration14 `30`). Same as `07` but SPLEEN one-vs-rest
    markers (`human_spleen_markers.major_lineage.top.tsv`) and `human_spleen.h5ad` var_names.
    Spleen has 8 cell types (PBMC-Azimuth transfer). Outputs `results/TF_dup4_spleen_DEG_OR.*`.
  - `09_human_cytokine_dict_standardize.py` (= iteration14 `07`). Standardizes the Oesinghaus
    2025 cytokine-dictionary pseudobulk table (`inputs/human_cytokine_dict/
    human_cytokine_dict_mini.csv`, soft-linked): gate |log_fc|>=1 & adj_p<0.05,
    well_biased=False. Output `results/intermediate/human_cytokine_dict_DEGs.tsv` (consumed
    by `10`).
  - `10_human_cytokine_dup4_DEG_OR.py` (= iteration14 `29`). Marker = cytokine RESPONSE DEG
    (response to >=1 cytokine in >=1 cell type, from `09`). Universe = catalog ∩ PBMC
    var_names (cytokine dictionary is PBMC). Same stats as `07`. Outputs
    `results/TF_dup4_cytokine_DEG_OR.{summary,by_celltype,per_TF}.tsv` + figure.

HEADLINE (reproduced from iteration14, expected here):
  - **Identity (PBMC + spleen): enrichment is the property of DISPERSED ohnologs specifically.**
    dispersed_ohnolog is the only FDR-significant enriched group (PBMC rate ~28.7%, OR ~1.79;
    spleen ~16.3%, OR ~2.87). Clustering ERASES it: clustered_cis_ohnolog collapses to the
    SSD-tandem rate; within-ohnolog clustered-vs-dispersed OR ~0.47 (PBMC) / ~0.33 (spleen).
  - **Cytokine response: enrichment follows ohnolog PROVENANCE and SURVIVES clustering.**
    dispersed_ohnolog most enriched (~37.5%, OR ~1.67, FDR-sig) BUT clustered_cis_ohnolog
    stays high (~29.2%, ns vs dispersed). CONCLUSION: clustering converts a 2R ohnolog from a
    candidate identity factor into a dedicated cytokine EFFECTOR — the ETS (ETS2/ETV3) + STAT
    clustered cis-ohnologs lost identity but keep cytokine response.

########################################## Part IV — Cytokine-dictionary dosage buffering (by category) ##########################################

  - `11_cytokine_dosage_buffering_by_category.py` (= iteration14 `32`). DOSAGE BUFFERING
    across the 4 categories, per within-DBD-family PAIR, on the cytokine-dictionary pseudobulk
    log_fc table (`human_cytokine_dict_mini.csv`; absent member imputed 0). Two metrics:
    (1) buffering index B = Var(x+y)/(Var x + Var y) over union-active conditions (B<1 buffer,
    B>1 amplify); (2) common-mode-removed PARTIAL r (regress out the mean response of all
    4-category TFs; partial r<0 = compensation). Each pair tagged NEIGHBOR (same cluster) vs
    DISTAL, benchmarked against random cross-family-pair null. Inputs:
    `inputs/human_cytokine_dict/human_cytokine_dict_mini.csv`, `TF_dup_2x2_classification.tsv`.
    Outputs: `results/cytokine_dosage_buffering_by_category.{pairs,summary}.tsv` +
    `results/figures/cytokine_dosage_buffering_by_category.{pdf,png}`.
    RESULT (from iteration14): buffering is NOT a category-level trait (every category median
    B~1.0 = null), but the two CLUSTERED categories DIVERGE among neighbors — SSD-tandem
    neighbors co-AMPLIFY (B~1.23; SAND SP100/110/140/140L B 1.4–1.8, partial r up to 0.82)
    while cis-ohnolog harbor the cleanest BUFFER, STAT1/STAT4 (B=0.78, marginal r=-0.27 ->
    partial r=-0.41 after common-mode removal). cis-vs-SSD neighbor B contrast p~0.038. So
    buffering lives in specific cis-ohnolog pairs (STAT); amplification is the SSD-tandem
    (SAND) signature.

########################################## Extensions added 2026-06-30 (beyond the core four parts) ##########################################

These scripts were added on top of the canonical 01–11 pipeline; they READ its outputs and
do not alter the four-part results above. Numbered with letter suffixes to stay out of the
core sequence.

  - `05b_gene_origin_age_2R.py` — GENE-ORIGIN (phylostratigraphic) age of SELECTED TFs vs the
    two WGD rounds (R1/R2). For each gene, finds the deepest Ensembl-Compara ortholog on a
    species ladder straddling 2R (yeast/worm/fly/tunicate | R1 | hagfish/lamprey | R2 |
    elephant-shark…mouse) and brackets origin: before_R1 (invertebrate ortholog) /
    R1_to_R2 (cyclostome only) / after_R2 (gnathostome only). Distinct from script 05's
    `youngest_paralog_node`, which is PARALOG-duplication age, not gene origin. Args = gene
    symbols (default BCL6 IRF1 KLF2 — all came out before_R1, i.e. pre-2R ancient singletons
    retained single-copy). Outputs `results/TF_gene_origin_age_2R.{tsv,full per-species matrix}`,
    figure, cache `results/intermediate/gene_origin_orthologs.tsv`.
  - `05c_TF_age_distribution_by_category.py` — distribution of DUPLICATION age
    (`youngest_paralog_node`, binned into evolutionary epochs) across the 4 dup_class4
    categories. Outputs `results/TF_age_distribution_by_category.tsv` + figure. Result:
    categories differ (Kruskal-Wallis p=0.006); ohnologs most ancient, clustered_SSD_tandem
    youngest.
  - `05d_gene_origin_age_full_catalog.py` — full-catalog version of 05b: the ortholog-ladder
    phylostratigraphy run on ALL 1,153 TFs (13 chunked biomart queries; own cache
    `results/intermediate/gene_origin_orthologs.full.tsv`), then gene-origin window
    distribution per dup_class4. Outputs `results/TF_gene_origin_age_2R.full.tsv`,
    `results/TF_gene_origin_age_by_category.tsv`, figure. Result: chi-square p=2.2e-7;
    dispersed_ohnolog most pre-R1 (72%), clustered_SSD_tandem youngest (52% pre-R1, 30%
    after R2), cis-ohnolog most enriched in the R1–R2 window (24%).
  - `10b_cytokine_DEG_heatmaps_by_category.py` — cytokine-dictionary DEG heatmaps per
    dup_class4 in the iteration14 two-panel style (TF×cell-type breadth + TF×cytokine breadth,
    discrete-count viridis, TF-subfamily sidebar). Outputs
    `results/figures/TF_dup4_cytokine_DEG_{summary,heatmap.<cat>}.cytokine.*` +
    `results/TF_dup4_cytokine_DEG_{breadth,category_means,ranking}.cytokine.tsv`.

########################################## Part V — Alpha vs beta ohnolog DEG enrichment (Zhu et al. 2026) ##########################################

  - `12_alpha_beta_ohnolog_DEG_OR.py`. THE LAST ANALYSIS. Following Zhu et al. (Nature 2026,
    `inputs/Zhu_Nature_2026.pdf`): the jawed-vertebrate 2R WGD was an ALLOPOLYPLOIDIZATION
    fusing two parental lineages, ALPHA and BETA; alpha-derived ohnologs were ~4x more
    retained and (in BRAIN) more marker-associated. We reproduce the alpha/beta split on our
    ohnolog TFs and ask the project's question in IMMUNE cells.
    SOURCE/MAPPING: Marletaz hagfish paralogons table, downloaded to
    `inputs/marletaz_paralogons/Vert_Evt_OGrrA.txt` (col `1R` = alpha1/alpha2/beta1/beta2 per
    orthogroup per species). Collapse alpha1/2→alpha, beta1/2→beta; assign each human ohnolog
    TF a lineage by (1) CHICKEN gene-symbol orthology (Zhu's stated method) then (2) an
    orthogroup-level fallback for the unmapped (the two agree 100% where both apply).
    Universe = ohnolog-provenance TFs ∩ PBMC var_names. Markers = PBMC identity (script 07
    source) and cytokine-response (script 09/10 source). Stats: alpha & beta marker-rates,
    the ALPHA-vs-BETA Fisher OR (headline), and alpha/beta one-vs-rest vs the detectable
    paralogous-TF background; dispersed_SSD shown as reference.
    Inputs: `TF_dup_2x2_classification.tsv`, `inputs/marletaz_paralogons/Vert_Evt_OGrrA.txt`,
    `human_pbmc.h5ad` (var_names), `results/intermediate/deg/human_markers.major_lineage.top.tsv`,
    `results/intermediate/human_cytokine_dict_DEGs.tsv`. Outputs:
    `results/TF_ohnolog_alpha_beta.tsv` (per-TF lineage + flags),
    `results/TF_alpha_beta_DEG_OR.summary.tsv`, `results/figures/TF_alpha_beta_DEG_OR.{pdf,png}`.
    RESULT: 451/629 detectable ohnolog TFs labeled = 343 alpha / 108 beta (3.2:1, reproduces
    the alpha-retention bias). Zhu's OHNOLOG>SSD enrichment REPLICATES in PBMC (both lineages
    beat dispersed_SSD: identity alpha 26% / beta 32% vs SSD 19%; cytokine 36–37% vs 26%), but
    the ALPHA>BETA asymmetry DOES NOT — alpha-vs-beta ns for both modalities. Reversed split
    instead: beta edges IDENTITY (beta-vs-rest OR 1.69, p=0.022; alpha-vs-rest ns) while alpha
    edges CYTOKINE response (alpha-vs-rest OR 1.36, p=0.026; beta-vs-rest ns). The robust
    claim is ohnolog>SSD; the within-ohnolog alpha/beta contrasts are trends (beta n small).

########################################## Part VI — Monogenic immune-disease (IEI) gene enrichment ##########################################

  - `13_IEI_disease_gene_enrichment.py`. Intersects the 4-group catalog with curated MONOGENIC
    immune-disease genes (Inborn Errors of Immunity) and asks which duplication category is
    enriched. IEI set = Genomics England PanelApp, distilled to
    `inputs/iei/panelapp_iei_genes.tsv` (panel 398 "Primary immunodeficiency or monogenic IBD"
    v9.18 + panel 1075 "Autoinflammatory disorders" v3.10; confidence 3=green/diagnostic-grade
    = primary set, also a MONOALLELIC/dominant = haploinsufficiency-compatible green subset).
    Match by HGNC ID (catalog symbol→hgnc_id via `inputs/hgnc_complete_set.txt` incl.
    alias/prev; symbol fallback). Universe = full classified catalog (1153); one-vs-rest Fisher
    OR WITHIN that universe; merges alpha/beta lineage (12) + gene-origin window (05d).
    Inputs: `TF_dup_2x2_classification.tsv`, `inputs/iei/panelapp_iei_genes.tsv`,
    `inputs/hgnc_complete_set.txt`, `TF_ohnolog_alpha_beta.tsv`, `TF_gene_origin_age_2R.full.tsv`.
    Outputs: `results/TF_IEI_disease_genes.tsv` (per-TF annot + IEI flags + MOI),
    `results/TF_IEI_enrichment.summary.tsv`, `results/figures/TF_IEI_disease_enrichment.{pdf,png}`.
    RESULT: 42/1153 paralogous TFs are green IEI genes; enrichment is driven by PROVENANCE not
    buffering. dispersed_ohnolog OR=3.53 (p=2e-4; the same identity-TF reservoir as Part III);
    ohnolog-vs-SSD OR=4.34 (p=1e-4); pre-R1-ancient-vs-younger OR=3.40 (p=0.003); clustered-vs-
    dispersed (buffered vs unbuffered) ns (OR=0.73). Strongest for the monoallelic/dominant
    subset (dispersed_ohnolog OR=4.55). P01 tie-in: SPI1/PU.1 + ELF4 are dispersed_ohnolog green
    IEI genes; ETS1/FLI1 is the clustered buffered pair and is NOT monogenic (polygenic lupus
    GWAS) — buffering may convert a monogenic dosage gene into a polygenic risk locus.
    CAVEAT: enrichment is WITHIN the paralogous-TF universe (not genome-wide); IEI=monogenic
    only (polygenic autoimmunity is in the GWAS layer); PanelApp green is curated but evolving.

########################################## Part VII — Polygenic (GWAS) complex-immune-disease enrichment ##########################################

  - `14_gwas_curate_immune_disease_genes.py`. Curates GWAS genes for 5 COMPLEX immune diseases
    (SLE_lupus, rheumatoid_arthritis, allergy, asthma, type_1_diabetes) from the NHGRI-EBI GWAS
    Catalog full ontology-annotated associations, read IN PLACE from
    `../iteration7_05.20.2026/inputs/gwas_catalog_full.tsv` (project convention, not copied).
    Keeps genome-wide-significant rows (PVALUE_MLOG>=7.3 = P<=5e-8), matches diseases on
    MAPPED_TRAIT + DISEASE/TRAIT, explodes MAPPED_GENE (nearest/overlapping; intergenic SNPs
    give both flanking genes), flags is_MHC + is_protein_coding (via HGNC). Outputs
    `results/gwas_immune_disease_genes.long.tsv` (gene x disease), `.by_disease.tsv`, and a light
    copy `inputs/gwas/gwas_immune_disease_genes.tsv`. Protein-coding non-MHC gene counts: SLE 412,
    RA 444, allergy 1608, asthma 639, T1D 271. Top loci recover textbook biology (SLE STAT4/IRF5/
    ETS1; asthma IL33/IL1RL1/TSLP; T1D INS/IL2RA/PTPN22).
  - `15_gwas_TF_enrichment.py`. Intersects the curated GWAS genes (protein-coding non-MHC) with
    the 4-group catalog and computes one-vs-rest Fisher OR per dup_class4 (per disease + union),
    + the same headline contrasts as 13, merging lineage (12) / origin (05d) / IEI flag (13) for
    a MONOGENIC-vs-POLYGENIC contrast. Inputs: `gwas_immune_disease_genes.long.tsv`,
    `TF_dup_2x2_classification.tsv`, `inputs/hgnc_complete_set.txt`, `TF_ohnolog_alpha_beta.tsv`,
    `TF_gene_origin_age_2R.full.tsv`, `TF_IEI_disease_genes.tsv`. Outputs
    `results/TF_GWAS_disease_genes.tsv`, `results/TF_GWAS_enrichment.summary.tsv`,
    `results/figures/TF_GWAS_disease_enrichment.{pdf,png}`.
    RESULT (the monogenic->polygenic FLIP): 245/1153 (21%) paralogous TFs are GWAS genes.
    clustered_cis_ohnolog is the MOST GWAS-enriched category (29.2%, OR=1.61, p=0.045) despite
    being at BACKGROUND for monogenic IEI (3.8%) — buffering converts a monogenic dosage gene
    into a polygenic risk locus. dispersed_ohnolog also enriched (26%, OR=1.67); ohnolog-vs-SSD
    OR=2.06 (p=2e-6); dispersed_SSD depleted (OR=0.50). The GWAS-not-IEI clustered_cis_ohnolog
    TFs are the ETS archetypes ETS1/FLI1/ERG (+STAT5A/BATF/EHF/ETV3) — lupus/allergy GWAS loci,
    not Mendelian = the Vahedi autoimmunity pole. CAVEAT: nearest-gene mapping (not fine-mapped),
    MHC excluded, enrichment within the paralogous-TF universe.

### How to run (from the iteration root)
```
PY=/mnt/alvand/apps/anaconda2/envs/py3/bin/python3
$PY scripts/01_curate_paralogous_TFs.py
$PY scripts/02_subdivide_C2H2_ZF.py
$PY scripts/03_remove_CSD.py
$PY scripts/04_cluster_paralogs.py
$PY scripts/05_TF_duplication_age.py            # OFFLINE via cached ensembl_paralog_nodes.tsv + OHNOLOGS (../iteration7)
$PY scripts/06_TF_dup_2x2_classification.py     # 4-group catalog (needs HGNC + Lambert intermediate + OHNOLOGS via 05)
$PY scripts/07_human_pbmc_dup4_DEG_OR.py        # PBMC identity OR, 4 groups
$PY scripts/08_human_spleen_dup4_DEG_OR.py      # spleen identity OR, 4 groups
$PY scripts/09_human_cytokine_dict_standardize.py
$PY scripts/10_human_cytokine_dup4_DEG_OR.py    # cytokine-response OR, 4 groups
$PY scripts/11_cytokine_dosage_buffering_by_category.py   # cytokine dosage buffering, 4 groups
# --- extensions (2026-06-30; read core outputs, need internet for biomart in 05b/05d) ---
$PY scripts/05b_gene_origin_age_2R.py [GENE ...]           # gene-origin age vs R1/R2 (selected genes)
$PY scripts/05c_TF_age_distribution_by_category.py         # duplication-age distribution by category (offline)
$PY scripts/05d_gene_origin_age_full_catalog.py            # gene-origin age, all 1153 TFs (biomart; ~10 min)
$PY scripts/10b_cytokine_DEG_heatmaps_by_category.py [TOP_N]   # cytokine DEG heatmaps per category
$PY scripts/12_alpha_beta_ohnolog_DEG_OR.py               # Part V: alpha/beta ohnolog DEG OR (Zhu 2026)
$PY scripts/13_IEI_disease_gene_enrichment.py            # Part VI: monogenic immune-disease (IEI) enrichment
$PY scripts/14_gwas_curate_immune_disease_genes.py       # Part VII: curate GWAS genes (5 complex diseases; reads iteration7 catalog)
$PY scripts/15_gwas_TF_enrichment.py                     # Part VII: GWAS enrichment by category + monogenic-vs-polygenic flip
```

### CAVEATS to carry
  (i) The OR denominator is the DETECTABLE paralogous-TF universe (catalog ∩ h5ad var_names),
      not an expression-matched background. The robust, FDR-backed result is the
      dispersed_ohnolog identity enrichment + the provenance-is-the-driver factorial; the
      within-clustered contrasts and the interaction are TRENDS (the two clustered groups are
      small, n≈73–106 catalog / fewer detectable).
  (ii) Spleen cell-type labels come from a blood-PBMC Azimuth transfer, so its cell-type
      universe is the same lymphoid+monocyte lineages as PBMC (good for a PBMC-vs-spleen
      reproducibility contrast, not added cell-type diversity). The cytokine dictionary is
      also PBMC — no tissue-resident subsets anywhere in this iteration.
  (iii) Provenance depends on OHNOLOGS v2, which is conservative; the cis-ohnolog rescue
      (`06`) fixes only CLUSTERED misses. A DISPERSED singleton OHNOLOGS miss (e.g. FEV) stays
      SSD — accepted tradeoff for avoiding manual curation. Use `dup_class4`, never per-gene
      `dup_mode`.
  (iv) `results/intermediate/h5ad/` is a directory symlink into iteration11's SHARED, REAL
      h5ad dir — read-only; never write inside it.
