############################ Provenance: set up 2026-07-17 from iteration16 ############################

Created from `../iteration16_07.02.2026/`. The purpose of this iteration is FOCUSED: keep the
four-group TF duplication catalog EXACTLY as built in iteration16, and add ONE new analysis —
whether CTCF binding shows distinct patterns between the two CLUSTERED groups.

**Reuse (do NOT recompute):**
  - The four-group catalog `results/TF_dup_2x2_classification.tsv` (per-TF `dup_class4`,
    `cluster_id`, `TF_subfamily`) and the clustered-locus coordinates
    `results/clustered_paralogs.noKRABHOX.tsv` (`cluster_id`, `chrom`, `cluster_start`,
    `cluster_end`, `span_bp`, `members`, `member_strands`, `member_tss`). Genome build = **hg38**
    (GENCODE v45 basic), same as iteration16.
  - Copy (or soft-link) these two files from `../iteration16_07.02.2026/results/` into this
    iteration's `inputs/` so the CTCF analysis runs standalone. Scripts `01`–`06` that BUILD the
    catalog do not need to be re-run; if a full rebuild is ever wanted, they live verbatim in
    iteration16.

The four duplication groups (`dup_class4`), kept AS-IS, with iteration16 counts (n=1153):
  - `dispersed_ohnolog`      524
  - `dispersed_SSD`          450
  - `clustered_cis_ohnolog`  106   <- CLUSTERED group (1)
  - `clustered_SSD_tandem`    73   <- CLUSTERED group (2)

### Environment
- Python: `/mnt/alvand/apps/anaconda2/envs/py3/bin/python3`. Run from the iteration root;
  scripts resolve paths via `Path(__file__).parents[1]`.
- Peak/interval ops: use `pybedtools`/`bedtools` if available, else pure-python interval
  overlap. Keep everything OFFLINE after the one ENCODE download (cache the peak file in
  `inputs/`).

########################################## Goal ##########################################

Test whether **CTCF binding architecture differs between clustered ohnolog loci and clustered
SSD (tandem) loci**. Biological rationale: ohnolog (whole-genome-duplication) clustered pairs
and small-scale-duplication tandem clustered pairs arose by different mechanisms and may sit in
different chromatin/insulator contexts. CTCF marks insulator / TAD-boundary structure, so its
density and distribution around the two kinds of clustered loci is a readable proxy for whether
the two groups occupy distinct genome-architectural niches.

The comparison of interest is strictly the two CLUSTERED categories:
  (1) `clustered_cis_ohnolog`  vs  (2) `clustered_SSD_tandem`.
(Dispersed categories are not clustered loci and are out of scope here.)

########################################## Steps ##########################################

### 1. CTCF peaks from ENCODE (one cell type)
Download a single ENCODE CTCF ChIP-seq peak call, hg38, from an immune-relevant cell type.
**Default choice: GM12878** (EBV-transformed B-lymphoblastoid — matches the lymphoid focus of
the project). Use the **IDR-thresholded / conservative peak set, bed narrowPeak, hg38**.
  - Fetch from the ENCODE portal (`https://www.encodeproject.org`). Record the exact experiment
    accession (ENCSR...) and file accession (ENCFF...), assembly, output-type, and download date
    in `log/` and as a header comment in the script — CONFIRM the accession on the portal at
    download time rather than trusting a hard-coded guess (candidate to verify: GM12878 CTCF
    experiment ENCSR000DZN, hg38 optimal/conservative IDR narrowPeak).
  - Save the raw peak file under `inputs/encode/` (gzip-ok). One cell type is sufficient per the
    request; keep it swappable via a CLI arg so a second cell type can be added later.
  - Sanity-check: number of peaks, chrom set = hg38 primary chroms, and that peaks are point-ish
    narrowPeak intervals; drop non-primary contigs.

### 2. Extend each clustered locus by +/-100 kb
For every cluster in `clustered_paralogs.noKRABHOX.tsv` whose `cluster_id` maps (via
`TF_dup_2x2_classification.tsv`) to `clustered_cis_ohnolog` or `clustered_SSD_tandem`:
  - window = [`cluster_start` - 100000, `cluster_end` + 100000], clipped to chromosome bounds
    (use an hg38 chrom-sizes file; add one to `inputs/` if not present).
  - Emit a BED/TSV of these extended clustered-locus windows tagged with `cluster_id`,
    `dup_class4`, `members`, `span_bp`, and window length.

### 3. Bin each window and quantify CTCF enrichment
  - Bin every extended window into **fixed-width bins (default 5 kb; expose `--binsize`)**.
    Because loci differ in native span, ALSO produce a length-normalized representation: rescale
    each window to a common relative-position axis (e.g. N equal fractional bins across
    [-100 kb ... locus body ... +100 kb], or a metagene-style profile with a fixed flank-bp
    scale + a scaled body) so per-bin CTCF can be averaged/compared across loci of different
    sizes. Document exactly which binning is used for which figure.
  - Per bin, compute CTCF enrichment. Report BOTH:
      (a) peak-overlap density = number of CTCF peaks (or peak-summits) overlapping the bin,
          and/or bp-fraction of the bin covered by peaks;
      (b) if a CTCF signal bigWig is also downloaded, mean signal per bin (optional, secondary).
    Peak-based density from the narrowPeak file is the primary metric (one download, robust).
  - Normalize for comparability: per-window CTCF peak count / window length (peaks per kb), and a
    genome-wide expected-rate baseline (total CTCF peaks / mappable genome length) so enrichment
    can be expressed as observed/expected. Emit per-cluster summaries and per-bin matrices.

### 4. Compare the two clustered groups
  - Per-locus CTCF metrics (peaks/kb over the +/-100 kb window; peaks/kb in the locus body vs the
    flanks; distance of nearest CTCF peak to each member TSS; number of CTCF peaks between the two
    paralog members) grouped by `dup_class4` (only the two clustered classes).
  - Statistics — name the test explicitly in the results text (reproducibility rule below):
      * **Mann-Whitney U** (two-sided) for clustered_cis_ohnolog vs clustered_SSD_tandem on each
        per-locus continuous metric; report U, p, n per group, medians, and effect size.
      * Mean CTCF metagene profiles per group with bootstrap/CI bands; optionally a per-bin
        two-group test (with BH-FDR across bins) to localize where along the window they differ.
  - Figures: (i) mean CTCF metagene/relative-position profile, the two clustered groups overlaid
    with CI bands; (ii) per-locus distribution (box/violin + points) of peaks-per-kb by group;
    (iii) optional heatmap of per-locus binned CTCF, rows = loci, sorted/split by group.

### Outputs (this iteration)
  - `results/clustered_loci_windows.pm100kb.tsv`  (extended windows + class tags)
  - `results/ctcf/ctcf_bins.<binsize>.matrix.tsv`  (per-locus x per-bin CTCF)
  - `results/ctcf/ctcf_per_locus_summary.tsv`      (per-locus metrics + `dup_class4`)
  - `results/ctcf/ctcf_group_comparison.tsv`       (MWU results, medians, effect sizes)
  - `results/figures/ctcf_metagene_clustered_groups.{pdf,png}`
  - `results/figures/ctcf_peaks_per_kb_by_group.{pdf,png}`
  - report in `reports/` (tex -> pdf) summarizing the ENCODE source, method, and the ohnolog-vs-
    SSD CTCF comparison.

########################################## ADD-ON (2026-07-18): ATAC-seq / ENCODE-rE2G in a T cell ##########################################
User request: "do similar analysis for ATAC-seq data from any T cell subtype in ENCODE, especially
the new E2G dataset. The question is enrichment within clustered ohnolog or clustered SSD with
respect to gene pair promoters." This mirrors the CTCF comparison but swaps the assay for chromatin
ACCESSIBILITY (ATAC-seq) and for the ENCODE-rE2G enhancer->gene regulatory-link resource, and
sharpens the readout to the two paralog PROMOTERS.

Cell type = **regulatory CD4 T cell (Treg)** — chosen for the project's autoimmunity/lupus focus.

Two ENCODE data products (hg38 / GRCh38), downloaded 2026-07-18, cached under `inputs/encode/`:
  (A) **ATAC-seq accessibility peaks** — experiment ENCSR159GFS (CD4+CD25+ ab regulatory T cell),
      file **ENCFF519QMA** = IDR-thresholded peaks, bed narrowPeak (57,636 peaks, summit in col10).
      This is the direct accessibility analog to the CTCF narrowPeak and is fed to the existing
      CTCF scripts (03/04/05/06) via `--peaks ... --subdir atac`, so the binning, per-locus,
      flank, and paired-promoter density code is drop-in reused.
  (B) **ENCODE-rE2G predictions** (the "new E2G dataset") — annotation ENCSR759OVP ("enhancer gene
      links in regulatory CD4 T cells", multiomic model using ATAC-seq + H3K27ac + HiC + RNA),
      file **ENCFF393GIF** = thresholded element-gene links, bed3+ (295,377 links; each row = a
      regulatory ELEMENT (chr,start,end) -> TargetGene, with ABC/rE2G Score, DistanceToTSS, and the
      per-element ATAC accessibility column `Open`, plus class = promoter/lncRNA/intergenic/genic).

### Steps
1. ATAC accessibility mirror of the CTCF analysis (scripts 03/04/05/06, `--subdir atac`):
   per-locus ATAC peaks/kb (window / body / flank), metagene, paired-promoter accessibility, flank
   asymmetry; Mann-Whitney U clustered_cis_ohnolog vs clustered_SSD_tandem on each metric.
2. E2G promoter-centric analysis (NEW script `07_e2g_promoters.py`): using ENCFF393GIF,
   - unique E2G elements deduped by (chr,start,end) treated as accessible regulatory elements;
     per-locus element density in the +/-100 kb window/body/flanks and observed/expected vs a
     genome-wide baseline; MWU between the two clustered groups.
   - **Gene-pair promoter readout (the core question):** for each paralog pair, count E2G
     enhancer->gene links whose TargetGene is a pair member; report distinct enhancer elements
     linked per member, links landing on each promoter, and **shared enhancer elements linked to
     BOTH paralogs** (a co-regulation / shared-regulatory-input proxy). Also count `class==promoter`
     E2G elements at the pair TSSs. Compare clustered_cis_ohnolog vs clustered_SSD_tandem (MWU).

### Outputs (add-on)
  - `results/ctcf/atac/*`  and `results/figures/atac/*`   (ATAC mirror; same schema as CTCF, subdir)
  - `results/e2g/e2g_per_locus_summary.tsv`               (per-locus E2G element density + dup_class4)
  - `results/e2g/e2g_pair_promoter_links.tsv`             (per-pair E2G links to member promoters + shared)
  - `results/e2g/e2g_group_comparison.tsv`                (MWU results, medians, effect sizes)
  - `results/figures/e2g/*`                               (E2G group-comparison figures)
  - report in `reports/` (tex -> pdf) covering ATAC + E2G source, method, and the two-group comparison.

### Reproducibility panel (2026-07-18)
The E2G enhancer-sharing result was re-run on 5 additional ATAC-based ENCODE-rE2G sets
(extended model, hg38) via `07_e2g_promoters.py --e2g ... --subdir ... --label ...`, then
aggregated by `08_e2g_reproducibility.py`. ENCODE's ATAC-based rE2G resource offers only 3
T-cell subtypes beyond Treg (memory CD4 = ENCFF440AKC, CD8 = ENCFF993UOE, Jurkat/T-ALL =
ENCFF705FWU); NK (ENCFF359ZSL) and B cell (ENCFF186ABX) were added as lymphoid controls to
complete the panel. Result reproduces in ALL 6 cell types: median shared enhancers per pair =
0 (ohnolog) vs 6-10.5 (SSD-tandem), Mann-Whitney U p from 1.7e-5 to 4.4e-3 (weakest = Jurkat).
Outputs: `results/e2g/<subdir>/`, `results/e2g/e2g_reproducibility_summary.tsv`,
`results/figures/e2g/e2g_reproducibility.{pdf,png}`.

### Add-on caveats
  (v) ATAC/E2G are cell-type specific (Treg here) — unlike constitutive CTCF, accessibility and
      enhancer-gene wiring vary across T-cell subsets; results are conditioned on Treg.
  (vi) ENCODE-rE2G thresholded links are a MODEL prediction (multiomic ABC-style), not a direct
       measurement; "shared enhancer" counts depend on the score threshold baked into the file.

########################################## Reproducibility rule (project-wide) ##########################################
If you run any statistical test, write the summary and test results in the results section,
naming the test and the section explicitly (e.g. "Mann-Whitney U, clustered_cis_ohnolog vs
clustered_SSD_tandem, +/-100 kb window peaks/kb"). Record the ENCODE accession + download date.

### CAVEATS
  (i) One cell type only (per request); CTCF architecture is largely constitutive but not
      identical across cell types — results are conditioned on the chosen line (default GM12878).
  (ii) The two clustered groups are small (n=106 / n=73 TFs -> fewer distinct clusters after
      collapsing members to `cluster_id`); per-locus tests may be underpowered. Report n_clusters
      per group, not n_TFs.
  (iii) Different loci have different native spans; always separate raw fixed-bp binning from
      length-normalized (metagene) binning when comparing groups, and state which is shown.
  (iv) Extended windows can overlap for nearby clusters and can run off chromosome ends — clip to
      hg38 chrom sizes and note any windows shortened by clipping.
