The goal for this project (in parent folder some work was done but don't use results/files without asking) is to study paralogous transcription factors in the immune system. In particular, I want to know what is unique about paralogous TFs from the same family which are clustered in proximity compared with paralogous TFs which are isolated from singletons. 
Recently Zhu_Nature_2026 (inputs folder) studied paralogous transcription factors in the brain. Some ideas can be borrowed from this study but my questions are in general different.

### Inputs and Outputs
Outputs are organized in the current directory as follows.

- anything used for inputs as `/inputs`
- all codes in a folder called `/scripts`
- all figures / spreadsheets in `/results` (intermediates under `/results/intermediate`)
- any log-related output in `/log`
- all reports written in LaTeX and compiled to PDF in `/reports` (we use `tectonic`)

########################################## Section 1: Curating paralogous TFs ##########################################

In 1-s2.0-S0092867418301065-mmc2 in the inputs folder, 'Table S1. Related to Figure 1B' tab, Lambert et al listed all human TFs.

Part 1- Filter this list and Include only those TFs that says 'Yes' to 'Is TF'. I want TFs from the same family to be referred as paralogous TFs. TFs from the same family are those with the same 'DBD' in the column. Curate a list of paralogous TFs with their DBD family as a column.
   - Done in `scripts/01_curate_paralogous_TFs.py`. The sheet has a two-row header (group label on row 0, e.g. 'Is TF?'; column names on row 1, e.g. ID/Name/DBD). Of 2,765 entries, 1,639 are `Is TF? == Yes`; 1,570 of those carry a usable DBD (69 'Unknown'-DBD entries dropped as unassignable). Defining family = shared DBD and paralogous = family size >= 2 gives 1,563 paralogous TFs across 57 multi-member DBD families (7 singletons).
   - Outputs: `results/paralogous_TFs.tsv` (gene_symbol, ensembl_id, DBD_family, family_size, + step-2 columns), `results/paralogous_TF_family_summary.tsv`, `results/intermediate/lambert_isTF_yes.tsv`.

Part 2- For C2H2 ZF TFs, can you group them more finely, e.g. identify all that have KRAB domain and called them KRAB-ZNF if it makes sense.
   - Done in `scripts/02_subdivide_C2H2_ZF.py`. Lambert has no domain info beyond DBD='C2H2 ZF', so the 747-member family is subdivided by HGNC gene-group membership for the three canonical C2H2-ZF accessory/effector domains: KRAB ('KRAB domain containing'), SCAN ('SCAN domain containing'), BTB/POZ ('BTB domain containing'). HGNC complete set is downloaded fresh to `inputs/hgnc_complete_set.txt` and joined to Lambert by Ensembl ID (745/747 matched; 2 unannotated clone IDs left as 'other'). A single fine label is assigned by priority KRAB > SCAN > BTB; per-gene boolean flags has_KRAB/has_SCAN/has_BTB are also kept (25 genes carry >=2 accessory domains).
   - Result: C2H2 ZF (747) -> KRAB-ZNF 354, C2H2 ZF (other) 316, BTB-ZNF 47, SCAN-ZNF 30. KRAB-ZNF becomes the single largest fine-grained paralog subfamily.
   - Outputs: `results/paralogous_TFs.tsv` updated with `DBD_family_fine`, `family_size_fine`, `has_KRAB/has_SCAN/has_BTB`; `results/C2H2_ZF_subfamily_summary.tsv`; `results/intermediate/c2h2_zf_subfamily_assignment.tsv`.
   - Remove CSD genes (LIN28A which is not a TF) and update paralogous_TFs.tsv and other files.
   - Done in `scripts/03_remove_CSD.py`. The CSD ('cold-shock domain') family in Lambert is LIN28A, LIN28B, YBX1, YBX2, YBX3 — RNA / single-stranded-nucleic-acid binders, not bona fide sequence-specific DNA-binding TFs. The whole family is dropped (idempotent: removes any row whose DBD_family is in EXCLUDE_FAMILIES={CSD}). Removing an entire family leaves all other families' sizes unchanged, so family_size/family_size_fine stay correct.
   - Result: paralogous TFs 1,563 -> 1,558; coarse DBD families 57 -> 56; fine-grouping family rows 60 -> 59.
   - Outputs: `results/paralogous_TFs.tsv` and `results/intermediate/lambert_isTF_yes.tsv` rewritten without CSD; `results/paralogous_TF_family_summary.tsv` regenerated; `results/intermediate/removed_genes.tsv` audits the 5 dropped genes.

Part 3- Curating clustered paralogs from singleton paralogs 

For human genome build hg38, find all paralogous-gene clusters on the same chromosome with the following criteria. Write them in file reflecting clustered paralogs. Everything else goes to singleton paralogs. 

1) Genes to be tested must come from the same family paralogous_TF_family_summary.tsv. 

2) Clustered TFs are defined as TFs from the same family and on the same chromosome and located in genomic proximity in clusters of 2 or more. Adjacency between same-family members A and B is established iff no *other* protein-coding gene has its gene body strictly contained in the inter-gene gap `(A.end, B.start)`. Other same-family members in the gap, and any genes that overlap A or B rather than lie strictly between them, do **not** block adjacency. Clusters are the connected components of the adjacency graph (per family). 

3) Non-coding genes between paralog members are tolerated. lncRNAs whose body intersects any inter-member gap are captured in the `lncRNAs_between` column.

4) Promoters are defined as the TSS (strand-aware gene start) drawn from the GENCODE v45 basic annotation, hg38 primary assembly. The catalog records each member's TSS in `member_tss`.

5) add a column showing member_strands e.g. +,-,-

   - Done in `scripts/04_cluster_paralogs.py`. Coordinates/TSS from GENCODE v45 basic (downloaded fresh to `inputs/gencode.v45.basic.annotation.gtf.gz`, main chroms only: 20,036 protein_coding, 19,370 lncRNA). Family = `DBD_family_fine`. Members matched to GENCODE by versionless Ensembl ID (1,555) then symbol (1); 2 unmapped (DUX1/DUX3, absent from GENCODE main assembly) -> `intermediate/unmapped_members.tsv`. Adjacency tested on consecutive same-family members by genomic position (sufficient: a blocker strictly between consecutive members also lies strictly between any wider spanning pair); the blocker set is all protein-coding gene bodies except the family's own members. lncRNAs overlapping any inter-member gap are listed in `lncRNAs_between`.
   - Result: 134 clusters (>=2 members) covering 507 members vs 1,049 singleton members (507+1,049 = 1,556 mapped). 88 clusters have >=1 lncRNA between members. Top clustering families: KRAB-ZNF (52), C2H2 ZF other (23), Homeodomain (21); largest single cluster = 22-member chr19 KRAB-ZNF array. ETS pairs recovered: ETS1-FLI1 (chr11, head-to-head), ERG-ETS2 (chr21), ELF5-EHF (chr11), ETV3L-ETV3 (chr1); SPI1/PU.1 and ELF4 fall out as singletons.
   - Outputs: `results/clustered_paralogs.tsv` (cluster_id, family, parent_DBD_family, chrom, n_members, members, member_strands, member_tss, member_starts/ends, cluster_start/end, span_bp, n_lncRNAs_between, lncRNAs_between); `results/singleton_paralogs.tsv` (gene_symbol, ensembl_id, DBD_family, family, chrom, start, end, strand, tss); `results/intermediate/gencode_v45_genes.tsv`, `results/intermediate/unmapped_members.tsv`.

IMPORTANT: Major outputs used for the downstream analyses are clustered_paralogs.tsv and singleton_paralogs.tsv

################################################## Section 2: Curating species for gene expression analyses ###############################################################

One goal of this project is to evaluate the expression of paralogous TFs in immune cells using available bulk or single-cell expression data across species. Can you search beyond human and mice, in which other species single cell or bulk RNA data exist for immune cells like blood, spleen or lymph nodes. Just search and report. 
