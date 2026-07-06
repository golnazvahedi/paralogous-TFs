% The Buffering Code
% Genetic redundancy as the hidden variable that decides whether a mutation causes disease — and how to retune it
% NIH Director's Pioneer Award — concept proposal

# The one-sentence vision

I propose that the fundamental unit of immune-disease genetics is not the gene but the **dynamically buffered paralog network**; that 500 million years of genome-duplication history wrote a latent **"buffering code"** that determines *whether, where, and in which clinical direction* perturbing a gene causes disease; and that learning to **read and rewrite this code** will let us predict a variant's clinical consequence from sequence alone and treat disease by **retuning redundancy** rather than fixing the broken gene.

# The paradigm it overturns

Human genetics is gene-centric: we ask "is this gene mutated?" and "is the disease Mendelian or complex?" as if those were fixed properties. The data motivating this proposal say otherwise. Across 1,153 paralogous transcription factors, **whether a gene causes monogenic versus polygenic immune disease is set not by the gene but by its paralog's buffering state**:

- Genes whose duplicate is **unbuffered** (dispersed whole-genome-duplication "ohnologs" — *SPI1*/PU.1, *ELF4*, *GATA2*, *IRF8*) are enriched for **monogenic** inborn errors of immunity (odds ratio ≈ 3.5; ohnolog-vs-small-scale-duplicate OR ≈ 4.3, p ≈ 10⁻⁴).
- The *same evolutionary class of gene*, when its duplicate is **buffered** by a co-regulated neighbor (clustered *cis*-ohnologs — *ETS1*/*FLI1*/*ERG*, the *STAT* cluster), flips to being the **single most polygenic** category (top GWAS enrichment, OR ≈ 1.6) while carrying **no Mendelian disease at all**.

That flip — the same class of gene landing in two entirely different diagnostic categories purely as a function of redundancy — is the spark. It implies a hidden variable, *buffering*, that no current framework measures, and that reorganizes disease genetics from "genes" to "buffered dyads."

# The central hypothesis

**Genetic buffering is a quantifiable, cell-type-specific, time-varying property that governs the genotype–phenotype map.** Formally, every paralog dyad has a **buffering coefficient *B*(cell type, state, time)** — the degree to which a partner absorbs perturbation of its sibling. I hypothesize that *B* is:

1. **Written by evolution.** It is predictable *a priori* from deep features already shown to matter — 2R-ohnolog provenance, *cis*-clustering, gene-origin age (ohnologs are pre-first-WGD ancient), and even the alpha/beta allopolyploid subgenome of origin. Evolution is a completed, genome-wide perturbation screen we can simply read.
2. **Dynamic and local.** *B* is not a constant. The same gene is buffered in a cell type/state where its paralog is co-expressed and **unbuffered where the partner falls silent** — so disease erupts specifically in the buffering-weakest window. This single idea could explain incomplete penetrance, tissue specificity, age-of-onset, and the episodic nature of autoimmune flares as **buffering failures in space and time**.
3. **The dial that sets disease direction.** I hypothesize *B* determines not just *whether* but *which* disease — where a perturbation lands on the immunodeficiency ↔ autoinflammation ↔ autoimmunity axis — because the framework already separates identity-determining (dispersed) from cytokine-effector (clustered) ohnologs.

# Three falsifiable, field-redefining predictions

- **The architecture prediction.** Given only a gene's evolutionary buffering features, I can predict whether its disease variants will be Mendelian or polygenic — *before any patient is sequenced.*
- **The cell-type prediction.** A buffered disease gene will become Mendelian-like in the specific immune cell type/state where its paralog is naturally extinguished — i.e., we can predict the *cellular site* of pathology from a buffering map.
- **The reversibility prediction.** Disease driven by loss of buffering can be corrected by **retuning the partner paralog**, not the mutated gene — a fundamentally new therapeutic target.

# The transformative program

**1. MEASURE — the first genome-wide, single-cell Buffering Atlas of the immune system.** Move beyond correlation to causation: titratable, paired perturbation (degron dose-series plus single- versus double-paralog knockdown) read by single-cell transcriptomics across primary human immune cell types and activation states, to assign an empirical *B* to every paralog dyad in every cellular context. This turns "buffering" from a concept into a measured coordinate of the genome.

**2. PREDICT — a generative "evolution-to-clinic" model.** Train a model that maps deep genomic history (duplication mode, age, subgenome, *cis*-architecture, paralog co-expression) → buffering → predicted disease architecture, direction, and cellular site. The audacious end state: a *paleogenomic disease oracle* that classifies a never-before-seen variant's clinical behavior from sequence and evolutionary context alone.

**3. REWIRE — paralog-rebalancing as a new therapeutic modality.** Treat the redundant paralog as the genome's built-in shock absorber and learn to dial it (CRISPR activation/interference, targeted protein degraders, or delivery of a "synthetic buffer" paralog) to restore network output in a buffering-failure disease. This is **not science fiction** — oncology has already proven paralogs are real, drug-like dependencies (*SMARCA2/4*, *ARID1A/B*, *MAGOH/MAGOHB* synthetic lethality). The pioneering leap is inverting that logic: in immunity, *restore* buffering to quiet disease rather than break it to kill tumors.

# Why this is Pioneer-appropriate (high-risk / high-reward)

- **Why it won't happen through normal grants.** The prevailing grant culture rewards one-gene–one-disease mechanism. A theory whose unit is the *dynamically buffered network*, and whose payoff is a predictive law, cuts across that grain; it needs a single investigator empowered to chase a unifying idea rather than a hypothesis-confirming aim.
- **The risk is real and stated honestly.** *B* may be too context-dependent to generalize; the evolution-to-clinic prior may not transfer beyond the gene families where the signal is cleanest; paralog-rebalancing may face the usual delivery and specificity hurdles. Any of these could falsify the strong form of the theory.
- **The reward is a paradigm.** If buffering is even partly the hidden variable, we gain (i) a mechanistic explanation for penetrance and for the immunodeficiency/autoinflammation/autoimmunity spectrum; (ii) variant interpretation that predicts clinical architecture from sequence; and (iii) an entirely new therapeutic axis. The principle is **general** — the same duplication-encoded buffering should shape robustness and disease in the brain (where 2R ohnologs were recently shown to drive cell-type evolution), in cancer, and across vertebrate physiology. Immunity is the ideal proving ground because it has both the cleanest genetics (curated monogenic and polygenic catalogs) and the cell-resolution tools.

# Why now

The convergence is new: complete paralog/ohnolog and alpha/beta subgenome maps, curated monogenic (inborn errors of immunity) and polygenic (GWAS) disease catalogs, single-cell perturbation at scale, and — from the work motivating this proposal — the first evidence that one evolutionary axis (buffering) reorganizes the entire monogenic-versus-polygenic landscape. The framework, the catalogs, and the integrating analyses already exist; the Pioneer Award would fund the leap from *observing* the buffering code to *measuring, predicting, and rewriting* it.

# What success would change

A successful program replaces "which gene is broken?" with "which buffered network failed, in which cell, and how do we re-balance it?" — delivering (1) a single-cell Buffering Atlas of human immunity, (2) a sequence-to-clinic predictor of disease architecture and direction, and (3) first-in-class paralog-rebalancing therapeutic concepts. Each is independently valuable; together they would establish genetic buffering as a measurable, predictable, and *druggable* layer of the genome.
