BCBM Contamination Study
Neural-like transcriptional signal in bulk breast cancer brain metastasis is largely explained by adjacent normal brain tissue
A reproducible bioinformatics pipeline that demonstrates a pervasive artefact in bulk transcriptomic studies of breast cancer brain metastasis (BCBM): neuronal and glial genes reported as "tumour-cell neural reprogramming" are, in the majority of cases, markers of admixed normal brain parenchyma rather than tumour-cell biology.

Key Findings
Per-gene contamination verdict — Of 1,254 differentially expressed genes tested in 69 brain-metastasis samples, 120 are contamination-driven (Spearman rho > 0.6 with deconvolution-derived brain-cell content, FDR < 0.05), 973 are tumour-intrinsic candidates (|rho| < 0.3), and 161 remain indeterminate.
C1QL1 reversal — The single most up-regulated gene (log2FC = 5.15) is itself the most completely contamination-derived (rho = 0.72, FDR = 6e-9). Simulation-based correction underestimates its contamination six-fold due to single-nucleus capture bias for secreted synaptic transcripts.
Multi-cohort replication — The gene-level association replicates across three independent BCBM cohorts (pooled rho = 0.918, 95% CI 0.879–0.945, n = 100 brain metastases), computed within brain-metastasis samples only.
Structural vulnerability — The widely used "de novo adaptation" category (genes low in primary, high in metastasis) is structurally enriched for brain-restricted transcripts (rho vs baseline expression r = -0.447, p = 1e-62).
Competence programme — Clean with respect to contamination (only 2 of 1,050 genes contamination-driven) but fails to predict brain metastasis in the expected direction across two independent BMFS cohorts (GSE2603 HR = 0.47/SD opposite to hypothesis; GSE12276 HR = 1.24/SD, p = 0.38).
Literature scan — 10 of 11 (91%) published BCBM bulk transcriptomic studies reporting neural/glial up-regulation report neither tumour purity nor deconvolution.

Directory Structure
D:\BCBM_Project\
├── code\                    # All analysis scripts (Python + R)
│   ├── common.py            # Shared GEO parsers, Ensembl→symbol mapping
│   ├── 00_download.py       # GEO data download
│   ├── 02_paired_de.py      # Paired differential expression (GSE184869)
│   ├── 05_classify.py       # Temporal decomposition (competence vs adaptation)
│   ├── 06_organ_specificity.py  # Brain vs liver/lung/bone metastasis
│   ├── 08-09_survival*.py   # Survival analysis (BMFS / OS)
│   ├── 10_singlecell.R      # Single-nucleus processing (GSE324453)
│   ├── 16_cellchat.R        # CellChat ligand-receptor
│   ├── 19-20_nichenet*.R    # NicheNet analysis
│   ├── 22_contamination.py  # Initial contamination diagnosis
│   ├── 26_bayesprism.R      # BayesPrism deconvolution
│   ├── 29_mixing_simulation.R  # In silico mixing (supporting)
│   ├── 32_dl_176078.R       # Download GSE176078 breast scRNA reference
│   ├── 33_bp_combined.R     # Merged-reference BayesPrism
│   ├── 39_Dfix.py           # Multi-dataset replication (BM-only)
│   ├── 40_C1QL1.py          # C1QL1专项分析
│   ├── 42_module_H.py       # Per-gene contamination verdict (core)
│   ├── 45_module_K.py       # Leave-marker-out negative control
│   ├── 46_J1.R              # Reference profile similarity
│   ├── 47_J4.py             # Neuron+ODG redefinition robustness
│   ├── 48-49_J2J3.py/R      # Class rebalance / astrocyte removal
│   ├── 51_tcga.R            # TCGA-BRCA deconvolution floor
│   ├── 52_batch_baseline.py # Batch effects + baseline-rho relationship
│   ├── 53_make_figures.py   # Main figures 1-6
│   ├── 55_supp_tables.py    # Supplementary tables
│   ├── 56-63_survival*.py   # Survival correction (event verification,
│   │                        #   GSE12276 characteristics, comprehensive)
│   ├── 60_neuron_competence.py  # Neuron-only brain content + competence
│   ├── 64-69_fix*.py        # Figure corrections and manuscript updates
│   └── 68_literature_scan.py    # Literature scan (Table S7)
├── data\raw\                # Downloaded GEO datasets
├── results\
│   ├── tables\              # All intermediate and final result tables
│   └── figures\             # Intermediate figures
├── manuscript\
│   ├── manuscript_v1.md     # Final manuscript (English)
│   ├── figures\             # Publication-quality figures (Fig 1-6, S1-S3)
│   └── supplementary\       # Tables S1-S7
├── report\                  # Chinese analysis report
└── pylibs\                  # Local Python packages (mygene, etc.)


Data Sources
Dataset
Role
Samples
Platform
GSE184869
Discovery: paired primary–BCBM
69 BM, 21 primary, 21 true pairs
RNA-seq (log2 TMM CPM)
GSE324453
Single-nucleus BCBM reference
56,268 nuclei, 12 cell types, 21 patients
snRNA-seq
GSE176078
Breast single-cell reference (deconvolution)
100,064 cells, 26 patients
scRNA-seq
GSE125989
Replication
16 BM
Microarray (GPL571)
GSE14017
Replication
15 BM
Microarray (GPL570)
GSE14018
Replication (excluded, underpowered)
7 BM
Microarray (GPL96)
GSE43837
Excluded (unannotated custom array)
19 BM
Custom (GPL1352)
GSE2603
BMFS survival
82 patients (14 events)
Microarray
GSE12276
BMFS (site-of-relapse endpoint)
196 patients (15 brain events)
Microarray
TCGA-BRCA
Deconvolution floor (descriptive)
50 primaries
STAR counts (recount3)


Analysis Pipeline
The pipeline proceeds in five phases. Scripts are numbered roughly in execution order, though some later scripts (40+) correspond to robustness analyses added during peer review.
Phase 1: Data acquisition and paired differential expression
00_download.py — Download all GEO series-matrix files.
02_paired_de.py — Paired t-tests on 21 true primary–metastasis pairs (FDR < 0.05, |log2FC| > 0.3).
05_classify.py — Temporal decomposition into competence / adaptation / intermediate by baseline expression.
Phase 2: Contamination diagnosis and deconvolution
22_contamination.py — Brain-cell marker scoring; initial evidence of contamination.
26_bayesprism.R — BayesPrism deconvolution with BCBM snRNA reference.
32_dl_176078.R — Download and process GSE176078 breast reference.
33_bp_combined.R — Merged-reference (BCBM + breast) BayesPrism.
51_tcga.R — TCGA-BRCA independent floor analysis.
Phase 3: Core evidence — per-gene contamination verdict
42_module_H.py — Core analysis: per-gene Spearman correlation between expression and brain-cell content within BM samples only; three-class verdict (contamination-driven / tumour-intrinsic / indeterminate).
40_C1QL1.py — C1QL1 专项: highest logFC gene is most contaminated.
39_Dfix.py — Multi-dataset replication within BM samples only; Fisher-z meta-analysis.
Phase 4: Robustness
30_leave_marker_out.py + 45_module_K.py — Stepwise marker removal + expression-matched random-gene negative control.
46_J1.R — Reference cell-type profile similarity.
48_J2.R / 49_J3.py — Class rebalance to 441 / astrocyte removal.
47_J4.py / 60_neuron_competence.py — Redefine brain content as neuron+ODG, then neuron-only; recompute all per-gene rho.
52_batch_baseline.py — Batch effects (within-MAYO, partial correlation), baseline-rho relationship, baseline-matched comparison.
29_mixing_simulation.R / 36_sim_full.R — In silico mixing (supporting evidence only; global fit statistics found uninformative).
Phase 5: Survival, figures, manuscript
56_check_events.py — Verify event indicators (critical correction: original code used event=1).
59_survival_fixed.py — GSE2603 BMFS with correct events; bootstrap CI; DFBETA.
62_gse12276_chars.py — Extract GSE12276 site-of-relapse annotation.
63_survival_comprehensive.py — Two-cohort BMFS; 500 random-gene null with correct coding.
53_make_figures.py — Main figures 1–6.
55_supp_tables.py — Supplementary tables S1–S6.
68_literature_scan.py — Literature scan (Table S7).
66_write_ms.py — Assemble final manuscript.

Dependencies
Python (3.10+)
pandas, numpy, scipy, matplotlib, seaborn
mygene (Ensembl→symbol mapping, installed locally in pylibs\)
statsmodels (optional, for partial correlations)
R (4.6.0)
BayesPrism 2.2.3
Seurat 5.5.1
recount3 1.22.0
Matrix 1.7.5
survival, survminer
CellChat, NicheNet (for earlier exploratory analyses)
Slingshot, SCENIC (for earlier exploratory analyses)
R path
C:\Program Files\R\R-4.6.0\bin\Rscript.exe

Python environment
$env:PYTHONPATH = "D:\BCBM_Project\pylibs"


How to Reproduce
# 1. Set Python path
$env:PYTHONPATH = "D:\BCBM_Project\pylibs"

# 2. Download data
python D:\BCBM_Project\code\00_download.py

# 3. Run R scripts (in order)
Rscript D:\BCBM_Project\code\10_singlecell.R
Rscript D:\BCBM_Project\code\26_bayesprism.R
Rscript D:\BCBM_Project\code\32_dl_176078.R
Rscript D:\BCBM_Project\code\33_bp_combined.R
# ... (see pipeline above for full order)

# 4. Run core Python analyses
python D:\BCBM_Project\code\02_paired_de.py
python D:\BCBM_Project\code\42_module_H.py
python D:\BCBM_Project\code\39_Dfix.py
python D:\BCBM_Project\code\45_module_K.py
# ... (see pipeline above for full order)

# 5. Generate figures and tables
python D:\BCBM_Project\code\53_make_figures.py
python D:\BCBM_Project\code\55_supp_tables.py

Note: GSE184869 is distributed as batch-corrected log2 TMM CPM (no raw counts). BayesPrism assumes a count model, so absolute compositional estimates from this cohort are not quantitatively interpretable. All primary inferences use rank correlation, which is invariant to monotone transformation.

Output Artifacts
Main Figures (manuscript\figures\)
File
Content
Figure1_study_design.png
Study design: contamination vs reprogramming hypotheses
Figure2_volcano.png
Paired DE volcano, coloured by contamination verdict
Figure3_per_gene_verdict.png
log2FC–rho scatter + competence density + overlap bar
Figure4_robustness.png
Batch / brain content by class / neuron-only / marker removal
Figure5_forest.png
Three-cohort replication forest plot (pooled rho = 0.918)
Figure6_baseline.png
Baseline–rho relationship + matched comparison + mixing f

Supplementary Figures
FigureS1_similarity.png — Reference profile similarity heatmap
FigureS2_deconv_comparison.png — Deconvolution across reference configurations
FigureS3_tcga_floor.png — TCGA-BRCA deconvolution floor
Supplementary Tables (manuscript\supplementary\)
File
Content
TableS1_sample_inventory.csv
90 samples by cohort and batch
TableS2_gene_contamination_verdict.csv
Core resource: 1,254 genes, rho/CI/FDR/three-class verdict
TableS3_deconvolution_absolute.csv
Absolute deconvolution estimates (input limitation noted)
TableS4_threshold_sensitivity.csv
Classification threshold sensitivity
TableS5_simulation_global_fit.csv
Simulation global fit (inapplicability noted)
TableS6_survival.csv
GSE2603 + GSE12276 BMFS, DFBETA, null distribution
TableS7_literature_scan.csv
11 studies, purity/deconvolution reporting


Key Methodological Caveats
Input type mismatch — GSE184869 is log-CPM, not raw counts. BayesPrism absolute estimates are descriptive only; rank-based inferences are primary.
Simulation limitation — In silico mixing using snRNA references systematically under-captures cytoplasmic/secreted transcripts (C1QL1 case). Per-gene correlation is the primary method; simulation is supporting only.
Survival sample size — GSE2603 has 14 events (EPV = 7 for two-covariable model, below conventional threshold of 10). The 500-gene-set null is itself shifted below HR = 1, indicating small-sample artefact. GSE12276 provides a second BMFS endpoint (15 events) via site-of-relapse annotation; no OS death-event field exists.
Indeterminate genes — 161 genes cannot be classified at n = 69. Larger cohorts or spatial transcriptomics needed.
Withdrawn results — Earlier analyses using event = 1 (all patients as events) for survival, and Spearman-based f* estimation from mixing simulation, are methodologically invalid and have been withdrawn.

Manuscript
Final manuscript: manuscript\manuscript_v1.md
Target journals: Briefings in Bioinformatics, npj Precision Oncology, Genome Biology (methods/correspondence).

License
Code and processed data available for academic reuse. Raw sequencing data are subject to GEO terms of use.# Contamination-drives-neural-signatures-in-BCBM
Analytical designs that define metastasis-acquired programmes by low expression in the primary tumour are structurally vulnerable to this artefact, and the temporal-decomposition framework as a whole does not yield clinically actionable molecular classes from bulk data.
