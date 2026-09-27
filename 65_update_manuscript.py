# 65_update_manuscript.py - update manuscript with all corrections
import re

f=r"D:\BCBM_Project\manuscript\manuscript_v1.md"
with open(f,encoding="utf-8") as fh: c=fh.read()

# 1. Replace 3.10 body (from after title to before 3.11)
old_310_start="Because the competence programme is the component unaffected by contamination, we asked whether it carries clinical information."
old_310_end="The GSE2603 BMFS association is hypothesis-generating only."
idx_start=c.find(old_310_start)
idx_end=c.find(old_310_end)+len(old_310_end)

new_310="""Because the competence programme is the component unaffected by contamination, we asked whether it carries clinical predictive value. We analysed two independent BMFS cohorts.

**GSE2603.** Using the correct brain-metastasis event indicator (`bm event` field, 0/1), excluding 22 MDA-MB-231 cell-line samples and specimens with missing follow-up, 82 evaluable patients remained (14 brain-metastasis events, 68 censored). The competence score was associated with reduced risk in this cohort: HR = 0.47 per standard deviation (p = 0.005; bootstrap 95% CI 0.26–0.80 from 1,000 resamples), and HR = 0.42 after adjustment for proliferation (p = 0.006). This direction is opposite to the programme's defining hypothesis — competence genes are selected for high expression in primary tumours that later metastasise to the brain, so higher scores should confer higher, not lower, risk.

We investigated whether this reversal reflects technical confounding. The competence score was higher in ER-negative than ER-positive primaries (mean z-score 0.28 vs −0.22), but ER-negative disease is itself a brain-metastasis risk factor; if the score were merely proxying ER status, the HR should be positive, not negative. The score contains GRB7 (HER2 amplicon), and HER2-positive disease also confers higher risk — again working against, not toward, the observed negative association. DFBETA influence analysis identified three high-leverage points (GSM50075, a censored patient with score −2.26; GSM50112 and GSM50066, early events at 1.2 years with scores +0.91 and −2.56), but leave-one-out estimates remained below 1.0.

Crucially, a null distribution of 500 size-matched random gene sets, run with the correct event coding, was itself shifted below unity (median null HR = 0.83, 99% of null sets had HR < 1). In this small-event dataset, gene-set scores of this size tend toward an apparent protective association regardless of gene identity. The observed HR of 0.47 is more extreme than the null (empirical p < 0.001), but the null's deviation from 1.0 indicates that the protective direction is a dataset-level property, not a specific property of the competence gene set.

**GSE12276.** We constructed a second BMFS endpoint from the "site of relapse (brain or other)" field (brain or brain+other = event; other/local or no relapse = censored), yielding 196 evaluable patients with 15 brain events. Here the competence score showed the expected positive direction but was not significant: HR = 1.24 per SD (p = 0.38). The two cohorts thus give opposite directions — significantly negative in GSE2603, nonsignificantly positive in GSE12276 — and neither supports the hypothesis that higher primary-tumour competence expression predicts brain metastasis. A previous overall-survival analysis of GSE12276 that assumed event = 1 for all patients is methodologically invalid and withdrawn; no death-event indicator exists in the GEO metadata.

We conclude that the competence programme does not provide a validated basis for CNS-surveillance stratification. Its failure to predict in the expected direction, alongside the contamination-driven failure of the adaptation programme, indicates that the temporal-decomposition framework as a whole does not yield clinically actionable molecular classes from bulk transcriptomic data."""

c=c[:idx_start]+new_310+c[idx_end:]

# 2. Update 3.8 to add clinical prediction failure note
old_38_end="their low fold changes (< 0.5) are consistent with this interpretation."
new_38_end=old_38_end+" However, this contamination-independence does not translate to clinical utility: as shown in Section 3.10, the competence score fails to predict brain metastasis in the expected direction across two independent BMFS cohorts. The programme is clean with respect to contamination but not valid as a predictive biomarker."
c=c.replace(old_38_end,new_38_end)

# 3. Update Methods 2.9
old_29="""Cox models used scores standardised per standard deviation. For GSE2603 (brain-metastasis-free survival), the event indicator was extracted from the `bm event` field (0/1) in the GEO series matrix; 22 MDA-MB-231 cell-line samples and specimens with missing follow-up or event data were excluded, leaving 82 evaluable patients (14 events). Models were adjusted for a proliferation score (average of MKI67, TOP2A, PCNA, CCNB1). For GSE12276, the GEO metadata contains survival time but no death-event indicator; the overall-survival analysis is therefore not reported, and a previous analysis that assumed event = 1 for all samples is withdrawn. Brain-metastasis-free survival and overall survival were analysed separately and not pooled."""
new_29="""Cox models used scores standardised per standard deviation. For GSE2603 (brain-metastasis-free survival), the event indicator was extracted from the `bm event` field (0/1) in the GEO series matrix; 22 MDA-MB-231 cell-line samples and specimens with missing follow-up or event data were excluded, leaving 82 evaluable patients (14 events). For GSE12276, a second BMFS endpoint was constructed from the "site of relapse (brain or other)" field (brain or brain+other = event; other/local or no relapse = censored), yielding 196 patients (15 events). Models were adjusted for a proliferation score (average of MKI67, TOP2A, PCNA, CCNB1). Bootstrap confidence intervals (1,000 resamples) and DFBETA influence statistics were computed for the GSE2603 model. A null distribution of 500 size-matched random gene sets was generated with the correct event coding. No death-event indicator exists in the GSE12276 GEO metadata; a previous overall-survival analysis that assumed event = 1 is withdrawn. Brain-metastasis-free survival endpoints were analysed separately by cohort and not pooled."""
c=c.replace(old_29,new_29)

# 4. Update dataset table GSE12276
c=c.replace("| GSE12276 | Overall survival (not analysable) | n = 196; no death-event indicator in GEO metadata |",
            "| GSE12276 | BMFS (site-of-relapse endpoint) | n = 196; 15 brain events; no death-event field for OS |")

# 5. Update Figure 4 legend (b) mean label
c=c.replace("(b) Mean brain-cell fraction by class, brain metastasis versus primary; only neuron fraction clearly distinguishes groups (1.5×), with fold-change ratios annotated.",
            "(b) Mean brain-cell fraction by class, BM versus primary; neuron 1.8×, astrocyte 1.2×, oligodendrocyte 1.1× (mean ratios).")

# 6. Update Figure 4 legend bottom section
c=c.replace("(b) Mean brain-cell fraction by class, BM versus primary; only neuron fraction distinguishes groups.",
            "(b) Mean brain-cell fraction by class, BM versus primary (mean ratios annotated).")

# 7. Update Table S6 description
c=c.replace("**Table S6.** Survival analysis: GSE2603 BMFS with correct event indicators (82 patients, 14 events; competence HR = 0.47/SD, p = 0.005; adjusted HR = 0.42, p = 0.006). GSE12276 OS not analysable (no death-event field in GEO metadata).",
            "**Table S6.** Survival analysis: GSE2603 BMFS (82 patients, 14 events; competence HR = 0.47/SD, p = 0.005; bootstrap CI [0.26, 0.80]; 500-gene-set null median HR = 0.83); GSE12276 BMFS from site-of-relapse (196 patients, 15 events; HR = 1.24, p = 0.38); DFBETA influence table.")

# 8. Update limitations - survival
old_lim="""Sixth, survival analysis is limited by data availability. GSE2603 BMFS has a small number of events (n = 14), making the protective association of the competence score hypothesis-generating. GSE12276 lacks a death-event indicator in its GEO metadata, so overall survival could not be properly analysed; a previous event = 1 assumption was invalid and has been withdrawn. Larger cohorts with complete event data are needed to assess the clinical utility of the competence score."""
new_lim="""Sixth, survival analysis is limited by data availability. GSE2603 BMFS has only 14 events (events-per-variable = 7 for the two-covariable model, below the conventional threshold of 10), and the null distribution of random gene sets is itself shifted below HR = 1, indicating that the apparent protective direction is partly a small-sample artefact. GSE12276 provides a second BMFS endpoint via site-of-relapse annotation (15 events) but lacks a death-event indicator for OS; a previous event = 1 assumption was invalid and has been withdrawn. The two cohorts give opposite directions for the competence score, precluding any clinical claim. Larger cohorts with complete event data are needed."""
c=c.replace(old_lim,new_lim)

with open(f,"w",encoding="utf-8") as fh: fh.write(c)
print("Manuscript updated successfully")
print(f"File length: {len(c)} chars")
