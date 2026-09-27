# 55_supp_tables.py - generate supplementary tables
import pandas as pd, numpy as np, re, os
D=r"D:\BCBM_Project"; ST=rf"{D}\results\tables\supplementary"; os.makedirs(ST,exist_ok=True)

# S1: Sample inventory
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
def kind(s):
    if re.match(r"^BM\d",s): return "BrainMet_named"
    if re.match(r"^\d+M_RCS",s): return "BrainMet_RCS"
    if re.match(r"^MAYO",s): return "BrainMet_MAYO"
    if re.match(r"^BP\d",s): return "Primary_named"
    if re.match(r"^\d+P_RCS",s): return "Primary_RCS"
    return "Other"
th["group"]=[kind(s) for s in th.index]
th["brain_total"]=th[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)
s1=th.reset_index().rename(columns={"index":"sample_id"})[["sample_id","group","Tumor","Neuron","Astrocyte","Oligodendrocyte","brain_total"]]
s1.to_csv(rf"{ST}\TableS1_sample_inventory.csv",index=False)
print("S1:",len(s1),"samples")
print(s1.group.value_counts().to_string())

# S4: Threshold sensitivity
H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
def classify(rho,hi,lo):
    if rho>hi: return "contamination_driven"
    if abs(rho)<lo: return "tumor_intrinsic_candidate"
    return "indeterminate"
rows=[]
for hi,lo in [(0.7,0.2),(0.6,0.3),(0.5,0.4)]:
    cls=H.rho_brain.apply(lambda r: classify(r,hi,lo))
    vc=cls.value_counts()
    rows.append({"threshold":f"rho>{hi}/|rho|<{lo}",
        "contamination_driven":vc.get("contamination_driven",0),
        "tumor_intrinsic":vc.get("tumor_intrinsic_candidate",0),
        "indeterminate":vc.get("indeterminate",0)})
pd.DataFrame(rows).to_csv(rf"{ST}\TableS4_threshold_sensitivity.csv",index=False)
print("\nS4:")
print(pd.DataFrame(rows).to_string())

# S5: simulation global fit stats (copy Anew)
if os.path.exists(rf"{D}\results\tables\Anew_metrics.csv"):
    a=pd.read_csv(rf"{D}\results\tables\Anew_metrics.csv")
    a.to_csv(rf"{ST}\TableS5_simulation_global_fit.csv",index=False)
    print("\nS5 copied:",len(a),"rows")

# S6: survival combined
s6=[]
for f,name in [(rf"{D}\results\tables\survival_gse2603.csv","GSE2603_BMFS"),(rf"{D}\results\tables\survival_gse12276.csv","GSE12276_OS")]:
    if os.path.exists(f):
        df=pd.read_csv(f); df["cohort"]=name; s6.append(df)
if s6:
    pd.concat(s6).to_csv(rf"{ST}\TableS6_survival.csv",index=False)
    print("S6 saved")

# S7: literature scan methodology (search was restricted; provide framework)
s7=pd.DataFrame([
    {"item":"Search strategy","detail":"PubMed: ('breast cancer' AND 'brain metastasis' AND ('neural' OR 'neuronal' OR 'glial' OR 'neuroendocrine' OR 'neural mimicry')) AND 'transcriptome' OR 'microarray' OR 'RNA-seq'"},
    {"item":"Inclusion criteria","detail":"Bulk transcriptomic studies of human BCBM reporting up-regulation of neural/glial genes"},
    {"item":"Exclusion criteria","detail":"Single-cell only, cell-line only, review articles, non-English"},
    {"item":"Extraction fields","detail":"First author, year, platform, N (BM samples), whether tumour purity reported, whether deconvolution performed, whether neural genes nominated as targets"},
    {"item":"Note","detail":"Systematic quantification to be completed at submission; this table documents the protocol. The central claim does not depend on a specific count."},
])
s7.to_csv(rf"{ST}\TableS7_literature_scan_protocol.csv",index=False)
print("S7 protocol saved")

# Copy S2 (verdict) and S3 (fractions)
import shutil
shutil.copy(rf"{D}\results\tables\gene_contamination_verdict.csv", rf"{ST}\TableS2_gene_contamination_verdict.csv")
shutil.copy(rf"{D}\results\tables\bp_combined_fractions.csv", rf"{ST}\TableS3_deconvolution_absolute.csv")
print("S2, S3 copied")
print("\nALL SUPPLEMENTARY TABLES DONE")
