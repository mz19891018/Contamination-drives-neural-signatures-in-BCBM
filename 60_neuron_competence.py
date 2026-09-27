# 60_neuron_only_competence.py - J.4 v2 neuron-only + competence overlap
import pandas as pd, numpy as np, openpyxl, re
from scipy.stats import spearmanr, mannwhitneyu
D=r"D:\BCBM_Project"

H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
def kind(s):
    if re.match(r"^BM\d",s): return "BM_named"
    if re.match(r"^\d+M_RCS",s): return "BM_RCS"
    if re.match(r"^MAYO",s): return "BM_MAYO"
    if re.match(r"^BP\d",s): return "P_named"
    if re.match(r"^\d+P_RCS",s): return "P_RCS"
    return "Other"
th["batch"]=[kind(s) for s in th.index]
bm=th[th.batch.str.startswith("BM")]
neuron=bm["Neuron"]
print("=== Neuron-only brain content ===")
print(f"BM neuron: median={neuron.median():.4f}, mean={neuron.mean():.4f}")
prim=th[th.batch.str.startswith("P")]
print(f"Primary neuron: median={prim['Neuron'].median():.4f}, mean={prim['Neuron'].mean():.4f}")
print(f"Astrocyte BM: {bm['Astrocyte'].mean():.4f}, Primary: {prim['Astrocyte'].mean():.4f}")
print(f"ODG BM: {bm['Oligodendrocyte'].mean():.4f}, Primary: {prim['Oligodendrocyte'].mean():.4f}")

# load expression
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("BM")]
bm_names=[header[i] for i in bm_idx]
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv"); emap=dict(zip(pm.ensembl,pm.symbol))
up_genes=pm[(pm.program!="Intermediate")&(pm.mean_logFC_BM_vs_BP>0.3)].symbol.tolist()
expr={}
for r in rows:
    sym=emap.get(r[0])
    if sym in up_genes: expr[sym]=[r[i] for i in bm_idx]
wb.close()
E=pd.DataFrame(expr,index=bm_names).T
common=[s for s in bm_names if s in neuron.index]
E=E[common]; neu=neuron.loc[common]

# recompute rho with neuron-only
res=[]
for g in E.index:
    y=E.loc[g].values.astype(float)
    if np.std(y)<1e-9: continue
    r,p=spearmanr(neu.values,y)
    res.append({"gene":g,"rho_neuron":r})
Rn=pd.DataFrame(res)
M=H.merge(Rn,on="gene",how="inner")
print(f"\n=== J.4 v2: Neuron-only rho vs original (n={len(M)}) ===")
core=["ATP1A2","TUBB4A","C1QL1","GNAO1","DEPTOR","ABCG2","CRYAB","GRB7","SCD","SPHK1","MAG","APLP1"]
for g in core:
    row=M[M.gene==g]
    if len(row): print(f"  {g}: original={row.rho_brain.values[0]:.3f}  neuron-only={row.rho_neuron.values[0]:.3f}")

# classification with neuron-only
def cls(r):
    if r>0.6: return "contamination"
    if abs(r)<0.3: return "tumor_intrinsic"
    return "indeterminate"
M["cls_neuron"]=[cls(r) for r in M.rho_neuron]
print(f"\nNeuron-only classification: {M.cls_neuron.value_counts().to_dict()}")
from sklearn.metrics import cohen_kappa_score
print(f"Cohen's kappa (original vs neuron-only): {cohen_kappa_score(M.verdict.map({'contamination_driven':'contamination','tumor_intrinsic_candidate':'tumor_intrinsic','indeterminate':'indeterminate'}), M.cls_neuron):.3f}")
print(f"Per-gene rho correlation: {spearmanr(M.rho_brain,M.rho_neuron)[0]:.3f}")

# adaptation score vs neuron content
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
adapt_in=[g for g in adapt if g in E.index]
adapt_score=E.loc[adapt_in].mean(0)
print(f"\nAdaptation score vs neuron-only: rho={spearmanr(neu,adapt_score)[0]:.3f}")

# ===== Competence overlap =====
print("\n=== Competence overlap with three classifications ===")
comp=[g.strip() for g in open(rf"{D}\results\tables\geneset_competence.txt",encoding="utf-8") if g.strip()]
comp_in_h=H[H.gene.isin(comp)]
print(f"Total competence genes: {len(comp)}")
print(f"In H verdict table: {len(comp_in_h)}")
print(comp_in_h.verdict.value_counts().to_string())
# list any competence genes that are contamination-driven
contam_comp=comp_in_h[comp_in_h.verdict=="contamination_driven"]
if len(contam_comp):
    print(f"\nCompetence genes classified as contamination-driven ({len(contam_comp)}):")
    print(contam_comp[["gene","rho_brain","logFC_obs"]].to_string())
else:
    print("\nNo competence genes classified as contamination-driven.")
indet_comp=comp_in_h[comp_in_h.verdict=="indeterminate"]
print(f"\nCompetence genes indeterminate ({len(indet_comp)}):")
if len(indet_comp)<=20: print(indet_comp[["gene","rho_brain"]].to_string())
else: print(indet_comp[["gene","rho_brain"]].head(20).to_string())

M.to_csv(rf"{D}\results\tables\J4v2_neuron_only.csv",index=False)
print("\nsaved J4v2_neuron_only.csv")
