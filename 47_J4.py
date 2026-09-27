# 47_J4_no_astrocyte.py - recompute module H with brain_frac = Neuron + ODG only
import pandas as pd, numpy as np, openpyxl, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
th["kind"]=["BrainMet" if re.match(r"^BM|^\d+M_RCS|^MAYO",s) else "Primary" if re.match(r"^BP|^\d+P_RCS",s) else "Other" for s in th.index]
bm=th[th.kind=="BrainMet"]
brain_all=bm[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)
brain_noastro=bm[["Neuron","Oligodendrocyte"]].sum(1)
print("BM n:",len(bm))
print("brain_all mean:",round(brain_all.mean(),3)," brain_noastro mean:",round(brain_noastro.mean(),3))
print("corr(all vs noastro):",round(spearmanr(brain_all,brain_noastro)[0],3))

# load expression
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and re.match(r"^BM|^\d+M_RCS|^MAYO",str(h))]
bm_names=[header[i] for i in bm_idx]
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv"); emap=dict(zip(pm.ensembl,pm.symbol))
obs_lfc=dict(zip(pm.symbol,pm.mean_logFC_BM_vs_BP))
up_genes=pm[(pm.program!="Intermediate")&(pm.mean_logFC_BM_vs_BP>0.3)].symbol.tolist()
expr={}
for r in rows:
    sym=emap.get(r[0])
    if sym in up_genes: expr[sym]=[r[i] for i in bm_idx]
wb.close()
E=pd.DataFrame(expr,index=bm_names).T
common=[s for s in bm_names if s in brain_all.index]
E=E[common]; ba=brain_all.loc[common]; bn=brain_noastro.loc[common]

def calc_rhos(bf):
    res=[]
    for g in E.index:
        y=E.loc[g].values.astype(float)
        if np.std(y)<1e-9: continue
        r,p=spearmanr(bf.values,y)
        res.append({"gene":g,"rho":r})
    return pd.DataFrame(res)

R_all=calc_rhos(ba); R_no=calc_rhos(bn)
M=R_all.merge(R_no,on="gene",suffixes=("_all","_noastro"))
print("\n=== core genes: old vs new rho ===")
core=["ATP1A2","TUBB4A","C1QL1","GNAO1","DEPTOR","ABCG2","CRYAB","GRB7","SCD","SPHK1","MAG","APLP1"]
for g in core:
    row=M[M.gene==g]
    if len(row): print(f"  {g}: all={row.rho_all.values[0]:.3f}  noastro={row.rho_noastro.values[0]:.3f}")

# classification comparison
def cls(r):
    if r>0.6: return "contamination"
    if abs(r)<0.3: return "tumor_intrinsic"
    return "indeterminate"
M["cls_all"]=[cls(r) for r in M.rho_all]
M["cls_no"]=[cls(r) for r in M.rho_noastro]
print("\n=== classification counts ===")
print("all:",M.cls_all.value_counts().to_dict())
print("noastro:",M.cls_no.value_counts().to_dict())

# Cohen's kappa
from sklearn.metrics import cohen_kappa_score
k=cohen_kappa_score(M.cls_all,M.cls_no)
print(f"\nCohen's kappa = {k:.3f}")
print(f"rho correlation (all vs noastro per-gene): {spearmanr(M.rho_all,M.rho_noastro)[0]:.3f}")

M.to_csv(rf"{D}\results\tables\J4_rho_comparison.csv",index=False)
print("saved")
