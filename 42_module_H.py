# 42_module_H.py - per-gene contamination correlation (core v3.0)
import pandas as pd, numpy as np, openpyxl, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

# ---- load bulk expression ----
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]
rows=ws.iter_rows(values_only=True)
header=list(next(rows))
# BM samples: BM* or *M_RCS
bm_idx=[i for i,h in enumerate(header) if h and (re.match(r"^BM",str(h)) or re.match(r"^\d+M_RCS$",str(h)) or re.match(r"^MAYO",str(h)))]
bp_idx=[i for i,h in enumerate(header) if h and (re.match(r"^BP",str(h)) or re.match(r"^\d+P_RCS$",str(h)))]
bm_names=[header[i] for i in bm_idx]
print("BM samples:",len(bm_names)," BP:",len(bp_idx))

# ensembl->symbol
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
emap=dict(zip(pm.ensembl,pm.symbol))
obs_lfc=dict(zip(pm.symbol,pm.mean_logFC_BM_vs_BP))
up_genes=pm[(pm.program!="Intermediate")&(pm.mean_logFC_BM_vs_BP>0.3)].symbol.tolist()
print("up genes:",len(up_genes))

# read expression for up genes only
expr={}
for r in rows:
    eid=r[0]
    sym=emap.get(eid)
    if sym in up_genes:
        expr[sym]=[r[i] for i in bm_idx]
wb.close()
E=pd.DataFrame(expr,index=bm_names).T  # genes x BM samples
print("expr matrix:",E.shape)

# ---- BayesPrism brain fractions (combined ref) ----
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
brain=th[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)
# match BM samples
common=[s for s in bm_names if s in brain.index]
print("BM samples with BP fraction:",len(common))
bf=brain.loc[common]
E=E[common]

# ---- Module H: per-gene rho ----
def boot_ci(x,y,n=2000):
    rng=np.random.default_rng(1); rs=[]
    for _ in range(n):
        idx=rng.choice(len(x),len(x),replace=True)
        rs.append(spearmanr(x[idx],y[idx])[0])
    return np.percentile(rs,2.5),np.percentile(rs,97.5)

res=[]
for g in E.index:
    x=bf.values; y=E.loc[g].values.astype(float)
    if np.std(y)<1e-9: continue
    r,p=spearmanr(x,y)
    lo,hi=boot_ci(x,y)
    res.append({"gene":g,"logFC_obs":obs_lfc.get(g,np.nan),
                "rho_brain":round(r,3),"rho_CI_lo":round(lo,3),"rho_CI_hi":round(hi,3),"rho_p":p})
H=pd.DataFrame(res)
H["rho_FDR"]=np.minimum(H.rho_p*len(H),1.0)  # BH approx
def verdict(r,fdr):
    if r>0.6 and fdr<0.05: return "contamination_driven"
    if abs(r)<0.3: return "tumor_intrinsic_candidate"
    return "indeterminate"
H["verdict"]=[verdict(r,f) for r,f in zip(H.rho_brain,H.rho_FDR)]
H.to_csv(rf"{D}\results\tables\gene_contamination_verdict.csv",index=False)
print("\n=== verdict counts (n=%d BM) ==="%len(common))
print(H.verdict.value_counts())
print("\n=== top contamination-driven (by rho) ===")
print(H[H.verdict=="contamination_driven"].sort_values("rho_brain",ascending=False).head(15)[["gene","logFC_obs","rho_brain","rho_FDR"]].to_string(index=False))
print("\n=== tumor_intrinsic_candidate (top by logFC) ===")
print(H[H.verdict=="tumor_intrinsic_candidate"].sort_values("logFC_obs",ascending=False).head(15)[["gene","logFC_obs","rho_brain"]].to_string(index=False))

# ---- G-new: C1QL1 ----
c1=H[H.gene=="C1QL1"]
print("\n=== G-new C1QL1 ===")
print(c1.to_string(index=False))

# ---- F-new step1: primary brain fraction decomposition ----
th2=th.copy()
th2["kind"]=["BrainMet" if re.match(r"^BM|^\d+M_RCS|^MAYO",s) else "Primary" if re.match(r"^BP|^\d+P_RCS",s) else "Other" for s in th2.index]
prim=th2[th2.kind=="Primary"]
print("\n=== F-new step1: Primary brain fraction by cell type (n=%d) ==="%len(prim))
print(prim[["Neuron","Astrocyte","Oligodendrocyte"]].mean().round(3))
print("total primary brain:",round(prim[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1).mean(),3))
bm_all=th2[th2.kind=="BrainMet"]
print("total BM brain:",round(bm_all[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1).mean(),3))
