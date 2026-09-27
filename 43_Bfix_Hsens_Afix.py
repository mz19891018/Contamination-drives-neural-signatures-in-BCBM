# 43_Bfix_Hsens_Afix.py
import pandas as pd, numpy as np
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
import re
th["kind"]=["BrainMet" if re.match(r"^BM|^\d+M_RCS|^MAYO",s) else "Primary" if re.match(r"^BP|^\d+P_RCS",s) else "Other" for s in th.index]
bm=th[th.kind=="BrainMet"]
brain=bm[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)

# need per-gene expression in BM samples - reload from module H output? We saved verdict but not expression.
# Recompute expression quickly from xlsx for adaptation genes
import openpyxl
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and re.match(r"^BM|^\d+M_RCS|^MAYO",str(h))]
bm_names=[header[i] for i in bm_idx]
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv"); emap=dict(zip(pm.ensembl,pm.symbol))
expr={}
for r in rows:
    sym=emap.get(r[0])
    if sym in adapt: expr[sym]=[r[i] for i in bm_idx]
wb.close()
E=pd.DataFrame(expr,index=bm_names).T
E=E[[s for s in bm_names if s in brain.index]]; bf=brain.loc[E.columns]
print("expr:",E.shape," BM:",len(bf))

def boot_rho(x,y,n=1000):
    rng=np.random.default_rng(0); rs=[]
    for _ in range(n):
        idx=rng.choice(len(x),len(x),replace=True); rs.append(spearmanr(x[idx],y[idx])[0])
    return np.percentile(rs,2.5),np.percentile(rs,97.5)

# B-fix: progressive removal of top brain markers (data-driven = module H contamination-driven, sorted by rho)
markers=H[H.verdict=="contamination_driven"].sort_values("rho_brain",ascending=False).gene.tolist()
adapt_in=[g for g in adapt if g in E.index]
print("adapt genes with expr:",len(adapt_in)," markers in adapt:",len(set(markers)&set(adapt_in)))
print("\n=== B-fix: progressive marker removal ===")
for nrem in [0,25,50,100,150]:
    rem=set(markers[:nrem])
    keep=[g for g in adapt_in if g not in rem]
    if len(keep)<10: print(f"  remove {nrem}: only {len(keep)} left, insufficient"); continue
    score=E.loc[keep].mean(0)
    r,p=spearmanr(bf.values,score.values); lo,hi=boot_rho(bf.values,score.values)
    print(f"  remove top-{nrem} markers: n_remain={len(keep)} rho={r:.3f} CI[{lo:.3f},{hi:.3f}]")

# H sensitivity: thresholds
print("\n=== H sensitivity: threshold variation ===")
for rhi,rlo in [(0.7,0.2),(0.6,0.3),(0.5,0.4)]:
    cd=((H.rho_brain>rhi)&(H.rho_FDR<0.05)).sum()
    ti=(H.rho_brain.abs()<rlo).sum()
    ind=len(H)-cd-ti
    print(f"  rho>{rhi} / |rho|<{rlo}: contamination={cd}, tumor_intrinsic={ti}, indeterminate={ind}")

# A-fix: RMSE/Deming on brain-specific subset (contamination-driven genes)
sim=pd.read_csv(rf"{D}\results\tables\sim_full_logFC.csv",index_col=0)
obs=dict(zip(pm.symbol,pm.mean_logFC_BM_vs_BP))
brain_sub=[g for g in markers if g in sim.index and g in obs]
print(f"\n=== A-fix: RMSE/Deming on {len(brain_sub)} brain-specific genes ===")
fs=[float(c) for c in sim.columns]
o=np.array([obs[g] for g in brain_sub]); S=sim.loc[brain_sub].values
rmse=np.sqrt(np.mean((S-o[:,None])**2,axis=0))
for f,r in zip(fs,rmse): print(f"  f={f:.3f} RMSE={r:.3f}")
print(f"  f* RMSE min = {fs[np.argmin(rmse)]}")
