# 54_regen_K_LMO.py - regenerate module K null and stepwise LMO
import pandas as pd, numpy as np, openpyxl, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
ti=H[H.verdict=="tumor_intrinsic_candidate"].gene.tolist()
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
def kind(s):
    if re.match(r"^BM\d",s): return "BM_named"
    if re.match(r"^\d+M_RCS",s): return "BM_RCS"
    if re.match(r"^MAYO",s): return "BM_MAYO"
    if re.match(r"^BP\d",s): return "P_named"
    if re.match(r"^\d+P_RCS",s): return "P_RCS"
    return "Other"
th["kind"]=[kind(s) for s in th.index]
bm=th[th.kind.str.startswith("BM")]
brain=bm[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)

wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("BM")]
bm_names=[header[i] for i in bm_idx]
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv"); emap=dict(zip(pm.ensembl,pm.symbol))
expr={}
for r in rows:
    sym=emap.get(r[0])
    if sym in set(ti)|set(adapt): expr[sym]=[r[i] for i in bm_idx]
wb.close()
E=pd.DataFrame(expr,index=bm_names).T
common=[s for s in bm_names if s in brain.index]
E=E[common]; bf=brain.loc[E.columns]

# Stepwise LMO
markers=H[H.verdict=="contamination_driven"].sort_values("rho_brain",ascending=False).gene.tolist()
lmo_rows=[]
for n in [0,25,50,100,150]:
    rem=set(markers[:n])
    keep=[g for g in adapt if g in E.index and g not in rem]
    s=E.loc[keep].mean(0)
    r=spearmanr(bf.values,s.values)[0]
    lmo_rows.append({"markers_removed":n,"genes_remaining":len(keep),"rho":round(r,3)})
    print(f"LMO n={n}: genes={len(keep)} rho={r:.3f}")
pd.DataFrame(lmo_rows).to_csv(rf"{D}\results\tables\leave_marker_stepwise.csv",index=False)

# K null (matched)
ti_in=[g for g in ti if g in E.index]
ti_mean=E.loc[ti_in].mean(1)
adapt_keep150=[g for g in adapt if g in E.index and g not in set(markers[:150])]
adapt_mean=E.loc[adapt_keep150].mean(1)
bins=np.quantile(adapt_mean,[0,0.2,0.4,0.6,0.8,1.0])
adapt_bin=np.digitize(adapt_mean,bins[1:-1])
ti_bin=np.digitize(ti_mean,bins[1:-1])
n_target=len(adapt_keep150)
rng=np.random.default_rng(42)
null=[]
for rep in range(1000):
    picked=[]
    for b in range(5):
        n_b=(adapt_bin==b).sum()
        pool=np.where(ti_bin==b)[0]
        if len(pool)==0: pool=np.arange(len(ti_in))
        picked.extend(rng.choice(pool,min(n_b,len(pool)),replace=False).tolist())
    genes=[ti_in[i] for i in picked[:n_target]]
    s=E.loc[genes].mean(0)
    null.append(spearmanr(bf.values,s.values)[0])
null=np.array(null)
print(f"\nK null: median={np.median(null):.3f} 95%=[{np.percentile(null,2.5):.3f},{np.percentile(null,97.5):.3f}]")
obs=spearmanr(bf.values,E.loc[adapt_keep150].mean(0).values)[0]
print(f"observed={obs:.3f} p={(np.sum(null>=obs)+1)/(len(null)+1):.4f}")
pd.DataFrame({"matched_null":null}).to_csv(rf"{D}\results\tables\moduleK_null.csv",index=False)
print("saved")
