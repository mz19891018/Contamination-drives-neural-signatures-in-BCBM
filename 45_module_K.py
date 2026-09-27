# 45_module_K.py - negative control for B-fix
import pandas as pd, numpy as np, openpyxl, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
ti=H[H.verdict=="tumor_intrinsic_candidate"].gene.tolist()  # 973
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
th["kind"]=["BrainMet" if re.match(r"^BM|^\d+M_RCS|^MAYO",s) else "Primary" if re.match(r"^BP|^\d+P_RCS",s) else "Other" for s in th.index]
bm=th[th.kind=="BrainMet"]
brain=bm[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)

# load expression for all genes in H
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and re.match(r"^BM|^\d+M_RCS|^MAYO",str(h))]
bm_names=[header[i] for i in bm_idx]
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv"); emap=dict(zip(pm.ensembl,pm.symbol))
expr={}
for r in rows:
    sym=emap.get(r[0])
    if sym in set(ti)|set(adapt): expr[sym]=[r[i] for i in bm_idx]
wb.close()
E=pd.DataFrame(expr,index=bm_names).T
E=E[[s for s in bm_names if s in brain.index]]; bf=brain.loc[E.columns]
print("expr:",E.shape," BM:",len(bf))

# observed: adaptation after removing top-150 markers
markers=H[H.verdict=="contamination_driven"].sort_values("rho_brain",ascending=False).gene.tolist()
rem=set(markers[:150])
adapt_keep=[g for g in adapt if g in E.index and g not in rem]
obs_score=E.loc[adapt_keep].mean(0)
obs_rho=spearmanr(bf.values,obs_score.values)[0]
print(f"observed (n={len(adapt_keep)} genes): rho={obs_rho:.3f}")

# expression-matched random sampling from tumor-intrinsic candidates
ti_in=[g for g in ti if g in E.index]
ti_mean=E.loc[ti_in].mean(1)
adapt_mean=E.loc[adapt_keep].mean(1)
# 5 expression bins
bins=np.quantile(adapt_mean,[0,0.2,0.4,0.6,0.8,1.0])
adapt_bin=np.digitize(adapt_mean,bins[1:-1])
ti_bin=np.digitize(ti_mean,bins[1:-1])
n_target=len(adapt_keep)

rng=np.random.default_rng(42)
null=[]
for rep in range(1000):
    picked=[]
    for b in range(5):
        n_b=(adapt_bin==b).sum()
        pool=np.where(ti_bin==b)[0]
        if len(pool)==0: pool=np.where(ti_bin!=99)[0]
        picked.extend(rng.choice(pool,min(n_b,len(pool)),replace=False).tolist())
    genes=[ti_in[i] for i in picked[:n_target]]
    s=E.loc[genes].mean(0)
    null.append(spearmanr(bf.values,s.values)[0])
null=np.array(null)
print(f"\n=== K null (1000 matched random sets from tumor-intrinsic, n={n_target}) ===")
print(f"null median={np.median(null):.3f}  95%=[{np.percentile(null,2.5):.3f},{np.percentile(null,97.5):.3f}]")
print(f"observed rho={obs_rho:.3f}  empirical p={(np.sum(null>=obs_rho)+1)/(len(null)+1):.4f}")

# unmatched version
null2=[]
for rep in range(1000):
    genes=rng.choice(ti_in,n_target,replace=False)
    s=E.loc[genes].mean(0)
    null2.append(spearmanr(bf.values,s.values)[0])
null2=np.array(null2)
print(f"\nunmatched null median={np.median(null2):.3f}  95%=[{np.percentile(null2,2.5):.3f},{np.percentile(null2,97.5):.3f}]")
print(f"unmatched empirical p={(np.sum(null2>=obs_rho)+1)/(len(null2)+1):.4f}")

pd.DataFrame({"matched_null":null,"unmatched_null":null2}).to_csv(rf"{D}\results\tables\moduleK_null.csv",index=False)
print("saved")
