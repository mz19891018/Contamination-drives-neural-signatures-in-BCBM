# 52_batch_baseline_paired.py
import pandas as pd, numpy as np, openpyxl, re
from scipy.stats import spearmanr, kruskal, pearsonr, mannwhitneyu, wilcoxon
from scipy.stats import rankdata
D=r"D:\BCBM_Project"

# load fractions
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
def kind(s):
    if re.match(r"^BM\d",s): return "BM_named"
    if re.match(r"^\d+M_RCS",s): return "BM_RCS"
    if re.match(r"^MAYO",s): return "BM_MAYO"
    if re.match(r"^BP\d",s): return "P_named"
    if re.match(r"^\d+P_RCS",s): return "P_RCS"
    return "Other"
th["batch"]=[kind(s) for s in th.index]
th["brain"]=th[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)
bm=th[th.batch.str.startswith("BM")]
print("=== Batch distribution of brain fraction ===")
print(bm.groupby("batch")["brain"].agg(["count","mean","median","std"]).round(3).to_string())
# Kruskal-Wallis
groups=[bm[bm.batch==b]["brain"].values for b in ["BM_named","BM_RCS","BM_MAYO"]]
kw=kruskal(*groups)
print(f"\nKruskal-Wallis across 3 BM batches: H={kw.statistic:.2f}, p={kw.pvalue:.4f}")

# load expression for all BM samples
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("BM")]
bm_names=[header[i] for i in bm_idx]
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv"); emap=dict(zip(pm.ensembl,pm.symbol))
# also primary baseline
prim_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("P")]
prim_names=[header[i] for i in prim_idx]
expr={}; baseline={}
for r in rows:
    sym=emap.get(r[0])
    if sym is None: continue
    expr[sym]=[r[i] for i in bm_idx]
    baseline[sym]=np.mean([r[i] for i in prim_idx])
wb.close()
E=pd.DataFrame(expr,index=bm_names).T
common=[s for s in bm_names if s in th.index]
E=E[common]; bf=th.loc[common,"brain"]
batches=th.loc[common,"batch"]
print(f"\nBM samples with expression+fractions: {len(common)}")

# rho per batch
print("\n=== rho per batch (adaptation score vs brain fraction) ===")
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
adapt_in=[g for g in adapt if g in E.index]
for b in ["BM_named","BM_RCS","BM_MAYO"]:
    mask=batches==b
    if mask.sum()<5: 
        print(f"  {b}: n={mask.sum()}, too few"); continue
    s=E.loc[adapt_in,mask].mean(0)
    r,p=spearmanr(bf[mask],s)
    print(f"  {b}: n={mask.sum()}, rho={r:.3f}, p={p:.2e}")
# all BM
s_all=E.loc[adapt_in].mean(0)
r_all,p_all=spearmanr(bf,s_all)
print(f"  ALL BM: n={len(bf)}, rho={r_all:.3f}, p={p_all:.2e}")

# partial correlation controlling for batch (rank-based: residualize ranks on batch dummies)
def partial_spearman(x,y,cov):
    # cov is categorical -> use dummy residuals
    import numpy as np
    xr=rankdata(x); yr=rankdata(y)
    dummies=pd.get_dummies(cov,drop_first=True).values.astype(float)
    X=np.column_stack([np.ones(len(x)),dummies])
    bx=np.linalg.lstsq(X,xr,rcond=None)[0]; resid_x=xr-X@bx
    by=np.linalg.lstsq(X,yr,rcond=None)[0]; resid_y=yr-X@by
    return pearsonr(resid_x,resid_y)
pr,pp=partial_spearman(bf.values,s_all.values,batches.values)
print(f"\nPartial Spearman (controlling batch): rho={pr:.3f}, p={pp:.2e}")

# rho vs primary baseline expression (per gene)
print("\n=== rho vs primary baseline expression ===")
H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
H["baseline"]=H.gene.map(baseline)
H=H.dropna(subset=["baseline"])
r_bl,p_bl=spearmanr(H.baseline,H.rho_brain)
print(f"Spearman(baseline, rho_brain) = {r_bl:.3f}, p={p_bl:.2e}, n={len(H)}")
# by baseline quintiles
H["bl_bin"]=pd.qcut(H.baseline,5,labels=False)
print(H.groupby("bl_bin")["rho_brain"].agg(["count","mean","median"]).round(3).to_string())

# baseline-matched comparison: brain-specific (contamination_driven) vs non-brain (tumor_intrinsic)
cont=H[H.verdict=="contamination_driven"]
ti=H[H.verdict=="tumor_intrinsic_candidate"]
# match on baseline: for each cont gene, find ti gene with closest baseline
matched_ti=[]
for _,row in cont.iterrows():
    diffs=(ti.baseline-row.baseline).abs()
    matched_ti.append(ti.loc[diffs.idxmin(),"rho_brain"])
print(f"\n=== Baseline-matched comparison ===")
print(f"contamination-driven (n={len(cont)}): mean rho={cont.rho_brain.mean():.3f}, median={cont.rho_brain.median():.3f}")
print(f"baseline-matched tumor-intrinsic (n={len(matched_ti)}): mean rho={np.mean(matched_ti):.3f}, median={np.median(matched_ti):.3f}")
u,p=mannwhitneyu(cont.rho_brain,matched_ti,alternative="greater")
print(f"Mann-Whitney U (cont > matched TI): p={p:.2e}")
print(f"baseline: cont mean={cont.baseline.mean():.2f}, matched TI mean={ti.loc[ti.gene.isin([ti.loc[(ti.baseline-row.baseline).abs().idxmin(),'gene'] for _,row in cont.iterrows()]),'baseline'].mean():.2f}")

# paired delta test (21 paired samples)
print("\n=== Paired BM vs Primary delta (21 pairs) ===")
pairs=[]
for i in range(1,16):
    bm_s=f"BM{i}"; p_s=f"BP{i}"
    if bm_s in th.index and p_s in th.index: pairs.append((bm_s,p_s))
for i in range(1,7):
    bm_s=f"{i}M_RCS"; p_s=f"{i}P_RCS"
    if bm_s in th.index and p_s in th.index: pairs.append((bm_s,p_s))
print(f"paired samples: {len(pairs)}")
bm_brain=[th.loc[b,"brain"] for b,p in pairs]
p_brain=[th.loc[p,"brain"] for b,p in pairs]
print(f"BM mean brain={np.mean(bm_brain):.3f}, Primary mean={np.mean(p_brain):.3f}")
w,pw=wilcoxon(bm_brain,p_brain)
print(f"Wilcoxon paired: W={w}, p={pw:.4f}")
deltas=[b-p for b,p in zip(bm_brain,p_brain)]
print(f"delta mean={np.mean(deltas):.3f}, median={np.median(deltas):.3f}, range=[{min(deltas):.3f},{max(deltas):.3f}]")

# save key results
res={
 "batch_kw_p":kw.pvalue,"rho_all":r_all,"rho_all_p":p_all,
 "partial_rho":pr,"partial_p":pp,
 "baseline_rho_corr":r_bl,"baseline_p":p_bl,
 "paired_n":len(pairs),"paired_wilcoxon_p":pw,"paired_delta_mean":np.mean(deltas)
}
pd.Series(res).to_csv(rf"{D}\results\tables\batch_baseline_results.csv")
print("\nsaved")
