# 34_bp_combined_analyze.py
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
def kind(s):
    if re.match(r"^BM",s): return "BrainMet"
    if re.match(r"^BP",s): return "Primary"
    return "Other"
th["kind"]=[kind(i) for i in th.index]
brain=["Neuron","Astrocyte","Oligodendrocyte"]
th["brain_nontumor"]=th[brain].sum(1)
th["tumor"]=th["Tumor"]
print(th.groupby("kind")[["tumor","brain_nontumor"]].mean().round(3).to_string())

cont=pd.read_csv(rf"{D}\results\tables\contamination_scores_bulk.csv",index_col=0)
cont=cont[[c for c in cont.columns if c in ("adapt_mod","comp_mod")]]
m=th.join(cont)
bm=m[m.kind=="BrainMet"]
r,p=spearmanr(bm.brain_nontumor,bm.adapt_mod)
print(f"\nBM: adaptation vs combined-ref brain fraction: rho={r:.3f} p={p:.2e}")
r2,p2=spearmanr(bm.brain_nontumor,bm.comp_mod)
print(f"BM: competence vs combined-ref brain fraction: rho={r2:.3f} p={p2:.3f}")
th.to_csv(rf"{D}\results\tables\bp_combined_with_groups.csv")
print("saved")
