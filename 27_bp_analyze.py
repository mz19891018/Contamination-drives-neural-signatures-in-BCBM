# 27_bp_analyze.py
import pandas as pd, numpy as np, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"
theta=pd.read_csv(rf"{D}\results\tables\bp_cell_fractions.csv",index_col=0)
print("fraction columns:",list(theta.columns))

# classify samples
def kind(s):
    if re.match(r"^BM",s): return "BrainMet"
    if re.match(r"^BP",s): return "Primary"
    return "Other"
theta["kind"]=[kind(i) for i in theta.index]
cols=[c for c in theta.columns if c!="kind"]
print("\n=== mean fraction by sample group ===")
g=theta.groupby("kind")[cols].mean().round(3)
print(g.to_string())

brain_cols=[c for c in cols if c in ("Neurons","Astrocytes","ODG")]
tumor_cols=[c for c in cols if "Tumor" in c]
theta["brain_nontumor"]=theta[brain_cols].sum(1)
theta["tumor_total"]=theta[tumor_cols].sum(1)

print("\n=== key fractions by group ===")
print(theta.groupby("kind")[["tumor_total","brain_nontumor"]].mean().round(3).to_string())

# correlation: adaptation score per sample vs brain_nontumor fraction
cont=pd.read_csv(rf"{D}\results\tables\contamination_scores_bulk.csv",index_col=0)
common=theta.index.intersection(cont.index)
m=theta.loc[common].join(cont[["adapt_mod","comp_mod"]])
print("\n=== in BrainMet: module vs BP-estimated brain fraction ===")
bm=m[m.kind=="BrainMet"]
for mod in ["adapt_mod","comp_mod"]:
    r,p=spearmanr(bm.brain_nontumor,bm[mod])
    print(f"{mod} vs BayesPrism brain(neuron+astro+ODG) fraction: rho={r:+.3f} p={p:.4g}")
theta.to_csv(rf"{D}\results\tables\bp_fractions_with_groups.csv")
print("\nsaved.")
