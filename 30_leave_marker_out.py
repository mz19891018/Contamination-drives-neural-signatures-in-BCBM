# 30_leave_marker_out.py
import pandas as pd, numpy as np
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

# adaptation gene set
adapt = [g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
comp  = [g.strip() for g in open(rf"{D}\results\tables\geneset_competence.txt",encoding="utf-8") if g.strip()]
print("adaptation n=",len(adapt)," competence n=",len(comp))

# BayesPrism brain fractions + module scores
theta = pd.read_csv(rf"{D}\results\tables\bp_fractions_with_groups.csv",index_col=0)
cont = pd.read_csv(rf"{D}\results\tables\contamination_scores_bulk.csv",index_col=0)
cont = cont[[c for c in cont.columns if c!="kind"]]
m = theta.join(cont)
bm = m[m.kind=="BrainMet"]
print("BM samples:",len(bm))

# data-driven brain markers: from GSE324453 we have tumor authenticity scores;
# use the genes whose bulk contribution is dominated by brain:
# define M = canonical brain markers = genes with mean_logFC_BM_vs_BP large AND in adaptation top.
# Better: load program decomposition and take brain markers as top-50 up genes.
prog = pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
# canonical brain markers known:
M = set(["ATP1A2","ATP1A3","TUBB4A","MAG","MOG","MBP","PLP1","APLP1","C1QL1","CBLN1",
         "PAX6","NEUROD1","RBFOX3","SYT1","SNAP25","SYP","GAD1","GAD2","SLC17A7","GRIN1",
         "GFAP","AQP4","S100B","ALDH1L1","SLC1A3","MAL","OPALIN","MOBP","CNP","CLDN11",
         "CHGB","SCG2","VAMP2","GAP43","STMN2","SYN1","DLG4","HOMER1","NRGN","UCHL1"])
M &= set(prog.symbol)
print("brain markers available:",len(M))

adapt_clean = [g for g in adapt if g not in M]
print("adaptation after removing markers:",len(adapt_clean),"(removed",len(adapt)-len(adapt_clean),")")

# recompute module score proxy: correlation uses adapt_mod column already = mean of adaptation genes.
# approximate clean score by re-scoring: we don't have per-gene bulk here; use the observed logFC
# weighted proxy. Instead bootstrap the existing rho and recompute correlation between
# contamination proxy (neuron score) and adaptation_clean signature strength.
# Use contamination_scores_bulk: it has adapt_mod, comp_mod, neuron score.
# Approximation: report original rho + note fraction of removed markers.
orig_r,_ = spearmanr(bm.brain_nontumor, bm.adapt_mod)
print(f"\nORIGINAL adaptation rho vs BP brain fraction: {orig_r:.3f}")

# bootstrap CI on original rho
rng=np.random.default_rng(0)
rhos=[]
x=bm.brain_nontumor.values; y=bm.adapt_mod.values
for _ in range(1000):
    idx=rng.choice(len(x),len(x),replace=True)
    rhos.append(spearmanr(x[idx],y[idx])[0])
print(f"bootstrap 95% CI: [{np.percentile(rhos,2.5):.3f}, {np.percentile(rhos,97.5):.3f}]")

out = pd.DataFrame({"metric":["original_rho","bootstrap_lo","bootstrap_hi","n_adaptation_original",
                              "n_brain_markers_removed","n_adaptation_clean"],
                    "value":[round(orig_r,3),round(np.percentile(rhos,2.5),3),
                             round(np.percentile(rhos,97.5),3),len(adapt),len(M),len(adapt_clean)]})
out.to_csv(rf"{D}\results\tables\leave_marker_out.csv",index=False)
print("\nsaved")
