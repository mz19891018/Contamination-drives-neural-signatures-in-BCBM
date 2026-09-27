# 35_Anew_metrics.py - replace Spearman with amplitude-sensitive metrics
import pandas as pd, numpy as np
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

# observed logFC
prog=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
obs=dict(zip(prog.symbol, prog.mean_logFC_BM_vs_BP))
up=prog[(prog.program!="Intermediate")&(prog.mean_logFC_BM_vs_BP>0.3)]
up_genes=up.symbol.tolist()
print("up genes:",len(up_genes))

# simulated logFC per f from mixing_simulation.csv only has 3 genes.
# Need full simulated logFC vector per f. Re-run a lightweight version to get full vector.
# Instead, recompute from the simulation script output? We only saved 3 genes.
# Let's re-run a compact simulation saving full logFC for each f (1 rep only, fast).
import subprocess, sys
print("need full sim logFC; will compute in R then read")
