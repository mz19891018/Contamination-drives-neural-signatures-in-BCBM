"""12_summary.py - gather key stats for the report."""
import os, sys
import pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import RESULTS

de = pd.read_csv(os.path.join(RESULTS,"tables","paired_de_gse184869.csv"))
prog = pd.read_csv(os.path.join(RESULTS,"tables","program_decomposition.csv"))

print("=== PAIRED DE (GSE184869, 45 pairs) ===")
print("genes tested:", len(de))
print("up FDR<0.05 logFC>0.3:", ((de.FDR<0.05)&(de.mean_logFC_BM_vs_BP>0.3)).sum())
print("down FDR<0.05 logFC<-0.3:", ((de.FDR<0.05)&(de.mean_logFC_BM_vs_BP<-0.3)).sum())
print("\n=== PROGRAM DECOMPOSITION ===")
print(prog.program.value_counts())
print("\n=== TOP ADAPTATION (de novo) ===")
a = prog[prog.program=="Adaptation_de_novo"].sort_values("mean_logFC_BM_vs_BP",ascending=False)
print(a[["symbol","mean_BP","mean_BM","mean_logFC_BM_vs_BP"]].head(15).to_string(index=False))
print("\n=== TOP COMPETENCE (pre-existing, high in primary) ===")
c = prog[prog.program=="Competence_preexisting"].sort_values("mean_BP",ascending=False)
print(c[["symbol","mean_BP","mean_BM","mean_logFC_BM_vs_BP"]].head(15).to_string(index=False))
