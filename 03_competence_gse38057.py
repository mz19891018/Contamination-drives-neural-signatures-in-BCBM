"""03_competence_gse38057.py
Identify genes in PRIMARY breast tumors that distinguish patients who later
developed brain metastasis (competence signal) from those who did not.
GSE38057: 87 HER2+ primaries, sample titles encode outcome.
"""
import os, sys, re
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

expr, meta = read_series_matrix(os.path.join(RAW, "GSE38057_series_matrix.txt.gz"))
print("expr shape:", expr.shape)
print("meta cols:", list(meta.columns)[:10])

# outcome from sample title
title_col = [c for c in meta.columns if "title" in c.lower()][0]
titles = meta[title_col]
def outcome(t):
    if "no brain metastasis" in t: return "noBM"
    if "brain metastasis" in t: return "BM"
    return np.nan
meta["outcome"] = titles.map(outcome)
print(meta["outcome"].value_counts())

# align
common = [s for s in expr.columns if s in meta.index]
expr = expr[common]
g = meta.loc[common, "outcome"]
bm_cols = g[g=="BM"].index
nb_cols = g[g=="noBM"].index
print(f"BM={len(bm_cols)} noBM={len(nb_cols)}")

# rows may be probe symbols; inspect
print("row names head:", list(expr.index[:10]))

X_bm = expr[bm_cols].values
X_nb = expr[nb_cols].values
# welch t-test
t, p = stats.ttest_ind(X_bm, X_nb, axis=1, equal_var=False, nan_policy="omit")
logfc = X_bm.mean(axis=1) - X_nb.mean(axis=1)
FDR = multipletests(p, method="fdr_bh")[1]
out = pd.DataFrame({"probe": expr.index, "mean_BM": X_bm.mean(axis=1),
                    "mean_noBM": X_nb.mean(axis=1),
                    "logFC_BMvsNoBM_primary": logfc, "p": p, "FDR": FDR})
out = out.sort_values("p")
out.to_csv(os.path.join(RESULTS, "tables", "competence_gse38057.csv"), index=False)
print("\n=== top competence genes (primary predicts later BM) ===")
print(out.head(30).to_string(index=False))
print(f"\nFDR<0.05: {(out.FDR<0.05).sum()}  | up in BM-developers: {((out.FDR<0.05)&(out.logFC_BMvsNoBM_primary>0)).sum()}")
