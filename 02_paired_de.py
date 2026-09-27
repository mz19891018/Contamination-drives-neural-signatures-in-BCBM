"""
02_paired_de.py
Paired differential expression: Primary (BP/P) vs Brain Metastasis (BM/M)
across 45 patients in GSE184869.

Output:
  results/tables/paired_de_gse184869.csv  (gene-level logFC, paired t p, FDR)
  results/figures/paired_volcano.png
"""
import os, sys, re
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, ensembl_to_symbol
from statsmodels.stats.multitest import multipletests

os.makedirs(os.path.join(RESULTS, "tables"), exist_ok=True)
os.makedirs(os.path.join(RESULTS, "figures"), exist_ok=True)

# ---------- load ----------
df = pd.read_excel(os.path.join(RAW, "GSE184869_expression.xlsx"))
df = df.rename(columns={"ensembl_gene_id": "ensembl"}).set_index("ensembl")
print("matrix:", df.shape)

# ---------- build pairs ----------
def classify(col):
    """Return (patient_id, side) where side in {'BP','BM'}."""
    c = col.strip()
    m = re.match(r"^MAYO_B([MP])_(\d+)$", c)
    if m:
        return f"mayo_{m.group(2)}", "BM" if m.group(1) == "M" else "BP"
    m = re.match(r"^B([MP])([\d\-]+)$", c)
    if m:
        return f"rm_{m.group(2)}", "BM" if m.group(1) == "M" else "BP"
    m = re.match(r"^(\d+)M_RCS$", c)
    if m:
        return f"rcs_{m.group(1)}", "BM"
    m = re.match(r"^(\d+)P_RCS$", c)
    if m:
        return f"rcs_{m.group(1)}", "BP"
    return None, None

pair_map = {}   # patient -> {'BP':col,'BM':col}
for col in df.columns:
    pid, side = classify(col)
    if pid is None:
        print("unmatched:", col); continue
    pair_map.setdefault(pid, {})[side] = col

# keep only complete pairs
complete = {p: v for p, v in pair_map.items() if "BP" in v and "BM" in v}
print(f"complete pairs: {len(complete)}")

bp_cols = [complete[p]["BP"] for p in sorted(complete)]
bm_cols = [complete[p]["BM"] for p in sorted(complete)]
patients = sorted(complete)

# ---------- paired test ----------
bp = df[bp_cols].values   # genes x patients
bm = df[bm_cols].values
delta = bm - bp           # log2 scale already; positive = up in BM

mean_d = delta.mean(axis=1)
# paired t-test
t, pval = stats.ttest_rel(bm, bp, axis=1, nan_policy="omit")
# also Wilcoxon for robustness
w_p = np.array([stats.wilcoxon(delta[i], zero_method="wilcox").pvalue
                if np.isfinite(delta[i]).sum() >= 3 else np.nan
                for i in range(delta.shape[0])])

FDR = multipletests(pval, method="fdr_bh")[1]
out = pd.DataFrame({
    "ensembl": df.index,
    "mean_logFC_BM_vs_BP": mean_d,
    "paired_t_p": pval,
    "wilcoxon_p": w_p,
    "FDR": FDR,
    "frac_up": (delta > 0).mean(axis=1),
})
# annotate symbols
sym = ensembl_to_symbol(list(out["ensembl"]))
out["symbol"] = out["ensembl"].map(sym)
out = out.sort_values("paired_t_p")
out.to_csv(os.path.join(RESULTS, "tables", "paired_de_gse184869.csv"), index=False)

print("\n=== top 25 up in brain metastasis ===")
up = out[(out.FDR < 0.05) & (out.mean_logFC_BM_vs_BP > 0.3)].sort_values("mean_logFC_BM_vs_BP", ascending=False)
print(up[["symbol","mean_logFC_BM_vs_BP","FDR"]].head(25).to_string(index=False))
print("\n=== top 25 down in brain metastasis ===")
dn = out[(out.FDR < 0.05) & (out.mean_logFC_BM_vs_BP < -0.3)].sort_values("mean_logFC_BM_vs_BP")
print(dn[["symbol","mean_logFC_BM_vs_BP","FDR"]].head(25).to_string(index=False))
print(f"\nTotal genes tested: {len(out)}")
print(f"Up (FDR<0.05 & logFC>0.3): {len(up)}")
print(f"Down (FDR<0.05 & logFC<-0.3): {len(dn)}")

# save pair list for downstream
pd.DataFrame({"patient": patients, "BP_col": bp_cols, "BM_col": bm_cols}).to_csv(
    os.path.join(RESULTS, "tables", "gse184869_pairs.csv"), index=False)
