"""
05_classify.py
Decompose brain-metastasis acquired genes into:
  - Competence: present in primary BP and retained/amplified in BM
  - Adaptation: low in primary BP, strongly induced upon brain colonization
Then assign the 2x2 biological classification.
"""
import os, sys, re
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS

# ---- load paired expression (same pairing as 02) ----
df = pd.read_excel(os.path.join(RAW, "GSE184869_expression.xlsx")).rename(
    columns={"ensembl_gene_id":"ensembl"}).set_index("ensembl")
de = pd.read_csv(os.path.join(RESULTS, "tables", "paired_de_gse184869.csv"))

def classify(col):
    c = col.strip()
    m = re.match(r"^MAYO_B([MP])_(\d+)$", c)
    if m: return f"mayo_{m.group(2)}", "BM" if m.group(1)=="M" else "BP"
    m = re.match(r"^B([MP])([\d\-]+)$", c)
    if m: return f"rm_{m.group(2)}", "BM" if m.group(1)=="M" else "BP"
    m = re.match(r"^(\d+)M_RCS$", c);
    if m: return f"rcs_{m.group(1)}", "BM"
    m = re.match(r"^(\d+)P_RCS$", c)
    if m: return f"rcs_{m.group(1)}", "BP"
    return None, None

pm = {}
for c in df.columns:
    p,s = classify(c)
    if p: pm.setdefault(p,{})[s]=c
complete = {p:v for p,v in pm.items() if "BP" in v and "BM" in v}
bp_cols = [complete[p]["BP"] for p in sorted(complete)]
bm_cols = [complete[p]["BM"] for p in sorted(complete)]

mean_bp = df[bp_cols].mean(axis=1)
mean_bm = df[bm_cols].mean(axis=1)

res = de.set_index("ensembl").copy()
res["mean_BP"] = mean_bp
res["mean_BM"] = mean_bm

# distribution of baseline expression in primary tumors
bp_median = mean_bp.median()
bp_q25 = mean_bp.quantile(0.25)
bp_q10 = mean_bp.quantile(0.10)
print(f"BP expression: median={bp_median:.2f} Q25={bp_q25:.2f} Q10={bp_q10:.2f}")

# acquired-up set
up = res[(res.FDR < 0.05) & (res.mean_logFC_BM_vs_BP > 0.3)].copy()
print(f"acquired-up genes: {len(up)}")

# Competence vs Adaptation:
#  Adaptation: low baseline in primary (mean_BP <= Q25) AND strong induction (logFC >= 0.7)
#  Competence: already present in primary (mean_BP >= median) AND up/retained in BM (logFC>0.3)
#  Intermediate: ambiguous
up["baseline_in_primary"] = pd.cut(up.mean_BP, bins=[-np.inf, bp_q25, bp_median, np.inf],
                                   labels=["low(Q1-25)","mid","high(>median)"])
def row_class(r):
    if r.mean_BP <= bp_q25 and r.mean_logFC_BM_vs_BP >= 0.7:
        return "Adaptation_de_novo"
    if r.mean_BP >= bp_median and r.mean_logFC_BM_vs_BP > 0.3:
        return "Competence_preexisting"
    return "Intermediate"
up["program"] = up.apply(row_class, axis=1)
print(up["program"].value_counts())

# also down-regulated in BM (loss of primary program)
down = res[(res.FDR<0.05) & (res.mean_logFC_BM_vs_BP < -0.3)].copy()
print(f"acquired-down (loss in BM): {len(down)}")

out = up.sort_values("mean_logFC_BM_vs_BP", ascending=False)
out.to_csv(os.path.join(RESULTS, "tables", "program_decomposition.csv"))

print("\n=== ADAPTATION (de novo induced; low in primary, high in BM) top 25 ===")
adapt = out[out.program=="Adaptation_de_novo"].sort_values("mean_logFC_BM_vs_BP", ascending=False)
print(adapt[["symbol","mean_BP","mean_BM","mean_logFC_BM_vs_BP","FDR"]].head(25).to_string(index=False))

print("\n=== COMPETENCE (pre-existing in primary, retained in BM) top 25 by logFC ===")
comp = out[out.program=="Competence_preexisting"].sort_values("mean_logFC_BM_vs_BP", ascending=False)
print(comp[["symbol","mean_BP","mean_BM","mean_logFC_BM_vs_BP","FDR"]].head(25).to_string(index=False))

# save gene lists
for name, sub in [("adaptation", adapt), ("competence", comp), ("intermediate", out[out.program=="Intermediate"])]:
    syms = sub.symbol.dropna().unique().tolist()
    with open(os.path.join(RESULTS, "tables", f"geneset_{name}.txt"), "w") as f:
        f.write("\n".join(syms))
    print(f"saved {len(syms)} symbols -> geneset_{name}.txt")
