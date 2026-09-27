"""debug pairing."""
import os, sys, re
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW

df = pd.read_excel(os.path.join(RAW, "GSE184869_expression.xlsx")).rename(
    columns={"ensembl_gene_id":"ensembl"}).set_index("ensembl")

def classify(col):
    c = col.strip()
    m = re.match(r"^MAYO_B([BP])_(\d+)$", c)
    if m: return f"mayo_{m.group(2)}", "BM" if m.group(1)=="M" else "BP"
    m = re.match(r"^B([BP])([\d\-]+)$", c)
    if m: return f"rm_{m.group(2)}", "BM" if m.group(1)=="M" else "BP"
    m = re.match(r"^(\d+)M_RCS$", c)
    if m: return f"rcs_{m.group(1)}", "BM"
    m = re.match(r"^(\d+)P_RCS$", c)
    if m: return f"rcs_{m.group(1)}", "BP"
    return None, None

pm = {}
for c in df.columns:
    p,s = classify(c)
    if p: pm.setdefault(p,{})[s]=c
complete = {p:v for p,v in pm.items() if "BP" in v and "BM" in v}
print("complete pairs:", len(complete))
print("unmatched patients:", [p for p,v in pm.items() if not ("BP" in v and "BM" in v)])

bp = df[[complete[p]["BP"] for p in sorted(complete)]]
bm = df[[complete[p]["BM"] for p in sorted(complete)]]
delta = (bm.values - bp.values)
print("delta shape:", delta.shape)
print("mean |delta| per gene:", np.abs(delta).mean(axis=1).describe() if hasattr(np.abs(delta).mean(axis=1),'describe') else "")
md = np.abs(delta).mean(axis=1)
print("delta abs: median=%.3f p90=%.3f p99=%.3f max=%.3f" % (np.median(md), np.percentile(md,90), np.percentile(md,99), md.max()))
# Is there ANY gene with consistent direction?
frac_up = (delta>0).mean(axis=1)
print("frac_up distribution: median=%.3f, genes with >70%% up: %d, >70%% down: %d" % (
    np.median(frac_up), (frac_up>0.7).sum(), (frac_up<0.3).sum()))
# overall mean difference
print("mean delta across all genes:", delta.mean())
# check a known brain-metastasis gene: e.g. GAPDH housekeeping vs a neural gene
sym_guess = df.index  # ensembl; print top variable genes
import numpy as np
v = df.var(axis=1)
print("most variable genes (ensembl):")
print(v.sort_values(ascending=False).head(10))
