# 24_permute_decomp.py - regression-to-mean + label-shuffle permutation
import pandas as pd, numpy as np, re
D = r"D:\BCBM_Project"
expr = pd.read_excel(rf"{D}\data\raw\GSE184869_expression.xlsx", sheet_name="log2TMMCPM").set_index("ensembl_gene_id")

# find BM/BP pairs
bm_cols = [c for c in expr.columns if re.match(r"^BM\d", c)]
bp_cols = [c for c in expr.columns if re.match(r"^BP\d", c)]
def key(c): return re.sub(r"^(BM|BP)","",c)
pairs = []
for bm in bm_cols:
    k = key(bm)
    match = [bp for bp in bp_cols if key(bp)==k]
    if match: pairs.append((bm,match[0]))
print("paired samples:", len(pairs), flush=True)

BM = expr[[p[0] for p in pairs]].values
BP = expr[[p[1] for p in pairs]].values
logfc = (BM - BP).mean(1)
base = BP.mean(1)

# 1) regression-to-mean: corr(baseline, |logFC|)
from scipy.stats import spearmanr
r,p = spearmanr(base, np.abs(logfc))
print(f"corr(primary baseline, |logFC|) = {r:.3f}  p={p:.2e}", flush=True)

# observed de-novo count by our rule: FDR<0.05 up, primary baseline low
# approximate: genes with logFC>0.7 AND baseline < 25th percentile
q25 = np.percentile(base,25)
obs_de_novo = int(np.sum((logfc>0.7) & (base<q25)))
obs_comp = int(np.sum((logfc>0.3) & (base>np.median(base))))
print(f"observed de-novo={obs_de_novo}, competence={obs_comp}", flush=True)

# 2) permutation: swap BM/BP within patient across random subsets -> null
rng = np.random.default_rng(0)
null_de = []
for _ in range(200):
    swap = rng.random(len(pairs)) < 0.5
    lf = np.where(swap[None,:], BP-BM, BM-BP).mean(1)
    bl = np.where(swap[None,:], BM, BP).mean(1)
    q = np.percentile(bl,25)
    null_de.append(np.sum((lf>0.7)&(bl<q)))
null_de = np.array(null_de)
print(f"permuted de-novo under random label: mean={null_de.mean():.1f} max={null_de.max()} p~{(np.sum(null_de>=obs_de_novo)+1)/201:.4f}", flush=True)

# 3) threshold sensitivity
for lfc in [0.5,0.7,1.0]:
    for qthr in [20,25,33]:
        q = np.percentile(base,qthr)
        dn = np.sum((logfc>lfc)&(base<q))
        print(f"  logFC>{lfc}, baseline<Q{qthr}: de_novo={dn}", flush=True)
