"""
06_organ_specificity.py
Project adaptation vs competence gene sets onto multi-organ metastasis cohorts
(GSE14017 U133plus2, GSE14018 U133A) and test brain specificity.
"""
import os, sys, gzip
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

def load_platform_map(gpl):
    """Parse GPL*.annot.gz -> DataFrame probe->symbol."""
    path = os.path.join(RAW, f"{gpl}.annot.gz")
    rows = {}
    with gzip.open(path, "rt", errors="replace", encoding="utf-8") as f:
        header = None
        for line in f:
            if line.startswith(("#","!","^")):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.strip() for h in parts]
                # find columns
                id_i = next((i for i,h in enumerate(header) if h.upper() in ("ID","PROBE_ID")), 0)
                sym_i = next((i for i,h in enumerate(header) if "GENE SYMBOL" in h.upper() or h.upper()=="SYMBOL"), None)
                continue
            if sym_i is None or len(parts) <= max(id_i, sym_i):
                continue
            pid = parts[id_i]
            syms = parts[sym_i]
            first_sym = syms.split("///")[0].strip()
            if first_sym:
                rows[pid] = first_sym
    return pd.Series(rows, name="symbol")

# load gene sets
adapt = pd.read_csv(os.path.join(RESULTS, "tables", "geneset_adaptation.txt"), header=None)[0].tolist()
comp  = pd.read_csv(os.path.join(RESULTS, "tables", "geneset_competence.txt"), header=None)[0].tolist()
print(f"adaptation set: {len(adapt)} genes; competence set: {len(comp)} genes")

def project_cohort(gse, gpl):
    expr, meta = read_series_matrix(os.path.join(RAW, f"{gse}_series_matrix.txt.gz"))
    pmap = load_platform_map(gpl)
    # map probes -> genes
    expr = expr.copy()
    expr["symbol"] = [pmap.get(i) for i in expr.index]
    expr = expr.dropna(subset=["symbol"])
    expr = expr.groupby("symbol").mean()  # gene-level
    # organ labels
    src = meta["!Sample_source_name_ch1"]
    organ = src.str.strip().str.lower()
    # z-score each gene across samples, then mean score per set
    Z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0,np.nan), axis=0)
    def score(genes):
        g = [x for x in genes if x in Z.index]
        return Z.loc[g].mean(axis=0), len(g)
    adapt_score, n1 = score(adapt)
    comp_score, n2 = score(comp)
    out = pd.DataFrame({"organ": organ, "adapt_score": adapt_score, "comp_score": comp_score})
    print(f"\n--- {gse} ({gpl}) n_genes adapted={n1}, comp={n2} ---")
    print(out.groupby("organ")[["adapt_score","comp_score"]].agg(["mean","count"]).round(3))
    # brain vs others t-test
    br = out[out.organ=="brain"]
    ot = out[out.organ!="brain"]
    for prog in ["adapt_score","comp_score"]:
        t,p = stats.ttest_ind(br[prog], ot[prog], equal_var=False)
        print(f"  {prog}: brain={br[prog].mean():.3f} vs non-brain={ot[prog].mean():.3f}  t={t:.2f} p={p:.4f}")
    return out

o1 = project_cohort("GSE14017","GPL570")
o2 = project_cohort("GSE14018","GPL96")
pd.concat([o1,o2], keys=["GSE14017","GSE14018"]).to_csv(
    os.path.join(RESULTS,"tables","organ_specificity_scores.csv"))
