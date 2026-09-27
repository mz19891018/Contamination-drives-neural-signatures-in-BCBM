"""
08_survival.py
Clinical validation:
  - Competence score -> BMFS / survival in primary cohorts (GSE2603, GSE12276)
  - Program scores projected with platform annotation.
"""
import os, sys, gzip, re
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

# gene sets
adapt = pd.read_csv(os.path.join(RESULTS,"tables","geneset_adaptation.txt"), header=None)[0].tolist()
comp  = pd.read_csv(os.path.join(RESULTS,"tables","geneset_competence.txt"), header=None)[0].tolist()

def load_map(gpl):
    rows={}
    with gzip.open(os.path.join(RAW,f"{gpl}.annot.gz"),"rt",errors="replace",encoding="utf-8") as f:
        header=None
        for line in f:
            if line.startswith(("#","!","^")): continue
            p=line.rstrip("\n").split("\t")
            if header is None: header=p; continue
            rows[p[0]] = p[2].split("///")[0].strip()
    return pd.Series(rows)

def parse_chars(meta):
    """Each sample's !Sample_characteristics_ch1 may be a list of 'key: value' strings."""
    out = {}
    for s, row in meta.iterrows():
        d = {}
        for c in meta.columns:
            if "characteristics" in c.lower():
                v = row[c]
                if isinstance(v,str) and ":" in v:
                    k,val = v.split(":",1)
                    d[k.strip()] = val.strip()
        out[s] = d
    return pd.DataFrame(out).T

def project(expr, meta, gpl):
    pmap = load_map(gpl)
    e = expr.copy()
    e["symbol"] = [pmap.get(i) for i in e.index]
    e = e.dropna(subset=["symbol"]).groupby("symbol").mean()
    Z = e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1).replace(0,np.nan),axis=0)
    cs = Z.loc[[g for g in comp if g in Z.index]].mean(axis=0)
    as_ = Z.loc[[g for g in adapt if g in Z.index]].mean(axis=0)
    sc = pd.DataFrame({"comp_score":cs, "adapt_score":as_})
    clin = parse_chars(meta)
    return sc.join(clin)

# ---------- GSE2603 : bmfs ----------
print("="*60); print("GSE2603 (bmfs)")
expr, meta = read_series_matrix(os.path.join(RAW,"GSE2603_series_matrix.txt.gz"))
df = project(expr, meta, "GPL96")
# keep primary human tumors (exclude cell lines)
src = meta["!Sample_source_name_ch1"]
df["source"] = src
print(df[["comp_score","adapt_score"]].describe())
# extract bmfs
def fnum(x):
    try: return float(x)
    except: return np.nan
df["bmfs_yr"] = df.get("bmfs (yr)", pd.Series(index=df.index,dtype=object)).map(fnum)
# event: if bmfs present but censored? GSE2603 provides bmfs time; infer event from 'bmfs event' if present
print("bmfs non-null:", df["bmfs_yr"].notna().sum())
print(df[[c for c in df.columns if 'bmf' in c.lower() or 'event' in c.lower() or 'status' in c.lower()]].head())

# ---------- GSE12276 : survival time ----------
print("="*60); print("GSE12276 (survival time months)")
expr2, meta2 = read_series_matrix(os.path.join(RAW,"GSE12276_series_matrix.txt.gz"))
df2 = project(expr2, meta2, "GPL570")
df2["time_mo"] = df2.get("survival time (months)", pd.Series(index=df2.index,dtype=object)).map(fnum)
# event status column?
print("columns:", [c for c in df2.columns])
print("time non-null:", df2["time_mo"].notna().sum())
print(df2.head())

# save
df.to_csv(os.path.join(RESULTS,"tables","survival_gse2603.csv"))
df2.to_csv(os.path.join(RESULTS,"tables","survival_gse12276.csv"))
