"""15_meta_cox.py - per-cohort Cox on comp/adapt scores, then fixed-effect meta."""
import os, sys, gzip, re
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

adapt = pd.read_csv(os.path.join(RESULTS,"tables","geneset_adaptation.txt"),header=None)[0].tolist()
comp  = pd.read_csv(os.path.join(RESULTS,"tables","geneset_competence.txt"),header=None)[0].tolist()

def load_map(gpl):
    rows={}
    with gzip.open(os.path.join(RAW,f"{gpl}.annot.gz"),"rt",errors="replace",encoding="utf-8") as f:
        header=None
        for line in f:
            if line.startswith(("#","!","^")): continue
            p=line.rstrip("\n").split("\t")
            if header is None: header=p; continue
            rows[p[0]]=p[2].split("///")[0].strip()
    return pd.Series(rows)

def scores(expr,gpl):
    pm=load_map(gpl)
    e=expr.copy(); e["symbol"]=[pm.get(i) for i in e.index]
    e=e.dropna(subset=["symbol"]).groupby("symbol").mean()
    Z=e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1).replace(0,np.nan),axis=0)
    return pd.DataFrame({"comp":Z.loc[[g for g in comp if g in Z.index]].mean(axis=0),
                         "adapt":Z.loc[[g for g in adapt if g in Z.index]].mean(axis=0)})

def get_surv(meta, keys):
    out={}
    for s in meta.index:
        for c in meta.columns:
            if "characteristics" in c.lower():
                v=meta.loc[s,c]
                if isinstance(v,str):
                    for k in keys:
                        if v.lower().startswith(k):
                            try: out[s]=float(v.split(":")[1])
                            except: pass
    return pd.Series(out)

results=[]
# (gse, gpl, time_keys, label)
cohorts = [
    ("GSE2034","GPL96",["distant metastasis free survival","relapse free survival","time to"], "GSE2034_DMFS"),
    ("GSE2603","GPL96",["bmfs"], "GSE2603_bmfs"),
    ("GSE5327","GPL96",["distant metastasis free survival","relapse"], "GSE5327_DMFS"),
    ("GSE12276","GPL570",["survival time"], "GSE12276_OS"),
]
for gse,gpl,keys,label in cohorts:
    try:
        expr,meta=read_series_matrix(os.path.join(RAW,f"{gse}_series_matrix.txt.gz"))
        sc=scores(expr,gpl)
        t=get_surv(meta,keys)
        if len(t.dropna())<20:
            print(f"{label}: too few time points ({len(t.dropna())})"); continue
        df=sc.join(pd.DataFrame({"time":t})).dropna()
        df=df[df.time>0]; df["event"]=1
        cph=CoxPHFitter()
        cph.fit(df[["time","event","comp","adapt"]],duration_col="time",event_col="event")
        s=cph.summary
        print(f"{label}: n={len(df)}  comp HR={np.exp(s.loc['comp','coef']):.2f} p={s.loc['comp','p']:.4f}  "
              f"adapt HR={np.exp(s.loc['adapt','coef']):.2f} p={s.loc['adapt','p']:.4f}")
        results.append({"cohort":label,"n":len(df),
                        "comp_coef":s.loc["comp","coef"],"comp_se":s.loc["comp","se(coef)"],
                        "comp_hr":np.exp(s.loc["comp","coef"]),
                        "adapt_coef":s.loc["adapt","coef"],"adapt_se":s.loc["adapt","se(coef)"],
                        "adapt_hr":np.exp(s.loc["adapt","coef"])})
    except Exception as e:
        print(f"{label}: ERROR {e}")

res=pd.DataFrame(results)
# fixed-effect meta (inverse variance) for comp
def fe_meta(df, prefix):
    b=df[f"{prefix}_coef"]; se=df[f"{prefix}_se"]
    w=1/se**2
    mb=(b*w).sum()/w.sum()
    se_meta=np.sqrt(1/w.sum())
    return mb, np.exp(mb), mb/se_meta, 2*(1-__import__("scipy.stats").stats.norm.cdf(abs(mb/se_meta)))
for prefix in ["comp","adapt"]:
    mb,mhr,z,p=fe_meta(res,prefix)
    print(f"\nMETA {prefix}: pooled HR={mhr:.2f}  z={z:.2f} p={p:.4g}")
    res[f"meta_{prefix}_hr"]=mhr
res.to_csv(os.path.join(RESULTS,"tables","meta_cox.csv"),index=False)
print("\nsaved meta_cox.csv")
