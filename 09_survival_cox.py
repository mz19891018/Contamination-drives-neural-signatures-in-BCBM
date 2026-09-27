"""09_survival_cox.py - Cox and KM for program scores."""
import os, sys, gzip
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

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

def scores(expr, gpl):
    pmap = load_map(gpl)
    e = expr.copy()
    e["symbol"]=[pmap.get(i) for i in e.index]
    e = e.dropna(subset=["symbol"]).groupby("symbol").mean()
    Z = e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1).replace(0,np.nan),axis=0)
    return pd.DataFrame({"comp":Z.loc[[g for g in comp if g in Z.index]].mean(axis=0),
                         "adapt":Z.loc[[g for g in adapt if g in Z.index]].mean(axis=0)})

results = {}

# ---- GSE2603: bmfs ----
expr, meta = read_series_matrix(os.path.join(RAW,"GSE2603_series_matrix.txt.gz"))
sc = scores(expr, "GPL96")
# bmfs per sample
bmfs = {}
for s in meta.index:
    for c in meta.columns:
        if "characteristics" in c.lower():
            v = meta.loc[s,c]
            if isinstance(v,str) and v.startswith("bmfs"):
                try: bmfs[s] = float(v.split(":")[1])
                except: pass
df = sc.join(pd.Series(bmfs, name="time_yr")).dropna()
# all event (recurrence cohort); exploratory
df["event"] = 1
df = df[df.time_yr > 0]
print(f"GSE2603 n={len(df)}")
cph = CoxPHFitter()
cph.fit(df[["time_yr","event","comp","adapt"]], duration_col="time_yr", event_col="event")
print(cph.summary[["coef","exp(coef)","p"]])
results["GSE2603_bmfs"] = cph.summary

# KM by median comp
med = df["comp"].median()
fig, ax = plt.subplots(1,2, figsize=(10,4))
for i,(prog,tl) in enumerate([("comp","Competence score"),("adapt","Adaptation score")]):
    km = KaplanMeierFitter()
    hi = df[df[prog] >= df[prog].median()]; lo = df[df[prog] < df[prog].median()]
    km.fit(hi.time_yr, hi.event, label="high")
    km.plot_survival_function(ax=ax[i], ci_show=False)
    km.fit(lo.time_yr, lo.event, label="low")
    km.plot_survival_function(ax=ax[i], ci_show=False)
    lr = logrank_test(hi.time_yr, lo.time_yr, hi.event, lo.event)
    ax[i].set_title(f"{tl}\nlogrank p={lr.p_value:.3f}")
    ax[i].set_xlabel("years"); ax[i].set_ylabel("BM-free survival")
plt.tight_layout()
plt.savefig(os.path.join(RESULTS,"figures","km_gse2603_bmfs.png"), dpi=120)
plt.close()

# ---- GSE12276: survival ----
expr2, meta2 = read_series_matrix(os.path.join(RAW,"GSE12276_series_matrix.txt.gz"))
sc2 = scores(expr2, "GPL570")
surv = {}
for s in meta2.index:
    for c in meta2.columns:
        if "characteristics" in c.lower():
            v = meta2.loc[s,c]
            if isinstance(v,str) and v.lower().startswith("survival time"):
                try: surv[s]=float(v.split(":")[1])
                except: pass
df2 = sc2.join(pd.Series(surv, name="time_mo")).dropna()
df2["event"]=1
df2 = df2[df2.time_mo>0]
print(f"\nGSE12276 n={len(df2)}")
cph2 = CoxPHFitter()
cph2.fit(df2[["time_mo","event","comp","adapt"]], duration_col="time_mo", event_col="event")
print(cph2.summary[["coef","exp(coef)","p"]])
results["GSE12276_survival"] = cph2.summary

import pickle
with open(os.path.join(RESULTS,"tables","cox_results.pkl"),"wb") as f:
    pickle.dump({k:v.reset_index() for k,v in results.items()}, f)
print("\nDone. KM saved.")
