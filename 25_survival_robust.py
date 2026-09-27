# 25_survival_robust.py - corrected survival: per-SD HR, proliferation-adjusted, separate endpoints, random gene-set null
import os, gzip, sys
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

adapt = pd.read_csv(os.path.join(RESULTS,"tables","geneset_adaptation.txt"), header=None)[0].tolist()
comp  = pd.read_csv(os.path.join(RESULTS,"tables","geneset_competence.txt"), header=None)[0].tolist()
prolif = ["MKI67","TOP2A","PCNA","CCNB1","CDK1","MCM2","BIRC5"]

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

def build(expr, gpl, gene_lists):
    pmap = load_map(gpl)
    e = expr.copy()
    e["symbol"]=[pmap.get(i) for i in e.index]
    e = e.dropna(subset=["symbol"]).groupby("symbol").mean()
    Z = e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1).replace(0,np.nan),axis=0)
    out = {}
    for name,genes in gene_lists.items():
        g = [x for x in genes if x in Z.index]
        out[name] = Z.loc[g].mean(axis=0)
    return pd.DataFrame(out), Z.index.tolist()

def run_cohort(name, sm_file, gpl, time_field, time_unit):
    expr, meta = read_series_matrix(os.path.join(RAW,sm_file))
    sc, allgenes = build(expr, gpl, {"comp":comp,"adapt":adapt,"prolif":prolif})
    # extract time + event
    time={}; ev={}
    for s in meta.index:
        for c in meta.columns:
            v = meta.loc[s,c]
            if not isinstance(v,str): continue
            lv=v.lower()
            if lv.startswith(time_field):
                try:
                    t=float(v.split(":")[1])
                    if t>0: time[s]=t
                except: pass
            # event / vital status
            if "vital status" in lv or "event" in lv or "status" in lv:
                ev[s] = 1 if ("death" in lv or "dead" in lv or "recurrence" in lv or "relapse" in lv) else 0
    df = sc.join(pd.Series(time,name="time")).dropna()
    # if no explicit event column, assume event=1 (exploratory; flag it)
    if name=="GSE2603":
        df["event"]=1  # recurrence-free cohort, all event=1 exploratory
    else:
        evs = pd.Series(ev)
        df = df.join(evs.rename("event"))
        df["event"]=df["event"].fillna(1)
    df = df[df.time>0]
    # standardize scores by cohort SD
    for c in ["comp","adapt","prolif"]:
        df[c+"_z"] = (df[c]-df[c].mean())/df[c].std()
    print(f"\n=== {name} (n={len(df)}) ===")
    # univariate comp (per SD)
    cph= CoxPHFitter()
    cph.fit(df[["time","event","comp_z"]], "time","event")
    hr = cph.summary.loc["comp_z","exp(coef)"]; p=cph.summary.loc["comp_z","p"]
    print(f"univariate comp: HR per SD = {hr:.2f}, p={p:.4g}")
    # multivariable: + proliferation + adapt
    try:
        cph2= CoxPHFitter()
        cph2.fit(df[["time","event","comp_z","adapt_z","prolif_z"]], "time","event")
        s=cph2.summary
        print("multivariable (comp+adapt+prolif):")
        print(s[["exp(coef)","p"]].round(3).to_string())
    except Exception as ex:
        print("multi failed", ex)
    return df, allgenes

df2603,genes2603 = run_cohort("GSE2603","GSE2603_series_matrix.txt.gz","GPL96","bmfs","yr")
df12276,genes12276 = run_cohort("GSE12276","GSE12276_series_matrix.txt.gz","GPL570","survival time","mo")

# random gene-set null on GSE12276: same size as comp, recompute HR
print("\n=== random gene-set null (GSE12276), n=500 ===")
rng=np.random.default_rng(1)
obs_hr = None
nul=[]
Z = None
# rebuild raw Z for gpl570
expr, meta = read_series_matrix(os.path.join(RAW,"GSE12276_series_matrix.txt.gz"))
pmap = load_map("GPL570")
e=expr.copy(); e["symbol"]=[pmap.get(i) for i in e.index]
e=e.dropna(subset=["symbol"]).groupby("symbol").mean()
Z=e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1).replace(0,np.nan),axis=0)
common = list(set(Z.columns) & set(df12276.index))
n = len([g for g in comp if g in Z.index])
nul=[]
for i in range(500):
    pick = rng.choice(Z.index, n, replace=False)
    score = Z.loc[pick, common].mean(axis=0)
    d = pd.DataFrame({"time":df12276.loc[common,"time"],
                      "event":df12276.loc[common,"event"],
                      "s":score})
    d["s"]=(d.s-d.s.mean())/d.s.std()
    try:
        c=CoxPHFitter(); c.fit(d,"time","event")
        nul.append(c.summary.loc["s","exp(coef)"])
    except: pass
nul=np.array(nul)
# observed comp HR on GSE12276 per SD, recomputed here
dobs = df12276[["time","event","comp_z"]].dropna()
cobs=CoxPHFitter(); cobs.fit(dobs,"time","event")
obs_hr = cobs.summary.loc["comp_z","exp(coef)"]
print(f"random-gene HR distribution: median={np.median(nul):.2f} 95%CI=[{np.percentile(nul,2.5):.2f},{np.percentile(nul,97.5):.2f}]")
print(f"observed comp HR (per SD) = {obs_hr:.2f}; empirical p~{(np.sum(nul>=obs_hr)+1)/len(nul):.4f}")
pd.Series(nul).to_csv(os.path.join(RESULTS,"tables","random_geneset_hr_null.csv"),index=False)
print("saved null.")
