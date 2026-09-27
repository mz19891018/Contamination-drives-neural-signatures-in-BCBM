# 63_survival_comprehensive.py - full survival diagnostics
import gzip, os, sys, numpy as np, pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import spearmanr, mannwhitneyu
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

adapt=[g.strip() for g in open(os.path.join(RESULTS,"tables","geneset_adaptation.txt"),encoding="utf-8") if g.strip()]
comp=[g.strip() for g in open(os.path.join(RESULTS,"tables","geneset_competence.txt"),encoding="utf-8") if g.strip()]

def load_map(gpl):
    rows={}
    with gzip.open(os.path.join(RAW,f"{gpl}.annot.gz"),"rt",errors="replace",encoding="utf-8") as f:
        for line in f:
            if line.startswith(("#","!","^")): continue
            p=line.rstrip("\n").split("\t")
            if len(p)>2: rows[p[0]]=p[2].split("///")[0].strip()
    return pd.Series(rows)

def scores(expr,gpl):
    pmap=load_map(gpl)
    e=expr.copy(); e["symbol"]=[pmap.get(i) for i in e.index]
    e=e.dropna(subset=["symbol"]).groupby("symbol").mean()
    Z=e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1).replace(0,np.nan),axis=0)
    return pd.DataFrame({"comp":Z.loc[[g for g in comp if g in Z.index]].mean(0),
                         "adapt":Z.loc[[g for g in adapt if g in Z.index]].mean(0)})

def parse_chars(path):
    sids=[]; cls=[]
    with gzip.open(path,"rt",errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                sids=[v.strip().strip('"') for v in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                cls.append([v.strip().strip('"') for v in line.split("\t")[1:]])
    data={s:{} for s in sids}
    for vals in cls:
        for i,s in enumerate(sids):
            if i<len(vals) and ":" in vals[i]:
                k,v=vals[i].split(":",1); data[s][k.strip()]=v.strip()
    return data

# ========== GSE2603: technical confounding + dfbeta + bootstrap ==========
print("="*60)
print("GSE2603: TECHNICAL CONFOUNDING + DFBETA + BOOTSTRAP")
print("="*60)
expr,meta=read_series_matrix(os.path.join(RAW,"GSE2603_series_matrix.txt.gz"))
sc=scores(expr,"GPL96")
ch=parse_chars(os.path.join(RAW,"GSE2603_series_matrix.txt.gz"))
rows=[]
for sid in sc.index:
    if sid not in ch: continue
    c=ch[sid]; tissue=c.get("tissue type","")
    is_cl=any(k in tissue.lower() for k in ["cell line","mda","subpopulation"])
    bmfs=c.get("bmfs (yr)","--"); bmevt=c.get("bm event","--")
    er=c.get("path er status","--"); her2=c.get("her2 status","--")
    if bmfs=="--" or bmevt=="--" or is_cl: continue
    try: t=float(bmfs); e=int(bmevt)
    except: continue
    if t<=0: continue
    rows.append({"sample":sid,"time":t,"event":e,"comp":sc.loc[sid,"comp"],
                 "er":er,"her2":her2})
df=pd.DataFrame(rows)
df["comp_sd"]=(df.comp-df.comp.mean())/df.comp.std()
print(f"n={len(df)}, events={df.event.sum()}, censored={(df.event==0).sum()}")
print(f"Event comp_sd median: {df[df.event==1].comp_sd.median():.3f}")
print(f"Censored comp_sd median: {df[df.event==0].comp_sd.median():.3f}")

# ER/HER2 confounding
print("\n--- ER status ---")
print(df.groupby("er")["comp_sd"].agg(["mean","std","count"]))
print("\n--- HER2 status ---")
print(df.groupby("her2")["comp_sd"].agg(["mean","std","count"]))

# Cox univariable
cph=CoxPHFitter()
cph.fit(df[["time","event","comp_sd"]],duration_col="time",event_col="event")
hr=cph.summary.loc["comp_sd","exp(coef)"]; p=cph.summary.loc["comp_sd","p"]
print(f"\nUnivariable: HR={hr:.3f}, p={p:.4f}")

# dfbeta influence
betas=[]
for i in range(len(df)):
    d2=df.drop(i)
    try:
        c=CoxPHFitter(); c.fit(d2[["time","event","comp_sd"]],duration_col="time",event_col="event")
        betas.append(c.summary.loc["comp_sd","coef"])
    except: betas.append(np.nan)
dfbeta=np.array(betas)-cph.summary.loc["comp_sd","coef"]
print(f"\n--- DFBETA (top 5 influential) ---")
inf=pd.DataFrame({"sample":df["sample"].values,"time":df["time"].values,"event":df["event"].values,
                   "comp_sd":df["comp_sd"].values,"dfbeta":dfbeta})
print(inf.reindex(inf.dfbeta.abs().sort_values(ascending=False).index).head(5).to_string())

# bootstrap CI
np.random.seed(42)
boot_hrs=[]
for _ in range(1000):
    idx=np.random.choice(len(df),len(df),replace=True)
    d2=df.iloc[idx]
    if d2.event.sum()<3: continue
    try:
        c=CoxPHFitter(); c.fit(d2[["time","event","comp_sd"]],duration_col="time",event_col="event")
        boot_hrs.append(c.summary.loc["comp_sd","exp(coef)"])
    except: pass
print(f"\nBootstrap HR: median={np.median(boot_hrs):.3f}, 95% CI=[{np.percentile(boot_hrs,2.5):.3f},{np.percentile(boot_hrs,97.5):.3f}] (n={len(boot_hrs)})")

# 500 random gene sets with CORRECT event coding
print("\n--- 500 random gene sets (correct event coding) ---")
pmap=load_map("GPL96")
e2=expr.copy(); e2["symbol"]=[pmap.get(i) for i in e2.index]
e2=e2.dropna(subset=["symbol"]).groupby("symbol").mean()
Z2=e2.sub(e2.mean(axis=1),axis=0).div(e2.std(axis=1).replace(0,np.nan),axis=0)
all_genes=[g for g in Z2.index if g not in comp]
np.random.seed(123)
null_hrs=[]; null_ps=[]
for i in range(500):
    gs=np.random.choice(all_genes,len(comp),replace=False)
    gs_in=[g for g in gs if g in Z2.index]
    s=Z2.loc[gs_in].mean(0)
    d2=df.copy(); d2["rand_raw"]=d2["sample"].map(s)
    d2["rand_sd"]=(d2["rand_raw"]-d2["rand_raw"].mean())/d2["rand_raw"].std()
    try:
        c=CoxPHFitter(); c.fit(d2[["time","event","rand_sd"]],duration_col="time",event_col="event")
        null_hrs.append(c.summary.loc["rand_sd","exp(coef)"])
        null_ps.append(c.summary.loc["rand_sd","p"])
    except: pass
print(f"Null HR: median={np.median(null_hrs):.3f}, IQR=[{np.percentile(null_hrs,25):.3f},{np.percentile(null_hrs,75):.3f}]")
print(f"Null HR < 1: {sum(1 for h in null_hrs if h<1)}/{len(null_hrs)} ({100*sum(1 for h in null_hrs if h<1)/len(null_hrs):.0f}%)")
print(f"Observed HR={hr:.3f}, empirical p (two-sided)={2*min(sum(1 for h in null_hrs if h<=hr),sum(1 for h in null_hrs if h>=hr))/len(null_hrs):.3f}")

# ========== GSE12276: BMFS from site of relapse ==========
print("\n"+"="*60)
print("GSE12276: BMFS FROM SITE OF RELAPSE")
print("="*60)
expr2,meta2=read_series_matrix(os.path.join(RAW,"GSE12276_series_matrix.txt.gz"))
sc2=scores(expr2,"GPL570")
ch2=parse_chars(os.path.join(RAW,"GSE12276_series_matrix.txt.gz"))
rows2=[]
for sid in sc2.index:
    if sid not in ch2: continue
    c=ch2[sid]
    site=c.get("site of relapse (brain or other)","")
    st=c.get("survival time (months)",c.get("Survival time (months)",""))
    if st=="" or st=="--": continue
    try: t=float(st)
    except: continue
    if t<=0: continue
    # event: brain or brain+other = 1; other/local = censored at relapse; '' = censored
    if site in ["brain","brain+other"]:
        event=1
    elif site in ["other","local"]:
        event=0  # competing event: relapsed elsewhere
    else:
        event=0  # no relapse
    rows2.append({"sample":sid,"time_mo":t,"event":event,"site":site,"comp":sc2.loc[sid,"comp"]})
df2=pd.DataFrame(rows2)
df2["comp_sd"]=(df2.comp-df2.comp.mean())/df2.comp.std()
print(f"n={len(df2)}, brain events={df2.event.sum()}, censored={(df2.event==0).sum()}")
print(f"Site distribution: {df2.site.value_counts().to_dict()}")
print(f"Event time range: {df2[df2.event==1].time_mo.min():.1f}-{df2[df2.event==1].time_mo.max():.1f} months")
if df2.event.sum()>=3:
    cph2=CoxPHFitter()
    cph2.fit(df2[["time_mo","event","comp_sd"]],duration_col="time_mo",event_col="event")
    hr2=cph2.summary.loc["comp_sd","exp(coef)"]; p2=cph2.summary.loc["comp_sd","p"]
    print(f"BMFS Cox: HR={hr2:.3f}/SD, p={p2:.4f}")
else:
    print("Too few events for Cox")
df2.to_csv(os.path.join(RESULTS,"tables","survival_gse12276_bmfs.csv"),index=False)

# Save GSE2603 diagnostics
inf.to_csv(os.path.join(RESULTS,"tables","survival_gse2603_dfbeta.csv"),index=False)
pd.DataFrame({"null_hr":null_hrs,"null_p":null_ps}).to_csv(os.path.join(RESULTS,"tables","survival_null_correct.csv"),index=False)
print("\nSaved all survival diagnostics.")
