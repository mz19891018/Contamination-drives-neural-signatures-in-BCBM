# 59_survival_fixed.py - redo survival with proper event indicators
import gzip, os, sys
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, read_series_matrix

adapt = [g.strip() for g in open(os.path.join(RESULTS,"tables","geneset_adaptation.txt"),encoding="utf-8") if g.strip()]
comp  = [g.strip() for g in open(os.path.join(RESULTS,"tables","geneset_competence.txt"),encoding="utf-8") if g.strip()]

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

def parse_characteristics(series_path):
    """Parse all !Sample_characteristics_ch1 lines into a dict of sample->field->value."""
    sample_ids=[]
    char_lines=[]
    with gzip.open(series_path,"rt",errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                sample_ids=[v.strip().strip('"') for v in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                char_lines.append([v.strip().strip('"') for v in line.split("\t")[1:]])
    # build dict
    data={sid:{} for sid in sample_ids}
    for vals in char_lines:
        for i,sid in enumerate(sample_ids):
            if i < len(vals) and ":" in vals[i]:
                k,v=vals[i].split(":",1)
                data[sid][k.strip()]=v.strip()
    return data

# ===== GSE2603 =====
print("=== GSE2603 (BMFS) ===")
expr, meta = read_series_matrix(os.path.join(RAW,"GSE2603_series_matrix.txt.gz"))
sc = scores(expr, "GPL96")
ch = parse_characteristics(os.path.join(RAW,"GSE2603_series_matrix.txt.gz"))
rows=[]
for sid in sc.index:
    if sid not in ch: continue
    c=ch[sid]
    tissue=c.get("tissue type","")
    is_cellline = any(k in tissue.lower() for k in ["cell line","mda","subpopulation"])
    bmfs=c.get("bmfs (yr)","--")
    bmevt=c.get("bm event","--")
    if bmfs=="--" or bmevt=="--": continue
    try:
        t=float(bmfs); e=int(bmevt)
    except: continue
    if t<=0: continue
    rows.append({"sample":sid,"time_yr":t,"event":e,"comp":sc.loc[sid,"comp"],"adapt":sc.loc[sid,"adapt"],"cell_line":is_cellline})
df=pd.DataFrame(rows)
df=df[~df.cell_line].drop(columns=["cell_line"])
print(f"n={len(df)}, events={df.event.sum()}, censored={(df.event==0).sum()}")
print(f"time: median={df.time_yr.median():.2f}, range=[{df.time_yr.min():.2f},{df.time_yr.max():.2f}]")

# Cox per SD
df["comp_sd"]=(df.comp-df.comp.mean())/df.comp.std()
cph=CoxPHFitter()
cph.fit(df[["time_yr","event","comp_sd"]],duration_col="time_yr",event_col="event")
print("\nUnivariable comp (per SD):")
print(cph.summary[["coef","exp(coef)","p"]].round(4))
# proliferation adjustment
df["prolif"]=df.comp_sd  # placeholder - need real proliferation
# Actually compute proliferation from expression
pmap=load_map("GPL96")
e2=expr.copy(); e2["symbol"]=[pmap.get(i) for i in e2.index]
e2=e2.dropna(subset=["symbol"]).groupby("symbol").mean()
Z2=e2.sub(e2.mean(axis=1),axis=0).div(e2.std(axis=1).replace(0,np.nan),axis=0)
prolif_genes=["MKI67","TOP2A","PCNA","CCNB1"]
prolif=Z2.loc[[g for g in prolif_genes if g in Z2.index]].mean(0)
df["prolif_raw"]=df["sample"].map(prolif)
df["prolif_sd"]=(df["prolif_raw"]-df["prolif_raw"].mean())/df["prolif_raw"].std()
cph2=CoxPHFitter()
cph2.fit(df[["time_yr","event","comp_sd","prolif_sd"]],duration_col="time_yr",event_col="event")
print("\nMultivariable comp + prolif (per SD):")
print(cph2.summary[["coef","exp(coef)","p"]].round(4))
df.to_csv(os.path.join(RESULTS,"tables","survival_gse2603_fixed.csv"),index=False)

# ===== GSE12276 =====
print("\n=== GSE12276 (OS) ===")
expr2, meta2 = read_series_matrix(os.path.join(RAW,"GSE12276_series_matrix.txt.gz"))
sc2 = scores(expr2, "GPL570")
ch2 = parse_characteristics(os.path.join(RAW,"GSE12276_series_matrix.txt.gz"))
# check for event fields
all_keys=set()
for sid,c in ch2.items(): all_keys.update(c.keys())
print("Available fields:", sorted(all_keys))
rows2=[]
for sid in sc2.index:
    if sid not in ch2: continue
    c=ch2[sid]
    st=c.get("Survival time (months)",c.get("survival time (months)","--"))
    if st=="--": continue
    try: t=float(st)
    except: continue
    if t<=0: continue
    # No explicit death event in GEO metadata
    rows2.append({"sample":sid,"time_mo":t,"event":np.nan,"comp":sc2.loc[sid,"comp"],"adapt":sc2.loc[sid,"adapt"]})
df2=pd.DataFrame(rows2)
print(f"n={len(df2)} with survival time, but NO event indicator available in GEO metadata")
print("NOTE: GSE12276 OS analysis cannot be properly redone without death event data from original publication.")
print("Previous event=1 assumption is invalid. Reporting as data limitation.")
df2.to_csv(os.path.join(RESULTS,"tables","survival_gse12276_noevent.csv"),index=False)

print("\nDone. GSE2603 fixed with proper events; GSE12276 event unavailable.")
