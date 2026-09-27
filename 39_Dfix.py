# 39_Dfix.py - only within BM samples
import pandas as pd, numpy as np, gzip, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"
BRAIN=["ATP1A2","TUBB4A","MAG","MOG","MBP","PLP1","APLP1","C1QL1","PAX6","NEUROD1",
       "RBFOX3","SYT1","SNAP25","GFAP","AQP4","S100B","SLC1A3","MAL","MOBP","CNP",
       "CHGB","SCG2","STMN2","SYN1","UCHL1","GRIN1","SYP","GAD1","HOMER1","NRGN"]
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]

def build_map():
    m={}
    with gzip.open(rf"{D}\data\raw\GPL570.annot.gz",'rt',errors='ignore') as f:
        in_tab=False; header=None
        for line in f:
            if line.startswith("!platform_table_begin"): in_tab=True; continue
            if line.startswith("!platform_table_end"): break
            if not in_tab: continue
            p=line.rstrip("\n").split("\t")
            if header is None: header=p; i_id=header.index("ID"); i_sym=header.index("Gene symbol"); continue
            if len(p)<=i_sym: continue
            m[p[i_id]]=p[i_sym].split(" /// ")[0]
    return m
PM=build_map()

def parse(path):
    rows={}; started=False; titles={}
    with gzip.open(path,'rt',encoding='utf-8',errors='ignore') as f:
        for line in f:
            if line.startswith('"ID_REF"') or line.startswith("ID_REF"): started=True
            if not started:
                if line.startswith("!Sample_title"):
                    p=[x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
                    titles=p[1:]
                continue
            p=[x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if not hasattr(parse,"ncol"): parse.ncol=len(p)-1
            p=p[:parse.ncol+1]+[""]*(parse.ncol+1-len(p))
            rows[p[0]]=p[1:]
    df=pd.DataFrame(rows).T
    df.columns=df.iloc[0]; df=df.iloc[1:]
    df=df.apply(pd.to_numeric,errors='coerce')
    df["sym"]=[PM.get(i,i) for i in df.index]
    df=df[df.sym!=""].groupby("sym").mean()
    return df, titles

def boot_rho(x,y,n=1000):
    rng=np.random.default_rng(0); rs=[]
    for _ in range(n):
        idx=rng.choice(len(x),len(x),replace=True)
        rs.append(spearmanr(x[idx],y[idx])[0])
    return np.percentile(rs,2.5),np.percentile(rs,97.5)

res=[]
for name,path in [("GSE125989",rf"{D}\data\raw\GSE125989_series_matrix.txt.gz"),
                  ("GSE14018",rf"{D}\data\raw\GSE14018_series_matrix.txt.gz"),
                  ("GSE14017",rf"{D}\data\raw\GSE14017_series_matrix.txt.gz")]:
    df,titles=parse(path)
    # identify BM samples by title
    bm_mask=[bool(re.search(r"brain|脑",t,re.I)) for t in titles]
    bm_cols=[c for c,m in zip(df.columns,bm_mask) if m]
    print(f"{name}: total {len(titles)}, BM {len(bm_cols)}")
    if len(bm_cols)<5:
        res.append({"dataset":name,"n_BM":len(bm_cols),"rho":np.nan,"ci_lo":np.nan,"ci_hi":np.nan,"note":"n<5"}); continue
    sub=df[bm_cols]
    bh=[g for g in BRAIN if g in sub.index]; ah=[g for g in adapt if g in sub.index]
    bs=sub.loc[bh].mean(0); ascore=sub.loc[ah].mean(0)
    r,p=spearmanr(bs,ascore)
    lo,hi=boot_rho(bs.values,ascore.values)
    res.append({"dataset":name,"n_BM":len(bm_cols),"rho":round(r,3),"ci_lo":round(lo,3),"ci_hi":round(hi,3),"p":f"{p:.2e}"})
    print(f"  rho={r:.3f} CI[{lo:.3f},{hi:.3f}] n={len(bm_cols)}")
pd.DataFrame(res).to_csv(rf"{D}\results\tables\Dfix_bm_only.csv",index=False)
print("saved")
