# 31_multi_dataset.py - contamination proxy across public BCBM microarray datasets
import pandas as pd, numpy as np, gzip, re
from scipy.stats import spearmanr
D=r"D:\BCBM_Project"

# brain marker set (contamination proxy)
BRAIN = ["ATP1A2","TUBB4A","MAG","MOG","MBP","PLP1","APLP1","C1QL1","PAX6","NEUROD1",
         "RBFOX3","SYT1","SNAP25","GFAP","AQP4","S100B","SLC1A3","MAL","MOBP","CNP",
         "CHGB","SCG2","STMN2","SYN1","UCHL1","GRIN1","SYP","GAD1","HOMER1","NRGN"]
# adaptation signature genes (from our set)
adapt = [g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]

def build_probe_map():
    m={}
    with gzip.open(rf"{D}\data\raw\GPL570.annot.gz",'rt',errors='ignore') as f:
        in_tab=False
        header=None
        for line in f:
            if line.startswith("!platform_table_begin"): in_tab=True; continue
            if line.startswith("!platform_table_end"): break
            if not in_tab: continue
            p=line.rstrip("\n").split("\t")
            if header is None:
                header=p; i_id=header.index("ID"); i_sym=header.index("Gene symbol"); continue
            if len(p)<=i_sym: continue
            sym=p[i_sym].split(" /// ")[0]
            m[p[i_id]]=sym
    return m
PROBE2SYM=build_probe_map()
print("probe map:",len(PROBE2SYM))

def read_gse_matrix(path):
    rows={}
    started=False
    with gzip.open(path,'rt',encoding='utf-8',errors='ignore') as f:
        for line in f:
            if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                started=True
            if not started: continue
            parts=[p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if not hasattr(read_gse_matrix,"ncol"):
                read_gse_matrix.ncol=len(parts)-1
            parts=parts[:read_gse_matrix.ncol+1]+[""]*(read_gse_matrix.ncol+1-len(parts))
            rows[parts[0]]=parts[1:]
    df=pd.DataFrame(rows).T
    df.columns=df.iloc[0]; df=df.iloc[1:]
    df=df.apply(pd.to_numeric,errors='coerce')
    df["sym"]=[PROBE2SYM.get(i,i) for i in df.index]
    df=df[df.sym!=""]
    df=df[~df.sym.str.startswith("---")]
    df=df.groupby("sym").mean()
    return df

datasets = {
 "GSE125989": rf"{D}\data\raw\GSE125989_series_matrix.txt.gz",
 "GSE43837":  rf"{D}\data\raw\GSE43837_series_matrix.txt.gz",
 "GSE14017":  rf"{D}\data\raw\GSE14017_series_matrix.txt.gz",
 "GSE14018":  rf"{D}\data\raw\GSE14018_series_matrix.txt.gz",
}

res=[]
for name,path in datasets.items():
    try:
        df=read_gse_matrix(path)
    except Exception as e:
        print(name,"FAIL",e); continue
    genes=list(df.index)
    brain_have=[g for g in BRAIN if g in genes]
    adapt_have=[g for g in adapt if g in genes]
    if len(brain_have)<5 or len(adapt_have)<5:
        print(name,"too few genes brain",len(brain_have),"adapt",len(adapt_have)); continue
    brain_score=df.loc[brain_have].mean(0)
    adapt_score=df.loc[adapt_have].mean(0)
    r,p=spearmanr(brain_score,adapt_score)
    res.append({"dataset":name,"n_samples":df.shape[1],"n_brain_markers":len(brain_have),
                "n_adapt_genes":len(adapt_have),"rho":round(r,3),"p":f"{p:.2e}"})
    print(f"{name}: n={df.shape[1]} rho={r:.3f} p={p:.2e}")

out=pd.DataFrame(res)
out.to_csv(rf"{D}\results\tables\multi_dataset_contamination.csv",index=False)
print("\nsaved")
