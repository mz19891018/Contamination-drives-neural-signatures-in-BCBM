"""14_subtype_sensitivity.py
Subtype-stratified paired analysis.
Assign HER2+/Luminal/TNBC to each patient using marker genes, then test whether
competence/adaptation programs differ by subtype.
"""
import os, sys, re
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, RESULTS, ensembl_to_symbol

df = pd.read_excel(os.path.join(RAW,"GSE184869_expression.xlsx")).rename(
    columns={"ensembl_gene_id":"ensembl"}).set_index("ensembl")
sym = ensembl_to_symbol(list(df.index))
df["symbol"] = pd.Series(sym)
# aggregate by symbol
d = df.dropna(subset=["symbol"]).groupby("symbol").mean()

def classify(col):
    c=col.strip()
    m=re.match(r"^MAYO_B([MP])_(\d+)$",c)
    if m: return f"mayo_{m.group(2)}","BM" if m.group(1)=="M" else "BP"
    m=re.match(r"^B([MP])([\d\-]+)$",c)
    if m: return f"rm_{m.group(2)}","BM" if m.group(1)=="M" else "BP"
    m=re.match(r"^(\d+)M_RCS$",c);
    if m: return f"rcs_{m.group(1)}","BM"
    m=re.match(r"^(\d+)P_RCS$",c)
    if m: return f"rcs_{m.group(1)}","BP"
    return None,None

pm={}
for c in d.columns:
    p,s=classify(c)
    if p: pm.setdefault(p,{})[s]=c
complete={p:v for p,v in pm.items() if "BP" in v and "BM" in v}

# subtype markers in PRIMARY tumors (BP)
markers = {"ESR1":"Luminal","PGR":"Luminal","FOXA1":"Luminal",
           "ERBB2":"HER2","GRB7":"HER2",
           "KRT5":"TNBC","KRT14":"TNBC","EGFR":"TNBC"}
# score per primary sample
def score_genes(genes):
    g=[x for x in genes if x in d.index]
    return d.loc[g].mean(axis=0) if g else pd.Series(0,index=d.columns)

lum = score_genes(["ESR1","PGR","FOXA1"])
her2 = score_genes(["ERBB2","GRB7"])
tnbc = score_genes(["KRT5","KRT14","EGFR"])

rows=[]
for p,v in complete.items():
    bp=v["BP"]
    L,H,T = lum[bp], her2[bp], tnbc[bp]
    if L > max(H,T)+0.5: sub="Luminal"
    elif H > max(L,T): sub="HER2+"
    elif T > max(L,H): sub="TNBC"
    else: sub="Unclassified"
    rows.append({"patient":p,"BP":bp,"BM":v["BM"],"subtype":sub,
                 "ESR1":L,"ERBB2":H,"basal":T})
sub = pd.DataFrame(rows).set_index("patient")
print(sub.subtype.value_counts())

# paired delta by subtype
adapt = pd.read_csv(os.path.join(RESULTS,"tables","geneset_adaptation.txt"),header=None)[0].tolist()
comp  = pd.read_csv(os.path.join(RESULTS,"tables","geneset_competence.txt"),header=None)[0].tolist()

def prog_score(col, genes):
    g=[x for x in genes if x in d.index]
    return d.loc[g,col].mean()

sub["BP_adapt"]=sub.BP.map(lambda c: prog_score(c,adapt))
sub["BM_adapt"]=sub.BM.map(lambda c: prog_score(c,adapt))
sub["adapt_delta"]=sub.BM_adapt-sub.BP_adapt
sub["BP_comp"]=sub.BP.map(lambda c: prog_score(c,comp))
sub["BM_comp"]=sub.BM.map(lambda c: prog_score(c,comp))
sub["comp_delta"]=sub.BM_comp-sub.BP_comp

print("\n=== program delta (BM-BP) by subtype ===")
print(sub.groupby("subtype")[["adapt_delta","comp_delta"]].agg(["mean","count"]).round(3))

# test adapt delta across subtypes
groups=[sub[sub.subtype==s].adapt_delta.values for s in ["Luminal","HER2+","TNBC"] if (sub.subtype==s).sum()>2]
if len(groups)>=2:
    f,p=stats.f_oneway(*groups)
    print(f"\nANOVA adapt_delta across subtypes: F={f:.2f} p={p:.4f}")
sub.to_csv(os.path.join(RESULTS,"tables","subtype_assignment.csv"))
print("\nsaved.")
