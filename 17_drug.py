"""17_drug_prioritization.py
Map adaptation program master regulators to druggable targets.
Uses top adaptation genes; annotates TF status and known BBB-penetrant pharmacology.
"""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RESULTS

prog = pd.read_csv(os.path.join(RESULTS,"tables","program_decomposition.csv"))
adapt = prog[prog.program=="Adaptation_de_novo"].sort_values("mean_logFC_BM_vs_BP",ascending=False)

# Curated annotation of top adaptation genes: function + druggability + BBB-relevance
annot = {
 "ATP1A2":("Na+/K+ ATPase alpha-2 (neuronal)","digoxin/ouabain-like inhibitors","BBB: CNS target","direct"),
 "ATP1A3":("Na+/K+ ATPase alpha-3 (neuronal)","digoxin; CNS-penetrant","BBB: CNS target","direct"),
 "C1QL1":("complement C1q-like, synaptic","synaptic; no direct drug","ligand","indirect"),
 "TUBB4A":("neuronal tubulin beta-4A","tubulin-targeting agents (taxanes cross BBB poorly)","cytoskeleton","direct"),
 "PAX6":("neural TF","TF - undruggable directly","master regulator","TF"),
 "ZIC2":("neural TF","undruggable","master regulator","TF"),
 "CHGB":("chromogranin B, neuroendocrine secretion","NE marker","secreted marker","biomarker"),
 "MAG":("myelin-associated glycoprotein","oligodendrocyte marker","ligand","indirect"),
 "APLP1":("amyloid precursor-like","no direct drug","surface","indirect"),
 "SLC1A3":("EAAT1/GLT-1 glutamate transporter","glutamate transporter modulators","BBB: CNS target","direct"),
 "ABCG2":("BCRP/ABCG2 efflux pump (BBB)","BCRP inhibitors (elacridar)","BBB efflux - modulates drug delivery","direct"),
 "PTPRZ1":("receptor tyrosine phosphatase Z (neural)","PTPRZ1 ligand pleiotrophin; emerging","receptor","direct"),
 "NRCAM":("neuronal cell adhesion","adhesion","surface","indirect"),
 "NCAN":("neurocan (CSPG)","CSPG","ECM","indirect"),
}

rows=[]
for _,r in adapt.head(40).iterrows():
    s=r.symbol
    if s in annot:
        func,drug,bbb,typ = annot[s]
    else:
        func,drug,bbb,typ="","","",""
    rows.append({"gene":s,"logFC_BMvsBP":round(r.mean_logFC_BM_vs_BP,2),
                 "mean_BP":round(r.mean_BP,2),"mean_BM":round(r.mean_BM,2),
                 "function":func,"druggability":drug,"BBB_relevance":bbb,"type":typ})
out=pd.DataFrame(rows)
print(out.to_string(index=False))
out.to_csv(os.path.join(RESULTS,"tables","drug_prioritization.csv"),index=False)

# Master regulators (TFs) in adaptation
tfs = out[out.type=="TF"]
print("\n=== Master regulator TFs in adaptation program ===")
print(tfs[["gene","function"]].to_string(index=False))
print("\n=== Direct druggable CNS targets ===")
print(out[out.type=="direct"][["gene","druggability"]].to_string(index=False))
