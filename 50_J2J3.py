# 50_J2J3_analyze.py
import pandas as pd, numpy as np, re
D=r"D:\BCBM_Project"
for name,path in [("J2_balanced",rf"{D}\results\tables\J2_balanced_fractions.csv"),
                  ("J3_no_astro",rf"{D}\results\tables\J3_no_astrocyte_fractions.csv"),
                  ("original",rf"{D}\results\tables\bp_combined_fractions.csv")]:
    th=pd.read_csv(path,index_col=0)
    th["kind"]=["BrainMet" if re.match(r"^BM|^\d+M_RCS|^MAYO",s) else "Primary" if re.match(r"^BP|^\d+P_RCS",s) else "Other" for s in th.index]
    brain_cols=[c for c in ["Neuron","Astrocyte","Oligodendrocyte"] if c in th.columns]
    th["brain"]=th[brain_cols].sum(1)
    print(f"\n=== {name} ===")
    print(th.groupby("kind")[["Tumor","brain"]+brain_cols].mean().round(3).to_string())
