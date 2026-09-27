"""07_inspect_survival.py - find survival/clinical endpoints in cohorts."""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, read_series_matrix

for gse in ["GSE12276","GSE2034","GSE2603","GSE5327","GSE43837","GSE125989"]:
    expr, meta = read_series_matrix(os.path.join(RAW, f"{gse}_series_matrix.txt.gz"))
    print(f"\n===== {gse}: {expr.shape} =====")
    # print characteristics fields
    for col in meta.columns:
        if "character" in col.lower() or "title" in col.lower() or "source" in col.lower():
            vals = meta[col].dropna().unique()[:4]
            print(f"  {col}: {list(vals)}")
