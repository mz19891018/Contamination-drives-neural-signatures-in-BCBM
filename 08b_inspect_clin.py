import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, read_series_matrix

for gse,gpl in [("GSE2603","GPL96"),("GSE12276","GPL570")]:
    expr, meta = read_series_matrix(os.path.join(RAW,f"{gse}_series_matrix.txt.gz"))
    print(f"\n=== {gse} all characteristic values (first 3 samples) ===")
    for col in meta.columns:
        if "characteristics" in col.lower():
            for s in meta.index[:3]:
                print(f"  {s}: {meta.loc[s,col]}")
