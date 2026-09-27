"""04_inspect_multiorgan.py - inspect GSE14017/18 organ labels."""
import os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, read_series_matrix

for gse in ["GSE14017", "GSE14018"]:
    expr, meta = read_series_matrix(os.path.join(RAW, f"{gse}_series_matrix.txt.gz"))
    print(f"\n===== {gse}: expr {expr.shape} =====")
    for col in ["!Sample_title", "!Sample_source_name_ch1", "!Sample_characteristics_ch1"]:
        if col in meta.columns:
            print(f"--- {col} ---")
            for v in meta[col].head(15):
                print("  ", v)
