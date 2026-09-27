"""
01_explore.py - inspect GSE184869 (45 paired primary-BM) and GSE38057 metadata.
"""
import pandas as pd, numpy as np, gzip, re, os

RAW = r"D:\BCBM_Project\data\raw"

# ---------- GSE184869 : Excel expression matrix ----------
xlsx = os.path.join(RAW, "GSE184869_expression.xlsx")
xl = pd.ExcelFile(xlsx)
print("=== GSE184869 sheets ===", xl.sheet_names)
for sh in xl.sheet_names:
    df = xl.parse(sh, nrows=5)
    print(f"\n--- sheet: {sh}  shape(head)={df.shape}")
    print(df.iloc[:5, :8].to_string())
