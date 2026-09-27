"""01b_explore.py - full column list of GSE184869."""
import pandas as pd, os
RAW = r"D:\BCBM_Project\data\raw"
df = pd.read_excel(os.path.join(RAW, "GSE184869_expression.xlsx"))
print("shape:", df.shape)
print("columns:")
for c in df.columns:
    print("  ", repr(c))
