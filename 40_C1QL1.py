# 40_C1QL1.py
import pandas as pd, numpy as np
D=r"D:\BCBM_Project"
# G1: primary expression of C1QL1 in GSE184869
import openpyxl
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]
rows=ws.iter_rows(values_only=True)
header=next(rows)
ens_col=header[0]
# find C1QL1 ensembl id
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
c1ql1_ens=pm[pm.symbol=="C1QL1"].ensembl.tolist()
print("C1QL1 ensembl ids:",c1ql1_ens)
bp_cols=[i for i,h in enumerate(header) if h and str(h).startswith("BP")]
bm_cols=[i for i,h in enumerate(header) if h and str(h).startswith("BM")]
print("BP samples:",len(bp_cols)," BM samples:",len(bm_cols))
for r in rows:
    if r[0] in c1ql1_ens:
        bp_vals=[r[i] for i in bp_cols if r[i] is not None]
        bm_vals=[r[i] for i in bm_cols if r[i] is not None]
        bp_det=np.mean(np.array(bp_vals)>1)
        print(f"C1QL1 {r[0]}: BP mean={np.mean(bp_vals):.2f} detect(>1CPM)={bp_det:.1%} n={len(bp_vals)}")
        print(f"             BM mean={np.mean(bm_vals):.2f} n={len(bm_vals)}")
wb.close()

# G2: single cell C1QL1 in tumor vs neurons - read from seurat via R export? 
# Use tumor_authenticity_scores if available
try:
    ta=pd.read_csv(rf"{D}\results\tables\tumor_authenticity_scores.csv")
    print("\ntumor authenticity cols:",list(ta.columns)[:10])
    if "C1QL1" in ta.columns:
        print("C1QL1 tumor mean:",ta["C1QL1"].mean())
except Exception as e:
    print("no tumor authenticity:",e)
