import openpyxl, re
wb=openpyxl.load_workbook(r"D:\BCBM_Project\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]
header=next(ws.iter_rows(values_only=True))
bm=[h for h in header if h and re.match(r"^BM",str(h))]
bp=[h for h in header if h and re.match(r"^BP",str(h))]
other=[h for h in header[1:] if h and not re.match(r"^B[MP]",str(h))]
print("BM:",len(bm), bm[:5])
print("BP:",len(bp), bp[:5])
print("other:",len(other), other[:10])
wb.close()
