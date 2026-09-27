import pandas as pd, sys
x = pd.ExcelFile(r'D:\BCBM_Project\data\raw\GSE184869_expression.xlsx')
print('sheets:', x.sheet_names, flush=True)
df = pd.read_excel(x, sheet_name=0, nrows=5)
print('cols:', list(df.columns)[:15], flush=True)
print(df.iloc[:3, :6].to_string(), flush=True)
