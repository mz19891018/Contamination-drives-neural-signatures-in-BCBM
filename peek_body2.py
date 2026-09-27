import gzip
with gzip.open(r"D:\BCBM_Project\data\raw\GSE14017_series_matrix.txt.gz",'rt',errors='ignore') as f:
    lines=f.readlines()
for i,l in enumerate(lines):
    if not l.startswith("!"):
        print("first data line:",repr(l[:150]))
        print("next:",repr(lines[i+1][:150]))
        break
