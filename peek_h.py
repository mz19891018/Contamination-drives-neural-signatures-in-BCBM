import gzip
with gzip.open(r"D:\BCBM_Project\data\raw\GSE14017_series_matrix.txt.gz",'rt',errors='ignore') as f:
    for l in f:
        if l.startswith("ID_REF") or "ID_REF" in l[:10]:
            print("HEADER:",repr(l[:200])); break
