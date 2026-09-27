import gzip
with gzip.open(r"D:\BCBM_Project\data\raw\GSE14017_series_matrix.txt.gz",'rt',errors='ignore') as f:
    started=False
    for line in f:
        if line.startswith("!series_matrix_table_begin"): started=True; continue
        if started:
            print(repr(line[:200])); break
