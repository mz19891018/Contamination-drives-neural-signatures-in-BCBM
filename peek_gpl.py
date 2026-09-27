import gzip
with gzip.open(r"D:\BCBM_Project\data\raw\GPL570.annot.gz",'rt',errors='ignore') as f:
    for i,line in enumerate(f):
        print(line.rstrip()[:200])
        if i>15: break
