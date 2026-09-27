import gzip
for n in ["GSE125989","GSE43837","GSE14017","GSE14018"]:
    with gzip.open(rf"D:\BCBM_Project\data\raw\{n}_series_matrix.txt.gz",'rt',errors='ignore') as f:
        for line in f:
            if "platform_id" in line.lower() or line.startswith("!Sample_platform"):
                print(n, line.strip()[:120]); break
