import gzip
path = r"D:\BCBM_Project\data\raw\GPL570.annot.gz"
with gzip.open(path,"rt",errors="replace",encoding="utf-8") as f:
    lines = f.readlines()
# find first non-comment
for i,l in enumerate(lines):
    if not l.startswith(("#","!","^")):
        print("HEADER line", i, ":", l[:300])
        print("NEXT data row:", lines[i+1][:300])
        break
