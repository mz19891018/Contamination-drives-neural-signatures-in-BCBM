import os, sys, gzip
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import RAW, read_series_matrix

expr, meta = read_series_matrix(os.path.join(RAW,"GSE14017_series_matrix.txt.gz"))
print("expr index head:", list(expr.index[:10]))

# platform map
path = os.path.join(RAW,"GPL570.annot.gz")
rows={}
with gzip.open(path,"rt",errors="replace",encoding="utf-8") as f:
    header=None
    for line in f:
        if line.startswith(("#","!","^")): continue
        parts=line.rstrip("\n").split("\t")
        if header is None:
            header=parts
            id_i=0; sym_i=2
            continue
        rows[parts[id_i]] = parts[sym_i].split(" /// ")[0].strip()
s=pd.Series(rows)
print("platform map size:", len(s))
print("platform head:", s.head())
print("overlap with expr index:", len(set(expr.index) & set(s.index)))
