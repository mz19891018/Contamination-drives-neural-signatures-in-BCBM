"""
common.py - shared loaders for GEO series-matrix files and annotation.
"""
import os, gzip, re, numpy as np, pandas as pd, sys

RAW = r"D:\BCBM_Project\data\raw"
RESULTS = r"D:\BCBM_Project\results"
PYLIBS = r"D:\BCBM_Project\pylibs"
if os.path.isdir(PYLIBS) and PYLIBS not in sys.path:
    sys.path.insert(0, PYLIBS)

def read_series_matrix(path):
    """Parse a GEO series_matrix.txt.gz into (expr_df, meta_df).
    expr_df: index=ID_REF, columns=samples.
    meta_df: index=samples, columns=parsed !Sample_* fields.
    """
    rows_meta = {}
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", errors="replace", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line.rstrip("\n"))
                continue
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                key = parts[0]
                vals = [v.strip().strip('"') for v in parts[1:]]
                rows_meta[key] = vals
    # meta: transpose -> rows=samples
    meta = pd.DataFrame(rows_meta).set_index("!Sample_geo_accession") if "!Sample_geo_accession" in rows_meta else pd.DataFrame(rows_meta)
    # expr table
    from io import StringIO
    expr = pd.read_csv(StringIO("\n".join(table_lines)), sep="\t", index_col=0)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr, meta


# ---------- Ensembl -> symbol mapping (cached) ----------
_MG = None
def _mg():
    global _MG
    if _MG is None:
        import mygene
        _MG = mygene.MyGeneInfo()
    return _MG

_cache_path = os.path.join(RESULTS, "ensembl_symbol_map.csv")
def ensembl_to_symbol(ens_ids):
    """Return dict ensg->symbol, cached on disk."""
    if os.path.exists(_cache_path):
        m = pd.read_csv(_cache_path)
        have = dict(zip(m.ensembl, m.symbol))
    else:
        have = {}
    missing = [e for e in ens_ids if e not in have]
    if missing:
        print(f"[map] querying {len(missing)} ensembl ids ...")
        res = _mg().querymany(missing, scopes="ensembl.gene",
                              fields="symbol,name", species="human", verbose=False)
        new = {}
        for e, r in zip(missing, res):
            sym = r.get("symbol") if isinstance(r, dict) else None
            new[e] = sym if sym else np.nan
        df = pd.DataFrame({"ensembl": list(new.keys()), "symbol": list(new.values())})
        if os.path.exists(_cache_path):
            df.to_csv(_cache_path, mode="a", header=False, index=False)
        else:
            df.to_csv(_cache_path, index=False)
        have.update(new)
    return {k: have.get(k) for k in ens_ids}
