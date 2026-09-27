"""
00_download.py
Download GEO series-matrix files for the BCBM competence-vs-adaptation study.

Datasets (per study design):
  GSE184869   45 paired primary breast cancer -> brain metastasis (centerpiece)
  GSE38057    87 HER2+ primary tumors with brain-metastasis-free survival
  GSE12276    BMFS cohort (validation)
  GSE2034     BMFS cohort (validation)
  GSE2603     BMFS cohort (validation)
  GSE5327     BMFS cohort (validation)
  GSE14020    BMFS cohort (validation)
  GSE14017    multi-organ metastasis (brain/liver/lung/bone)
  GSE14018    multi-organ metastasis (platform companion)
  GSE125989   external BCBM / primary validation
  GSE43837    external BM validation
  GSE322943   156 BCBM clinical cohort
  GSE324453   2026 BCBM snRNA-seq (21 samples) - processed count matrix if available
"""
import os, gzip, shutil, urllib.request, ssl, sys

DATA_DIR = r"D:\BCBM_Project\data\raw"
os.makedirs(DATA_DIR, exist_ok=True)

# series matrix URL: https://ftp.ncbi.nlm.nih.gov/geo/series/GSExxxnnn/GSEyyy/matrix/GSEyyy_series_matrix.txt.gz
def geo_matrix_url(gse):
    n = gse[:-3] + "nnn"          # e.g. GSE184nnn
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{n}/{gse}/matrix/{gse}_series_matrix.txt.gz"

DATASETS = [
    "GSE184869", "GSE38057", "GSE12276", "GSE2034", "GSE2603",
    "GSE5327", "GSE14020", "GSE14017", "GSE14018", "GSE125989",
    "GSE43837", "GSE322943",
]

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"[skip] {os.path.basename(dest)} exists ({os.path.getsize(dest)} bytes)")
        return True
    print(f"[get ] {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, context=ctx, timeout=120) as r, open(dest, "wb") as f:
            shutil.copyfileobj(r, f, length=1 << 16)
        print(f"[ok  ] {os.path.basename(dest)} -> {os.path.getsize(dest)} bytes")
        return True
    except Exception as e:
        print(f"[FAIL] {url}: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False

if __name__ == "__main__":
    results = {}
    for gse in DATASETS:
        dest = os.path.join(DATA_DIR, f"{gse}_series_matrix.txt.gz")
        results[gse] = fetch(geo_matrix_url(gse), dest)
    print("\n=== summary ===")
    for k, v in results.items():
        print(f"{k}: {'OK' if v else 'FAILED'}")
