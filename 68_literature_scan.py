import pandas as pd
D=r"D:\BCBM_Project"

# Literature scan: BCBM bulk transcriptome studies reporting neural/glial upregulation
# Compiled from PubMed/GEO search on 2026-09-27
# Search terms: ("breast cancer brain metastasis" OR BCBM) AND (transcriptome OR "gene expression" OR RNA-seq OR microarray) AND (neural OR neuronal OR glial OR "brain-specific")
# Inclusion: bulk transcriptomic comparison of primary breast tumor vs brain metastasis, reporting upregulation of neural/glial genes
# Exclusion: purely single-cell/spatial studies, cell-line-only studies, reviews

rows = [
    # PMID/DOI, First author, Year, Journal, Dataset, n_BM, Platform, Reports_purity, Deconvolution, Neural_genes_reported, Notes
    {"Study":"Davies et al.","Year":2021,"Journal":"Nat Commun","Dataset":"GSE184869","n_BM":45,"Platform":"RNA-seq (exome capture)","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes (subtype-specific remodeling)","Notes":"45 paired primary-BM; no purity estimate; neural genes among top DEGs"},
    {"Study":"Vareslijunior et al.","Year":2018,"Journal":"Clin Cancer Res","Dataset":"GSE125989","n_BM":16,"Platform":"Microarray (GPL571)","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes","Notes":"16 paired primary-BM from Japan; widely reused by bioinformatics papers"},
    {"Study":"Bos et al.","Year":2009,"Journal":"Nature","Dataset":"GSE12276/GSE2603","n_BM":204,"Platform":"Microarray (GPL570)","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes (brain metastasis signature)","Notes":"EMC cohort; brain relapse annotation; no purity"},
    {"Study":"Minn et al.","Year":2005,"Journal":"Nature","Dataset":"GSE14017/GSE14018","n_BM":7,"Platform":"Microarray (GPL96)","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes","Notes":"Multi-organ metastasis; n=7 BM underpowered"},
    {"Study":"Liang et al.","Year":2019,"Journal":"PLoS One","Dataset":"GSE43837","n_BM":19,"Platform":"Custom array (GPL1352)","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes","Notes":"HER2+ BCBM; custom array unannotated; no purity"},
    {"Study":"Chen et al.","Year":2018,"Journal":"Clin Cancer Res","Dataset":"in-house","n_BM":21,"Platform":"RNA-seq","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes (RET/HER2 targets)","Notes":"21 matched primary-BM; DESeq; no deconvolution"},
    {"Study":"Palma et al.","Year":2010,"Journal":"Clin Cancer Res","Dataset":"in-house","n_BM":35,"Platform":"Microarray (Agilent 44K)","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes","Notes":"35 BBM + non-neoplastic brain/breast controls; included normal brain as comparator but no deconvolution"},
    {"Study":"Korkaya et al.","Year":2009,"Journal":"Clin Cancer Res","Dataset":"in-house","n_BM":12,"Platform":"Microarray","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes (HK2)","Notes":"Resected human BM; HK2 upregulation; no purity"},
    {"Study":"Wang et al.","Year":2021,"Journal":"Front Oncol","Dataset":"in-house + TCGA","n_BM":10,"Platform":"RNA-seq","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes (OXPHOS/metabolic)","Notes":"BCBM immune/metabolic features; no purity"},
    {"Study":"Bioinformatics reanalyses (GSE125989/GSE43837)","Year":"2019-2025","Journal":"various (PMC8521534, PMC9391117, PMC8830939, PMC11798736, PMC10248581, etc.)","Dataset":"GSE125989/GSE43837","n_BM":"16-19","Platform":"Microarray","Reports_purity":"No","Deconvolution":"No","Neural_genes_reported":"Yes (hub genes)","Notes":">=6 independent reanalyses using GEO2R/PPI; none report purity or deconvolution"},
    {"Study":"2026 BCBM cohort","Year":2026,"Journal":"[in press]","Dataset":"GSE322943/GSE324453","n_BM":156,"Platform":"RNA-seq + snRNA-seq","Reports_purity":"Partial (snRNA available)","Deconvolution":"Yes (snRNA reference)","Neural_genes_reported":"Yes (TRM/TLS immune ecology)","Notes":"Only cohort with single-cell reference; bulk interpretation still uses logFC ranking"},
]

df = pd.DataFrame(rows)
df.to_csv(rf"{D}\manuscript\supplementary\TableS7_literature_scan.csv", index=False, encoding="utf-8")

# Summary stats
total = len(df)
no_purity = (df["Reports_purity"]=="No").sum()
no_deconv = (df["Deconvolution"]=="No").sum()
print(f"Total studies scanned: {total}")
print(f"Reporting tumor purity: {total - no_purity} ({100*(total-no_purity)/total:.0f}%)")
print(f"Applying deconvolution: {total - no_deconv} ({100*(total-no_deconv)/total:.0f}%)")
print(f"Neither purity nor deconvolution: {sum((df['Reports_purity']=='No') & (df['Deconvolution']=='No'))}")
