# 22_contamination.py - normal brain contamination diagnosis on GSE184869
import pandas as pd, numpy as np

D = r"D:\BCBM_Project"
expr = pd.read_excel(rf"{D}\data\raw\GSE184869_expression.xlsx", sheet_name="log2TMMCPM")
expr = expr.set_index("ensembl_gene_id")

# ensembl -> symbol map
pm = pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")[["ensembl","symbol"]].dropna()
emap = dict(zip(pm.ensembl, pm.symbol))
sym = expr.index.map(emap)
expr = expr.assign(symbol=sym).dropna(subset=["symbol"])
expr = expr.groupby("symbol").mean()

# canonical normal-brain cell markers
markers = {
 "Neuron":   ["SNAP25","SYT1","SYP","RBFOX3","STMN2","SYN1","GRIN1","GAP43","DNM1","VAMP2"],
 "Astrocyte":["GFAP","AQP4","ALDH1L1","SLC1A2","GJA1","AGT","SLC1A3","GLUL"],
 "Oligodendro":["MBP","PLP1","MAG","MOG","MOBP","OPALIN","CNP"],
 "Microglia":["P2RY12","CX3CR1","AIF1","TMEM119","CSF1R","ITGAM"],
 "Endo":     ["CLDN5","VWF","PECAM1","KDR","FLT1"],
}
available = {k:[g for g in v if g in expr.index] for k,v in markers.items()}
print("marker coverage:", {k:len(v) for k,v in available.items()}, flush=True)

# per-sample signature score: mean z-score of markers
X = expr.values
Z = (X - X.mean(1,keepdims=True)) / (X.std(1,keepdims=True)+1e-9)
Zdf = pd.DataFrame(Z, index=expr.index, columns=expr.columns)
contam = pd.DataFrame({k: Zdf.loc[v].mean() for k,v in available.items()})
contam["brain_total"] = contam[["Neuron","Astrocyte","Oligodendro"]].mean(1)

# classify samples: BM vs BP vs RCS
cols = list(expr.columns)
def kind(c):
    if c.startswith("BM"): return "BrainMet"
    if c.startswith("BP"): return "Primary"
    return "Other"
contam["kind"] = [kind(c) for c in cols]

print("\n=== mean contamination score by sample group ===", flush=True)
print(contam.groupby("kind")[["Neuron","Astrocyte","Oligodendro","brain_total"]].mean().round(3).to_string(), flush=True)

# adaptation / competence module score per sample
adapt = pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
adapt_genes = set(adapt.loc[adapt.program=="Adaptation_de_novo","symbol"]) & set(expr.index)
comp_genes  = set(adapt.loc[adapt.program=="Competence_preexisting","symbol"]) & set(expr.index)
contam["adapt_mod"] = Zdf.loc[list(adapt_genes)].mean()
contam["comp_mod"]  = Zdf.loc[list(comp_genes)].mean()

# CORRELATION: does adaptation module score track normal-brain contamination?
import scipy.stats as st
bm = contam[contam.kind=="BrainMet"]
print("\n=== in BrainMet samples: correlation adapt_mod vs contamination ===", flush=True)
for c in ["Neuron","Astrocyte","Oligodendro","Microglia","Endo","brain_total"]:
    r,p = st.spearmanr(bm.adapt_mod, bm[c])
    print(f"adapt_mod vs {c:12s}: rho={r:+.3f} p={p:.4g}", flush=True)
print("\n--- comp_mod vs contamination (should be weaker) ---", flush=True)
for c in ["Neuron","Astrocyte","Oligodendro","brain_total"]:
    r,p = st.spearmanr(bm.comp_mod, bm[c])
    print(f"comp_mod  vs {c:12s}: rho={r:+.3f} p={p:.4g}", flush=True)

# key: top adaptation genes -- are they also top ODG/astro markers?
print("\n=== are top adaptation genes canonical normal-brain markers? ===", flush=True)
top_adapt = adapt.sort_values("mean_logFC_BM_vs_BP",ascending=False).symbol.head(30).tolist()
allbrain = set(sum(available.values(),[]))
overlap = [g for g in top_adapt if g in allbrain]
print("top30 adaptation genes that are canonical brain-cell markers:", overlap, flush=True)

contam.to_csv(rf"{D}\results\tables\contamination_scores_bulk.csv")
print("\nsaved.", flush=True)
