import pandas as pd, numpy as np, re, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
D=r"D:\BCBM_Project"; FIG=rf"{D}\manuscript\figures"

H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
DE=pd.read_csv(rf"{D}\results\tables\paired_de_gse184869.csv")
DE=DE.dropna(subset=['symbol'])
brain_markers=set(H[H.verdict=="contamination_driven"].gene.tolist())
competence_genes=set([g.strip() for g in open(rf"{D}\results\tables\geneset_competence.txt",encoding="utf-8") if g.strip()])

fig,ax=plt.subplots(1,1,figsize=(8,6))
de=DE.copy()
de['neglog10fdr']=-np.log10(de['FDR'].clip(lower=1e-300))
def vclass(g):
    if g in brain_markers: return 'contamination'
    if g in competence_genes: return 'competence'
    return 'other'
de['cls']=de['symbol'].map(vclass)
colors={'other':'#CCCCCC','contamination':'#C44E52','competence':'#4C72B0'}
for c in ['other','competence','contamination']:
    sub=de[de.cls==c]
    # Use full set size for competence label (1050), not intersection with DE (1048)
    nlabel = len(competence_genes) if c=='competence' else len(sub)
    ax.scatter(sub['mean_logFC_BM_vs_BP'],sub['neglog10fdr'],s=8,c=colors[c],alpha=0.6 if c=='other' else 0.8,label=f"{c} (n={nlabel})",zorder=2 if c!='other' else 1)
top=de.nlargest(8,'neglog10fdr')
for _,r in top.iterrows():
    if r['symbol'] in ['C1QL1','ATP1A2','TUBB4A','MAG','APLP1','GNAO1','DEPTOR','TRIM9']:
        ax.annotate(r['symbol'],(r['mean_logFC_BM_vs_BP'],r['neglog10fdr']),fontsize=8,xytext=(5,5),textcoords='offset points')
ax.axhline(-np.log10(0.05),ls='--',c='gray',lw=0.8)
ax.axvline(0.3,ls='--',c='gray',lw=0.8); ax.axvline(-0.3,ls='--',c='gray',lw=0.8)
ax.set_xlabel("log2FC (Brain metastasis vs Primary)")
ax.set_ylabel("-log10(FDR)")
ax.set_title("Paired differential expression (n=21 pairs)\nColoured by per-gene contamination verdict")
ax.legend(loc='upper left',fontsize=9,markerscale=3)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure2_volcano.png",dpi=200,bbox_inches='tight')
plt.close()
print(f"Figure2 regenerated. competence full set n={len(competence_genes)}, in DE n={len(de[de.cls=='competence'])}")
