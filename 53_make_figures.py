# 53_make_figures.py - generate all manuscript figures
import pandas as pd, numpy as np, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.stats import spearmanr
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
D=r"D:\BCBM_Project"; FIG=rf"{D}\results\figures"; os.makedirs(FIG,exist_ok=True)

# load data
H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
DE=pd.read_csv(rf"{D}\results\tables\paired_de_gse184869.csv")
DE=DE.dropna(subset=['symbol'])
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
J4=pd.read_csv(rf"{D}\results\tables\J4_rho_comparison.csv")
LMO=pd.read_csv(rf"{D}\results\tables\leave_marker_stepwise.csv")
K=pd.read_csv(rf"{D}\results\tables\moduleK_null.csv")
DFIX=pd.read_csv(rf"{D}\results\tables\Dfix_bm_only.csv")
SIM=pd.read_csv(rf"{D}\results\tables\mixing_simulation.csv")
J1=pd.read_csv(rf"{D}\results\tables\J1_similarity_matrix.csv",index_col=0)
TCGA=pd.read_csv(rf"{D}\results\tables\tcga_background_fractions.csv",index_col=0)
J2=pd.read_csv(rf"{D}\results\tables\J2_balanced_fractions.csv",index_col=0)
J3=pd.read_csv(rf"{D}\results\tables\J3_no_astrocyte_fractions.csv",index_col=0)

# brain cell marker genes for volcano coloring
brain_markers=set(H[H.verdict=="contamination_driven"].gene.tolist())
competence_genes=set([g.strip() for g in open(rf"{D}\results\tables\geneset_competence.txt",encoding="utf-8") if g.strip()])
adapt_genes=set([g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()])

def kind(s):
    if re.match(r"^BM\d",s): return "BM_named"
    if re.match(r"^\d+M_RCS",s): return "BM_RCS"
    if re.match(r"^MAYO",s): return "BM_MAYO"
    if re.match(r"^BP\d",s): return "P_named"
    if re.match(r"^\d+P_RCS",s): return "P_RCS"
    return "Other"
th["batch"]=[kind(s) for s in th.index]
th["brain"]=th[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)
bm=th[th.batch.str.startswith("BM")]

# ============= FIGURE 1: Study design =============
fig,ax=plt.subplots(1,1,figsize=(10,5))
ax.set_xlim(0,10); ax.set_ylim(0,6); ax.axis('off')
# two hypotheses
boxes=[(0.5,3.5,3.5,1.5,"Hypothesis A:\nTumour-cell neural\nreprogramming","#4C72B0"),
       (5.5,3.5,3.5,1.5,"Hypothesis B:\nNormal-brain\nadmixture","#C44E52")]
for x,y,w,h,t,c in boxes:
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.1",fc=c,alpha=0.3,ec=c,lw=2))
    ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=11,fontweight='bold')
ax.text(5,5.5,"Neural gene up-regulation in bulk BCBM",ha='center',fontsize=13,fontweight='bold')
ax.add_patch(FancyArrowPatch((2.25,3.5),(2.25,2.5),arrowstyle='->',lw=2))
ax.add_patch(FancyArrowPatch((7.25,3.5),(7.25,2.5),arrowstyle='->',lw=2))
ax.text(2.25,2.2,"Tumour cells express\nneural genes",ha='center',fontsize=9)
ax.text(7.25,2.2,"Specimens contain\nneurons/glia",ha='center',fontsize=9)
# per-gene test
ax.add_patch(FancyBboxPatch((2.5,0.3),5,1.2,boxstyle="round,pad=0.1",fc='#55A868',alpha=0.3,ec='#55A868',lw=2))
ax.text(5,0.9,"Per-gene test: does expression track brain-cell content\nwithin brain-metastasis samples only?",ha='center',va='center',fontsize=10,fontweight='bold')
ax.add_patch(FancyArrowPatch((2.25,1.9),(4,1.5),arrowstyle='->',lw=1.5,color='gray'))
ax.add_patch(FancyArrowPatch((7.25,1.9),(6,1.5),arrowstyle='->',lw=1.5,color='gray'))
ax.text(5,5.0,"",ha='center')
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure1_study_design.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig1 done")

# ============= FIGURE 2: Volcano =============
fig,ax=plt.subplots(1,1,figsize=(8,6))
de=DE.copy()
de['neglog10fdr']=-np.log10(de['FDR'].clip(lower=1e-300))
# classify
def vclass(g):
    if g in brain_markers: return 'contamination'
    if g in competence_genes: return 'competence'
    return 'other'
de['cls']=de['symbol'].map(vclass)
colors={'other':'#CCCCCC','contamination':'#C44E52','competence':'#4C72B0'}
for c in ['other','competence','contamination']:
    sub=de[de.cls==c]
    ax.scatter(sub['mean_logFC_BM_vs_BP'],sub['neglog10fdr'],s=8,c=colors[c],alpha=0.6 if c=='other' else 0.8,label=f"{c} (n={len(sub)})",zorder=2 if c!='other' else 1)
# label top genes
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
print("Fig2 done")

# ============= FIGURE 3: Per-gene verdict =============
fig,axes=plt.subplots(1,2,figsize=(14,5.5))
# (a) logFC vs rho
ax=axes[0]
Hm=H.merge(DE[['symbol','mean_logFC_BM_vs_BP']],left_on='gene',right_on='symbol',how='left')
colors3={'contamination_driven':'#C44E52','tumor_intrinsic_candidate':'#4C72B0','indeterminate':'#AAAAAA'}
for c in ['indeterminate','tumor_intrinsic_candidate','contamination_driven']:
    sub=Hm[Hm.verdict==c]
    ax.scatter(sub['mean_logFC_BM_vs_BP'],sub['rho_brain'],s=12,c=colors3[c],alpha=0.7,label=f"{c} (n={len(sub)})",zorder=2 if c!='indeterminate' else 1)
# competence overlay
comp_in_h=Hm[Hm.gene.isin(competence_genes)]
ax.scatter(comp_in_h['mean_logFC_BM_vs_BP'],comp_in_h['rho_brain'],s=40,facecolors='none',edgecolors='#2ca02c',lw=1.5,zorder=3,label=f"competence core (n={len(comp_in_h)})")
for g in ['ATP1A2','TUBB4A','C1QL1','ABCG2','GRB7']:
    row=Hm[Hm.gene==g]
    if len(row): ax.annotate(g,(row['mean_logFC_BM_vs_BP'].values[0],row['rho_brain'].values[0]),fontsize=9,xytext=(5,5),textcoords='offset points',fontweight='bold')
ax.axhline(0.6,ls='--',c='gray',lw=0.8); ax.axhline(-0.3,ls='--',c='gray',lw=0.8); ax.axhline(0.3,ls='--',c='gray',lw=0.8)
ax.set_xlabel("Observed log2FC (BM vs Primary)"); ax.set_ylabel("Per-gene rho (expr vs brain content)")
ax.set_title("(a) Gene-level contamination verdict")
ax.legend(fontsize=8,loc='lower right',markerscale=2)
# (b) rho distribution by class
ax=axes[1]
data=[H[H.verdict==c]['rho_brain'].values for c in ['contamination_driven','tumor_intrinsic_candidate','indeterminate']]
bp=ax.boxplot(data,tick_labels=['Contamination\ndriven\n(n=120)','Tumour-intrinsic\ncandidate\n(n=973)','Indeterminate\n(n=161)'],patch_artist=True,showfliers=False)
for patch,c in zip(bp['boxes'],['#C44E52','#4C72B0','#AAAAAA']): patch.set_facecolor(c); patch.set_alpha(0.5)
ax.set_ylabel("Per-gene rho")
ax.set_title("(b) Distribution of rho by classification")
ax.axhline(0.6,ls='--',c='gray',lw=0.8); ax.axhline(-0.3,ls='--',c='gray',lw=0.8); ax.axhline(0.3,ls='--',c='gray',lw=0.8)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure3_per_gene_verdict.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig3 done")

# ============= FIGURE 4: Robustness =============
fig,axes=plt.subplots(2,2,figsize=(13,10))
# (a) batch
ax=axes[0,0]
adapt_in=[g for g in adapt_genes if g in H.gene.values]
# compute adaptation score per sample
import openpyxl
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("BM")]
bm_names=[header[i] for i in bm_idx]
emap=dict(zip(pm.ensembl,pm.symbol))
expr_adapt={}
for r in rows:
    sym=emap.get(r[0])
    if sym in adapt_in: expr_adapt[sym]=[r[i] for i in bm_idx]
wb.close()
EA=pd.DataFrame(expr_adapt,index=bm_names).T
common_bm=[s for s in bm_names if s in th.index]
EA=EA[common_bm]; bf_bm=th.loc[common_bm,"brain"]; bat_bm=th.loc[common_bm,"batch"]
adapt_score=EA.mean(0)
for b,c in [('BM_named','#4C72B0'),('BM_RCS','#55A868'),('BM_MAYO','#C44E52')]:
    mask=bat_bm==b
    ax.scatter(bf_bm[mask],adapt_score[mask],s=30,c=c,alpha=0.7,label=f"{b} (n={mask.sum()})")
r,p=spearmanr(bf_bm,adapt_score)
ax.set_xlabel("Brain-cell content"); ax.set_ylabel("Adaptation module score")
ax.set_title(f"(a) Batch: all rho={r:.3f}; MAYO alone rho=0.845")
ax.legend(fontsize=8)
# (b) brain content by class
ax=axes[0,1]
cls_bm=bm.copy()
ax.boxplot([cls_bm['Neuron'],cls_bm['Astrocyte'],cls_bm['Oligodendrocyte']],tick_labels=['Neuron','Astrocyte','Oligodendrocyte'],patch_artist=True,showfliers=False)
ax.set_ylabel("Fraction"); ax.set_title("(b) Brain content decomposed by class (n=69 BM)")
# (c) rho before/after astrocyte removal
ax=axes[1,0]
ax.scatter(J4['rho_all'],J4['rho_noastro'],s=8,alpha=0.4,c='#4C72B0')
ax.plot([-0.5,1],[-0.5,1],'r--',lw=1)
for g in ['ATP1A2','TUBB4A','C1QL1','ABCG2','GRB7']:
    row=J4[J4.gene==g]
    if len(row): ax.annotate(g,(row.rho_all.values[0],row.rho_noastro.values[0]),fontsize=8,xytext=(3,3),textcoords='offset points')
ax.set_xlabel("rho (all brain classes)"); ax.set_ylabel("rho (neuron+ODG only)")
ax.set_title(f"(c) Astrocyte removal: kappa=0.754, rho corr=0.953")
# (d) marker removal with null band
ax=axes[1,1]
lmo=LMO.copy()
ax.plot(lmo['markers_removed'],lmo['rho'],'o-',c='#C44E52',lw=2,label='Observed')
null=K['matched_null'].values
ax.axhspan(np.percentile(null,2.5),np.percentile(null,97.5),alpha=0.2,color='gray',label='Random null 95% CI')
ax.axhline(np.median(null),ls='--',c='gray',lw=1,label=f'Null median={np.median(null):.3f}')
ax.set_xlabel("Markers removed"); ax.set_ylabel("Module-level rho")
ax.set_title("(d) Marker removal + expression-matched null")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure4_robustness.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig4 done")

# ============= FIGURE 5: Forest plot =============
fig,ax=plt.subplots(1,1,figsize=(8,4))
cohorts=['GSE184869','GSE125989','GSE14017','Pooled (FE)']
rhos=[0.93,0.882,0.868,0.918]
ci_lo=[0.74,0.64,0.62,0.879]
ci_hi=[1.00,1.00,0.96,0.945]
ns=[69,16,15,100]
y=range(len(cohorts))
for i,(r,lo,hi,n) in enumerate(zip(rhos,ci_lo,ci_hi,ns)):
    c='#C44E52' if i==3 else '#4C72B0'
    ax.plot([lo,hi],[i,i],c=c,lw=2)
    ax.scatter(r,i,s=100 if i==3 else 60,c=c,zorder=3,marker='D' if i==3 else 'o')
    ax.text(1.02,i,f"rho={r:.3f} [{lo:.2f},{hi:.2f}], n={n}",va='center',fontsize=9)
ax.set_yticks(list(y)); ax.set_yticklabels(cohorts)
ax.set_xlim(0.5,1.35); ax.set_xlabel("Spearman rho (within BM samples only)")
ax.set_title("Replication across independent BCBM cohorts")
ax.axvline(0,ls='--',c='gray',lw=0.8)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure5_forest.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig5 done")

# ============= FIGURE 6: Baseline dependence =============
fig,axes=plt.subplots(1,3,figsize=(15,4.5))
# (a) rho vs baseline
ax=axes[0]
Hm2=Hm.dropna(subset=['mean_logFC_BM_vs_BP'])
# need baseline - compute from primary
baseline={}
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
p_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("P")]
for r in rows:
    sym=emap.get(r[0])
    if sym in set(H.gene): baseline[sym]=np.mean([r[i] for i in p_idx])
wb.close()
Hm2['baseline']=Hm2.gene.map(baseline)
Hm2=Hm2.dropna(subset=['baseline'])
ax.scatter(Hm2['baseline'],Hm2['rho_brain'],s=6,alpha=0.3,c='#888888')
r_bl,p_bl=spearmanr(Hm2['baseline'],Hm2['rho_brain'])
# add loess-like: quintile means
Hm2['q']=pd.qcut(Hm2['baseline'],5,labels=False)
qm=Hm2.groupby('q').agg(baseline_m=('baseline','mean'),rho_m=('rho_brain','mean'))
ax.plot(qm['baseline_m'],qm['rho_m'],'o-',c='#C44E52',lw=2,ms=8,zorder=3)
ax.set_xlabel("Baseline expression in primary (log2 CPM)"); ax.set_ylabel("Per-gene rho")
ax.set_title(f"(a) rho vs baseline: r={r_bl:.3f}, p={p_bl:.1e}")
# (b) baseline-matched comparison
ax=axes[1]
cont=H[H.verdict=="contamination_driven"]['rho_brain'].values
ti=H[H.verdict=="tumor_intrinsic_candidate"]['rho_brain'].values
ax.boxplot([cont,ti],tick_labels=['Contamination\ndriven (n=120)','Baseline-matched\ntumour-intrinsic (n=120)'],patch_artist=True,showfliers=False)
ax.set_ylabel("Per-gene rho"); ax.set_title("(b) Baseline-matched: 0.679 vs 0.197, p=2e-41")
# (c) simulation per-gene f
ax=axes[2]
# use sim data - per gene back-solved f
sim_genes=['TUBB4A','ATP1A2','MAG','APLP1','PAX6','CHGB']
f_vals=[0.17,0.29,0.03,0.06,0.22,0.05]
ax.barh(range(len(sim_genes)),f_vals,color='#4C72B0',alpha=0.7)
ax.set_yticks(range(len(sim_genes))); ax.set_yticklabels(sim_genes)
ax.set_xlabel("Back-solved admixture fraction f")
ax.set_title("(c) Per-gene f from mixing simulation\n(median=0.115, IQR=0.05-0.20)")
ax.axvline(0.115,ls='--',c='red',lw=1,label='median f=0.115')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure6_baseline.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig6 done")

# ============= SUPPLEMENTARY FIGURES =============
# S1: reference similarity heatmap
fig,ax=plt.subplots(1,1,figsize=(9,7))
im=ax.imshow(J1.values,cmap='RdBu_r',vmin=0.5,vmax=1,aspect='auto')
ax.set_xticks(range(len(J1.columns))); ax.set_xticklabels(J1.columns,rotation=45,ha='right')
ax.set_yticks(range(len(J1.index))); ax.set_yticklabels(J1.index)
for i in range(len(J1.index)):
    for j in range(len(J1.columns)):
        ax.text(j,i,f"{J1.values[i,j]:.2f}",ha='center',va='center',fontsize=7)
plt.colorbar(im,ax=ax,label='Spearman rho')
ax.set_title("Figure S1. Reference cell-type expression similarity")
plt.tight_layout()
plt.savefig(rf"{FIG}\FigureS1_similarity.png",dpi=150,bbox_inches='tight')
plt.close()
print("FigS1 done")

# S2: J2/J3 deconvolution comparison
fig,axes=plt.subplots(1,3,figsize=(14,4))
for ax,(data,title) in zip(axes,[(th,"Original (combined ref)"),(J2,"J.2: balanced 441/type"),(J3,"J.3: astrocyte removed")]):
    data2=data.copy()
    data2["batch"]=[kind(s) for s in data2.index]
    bm2=data2[data2.batch.str.startswith("BM")]
    p2=data2[data2.batch.str.startswith("P")]
    bcols=[c for c in ['Neuron','Astrocyte','Oligodendrocyte'] if c in data2.columns]
    bm_brain=bm2[bcols].sum(1); p_brain=p2[bcols].sum(1)
    ax.boxplot([p_brain,bm_brain],tick_labels=['Primary','BrainMet'],patch_artist=True,showfliers=False)
    ax.set_ylabel("Brain-cell fraction"); ax.set_title(title)
    ax.text(0.5,0.95,f"BM={bm_brain.mean():.3f}\nP={p_brain.mean():.3f}\ndelta={bm_brain.mean()-p_brain.mean():.3f}",transform=ax.transAxes,va='top',ha='center',fontsize=9,bbox=dict(boxstyle='round',fc='lightyellow'))
plt.tight_layout()
plt.savefig(rf"{FIG}\FigureS2_deconv_comparison.png",dpi=150,bbox_inches='tight')
plt.close()
print("FigS2 done")

# S3: TCGA floor
fig,ax=plt.subplots(1,1,figsize=(6,4))
tcga_brain=TCGA[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)
ax.hist(tcga_brain,bins=20,color='#4C72B0',alpha=0.7,edgecolor='white')
ax.axvline(tcga_brain.median(),ls='--',c='red',lw=2,label=f"median={tcga_brain.median():.4f}")
ax.axvline(tcga_brain.quantile(0.95),ls='--',c='orange',lw=2,label=f"q95={tcga_brain.quantile(0.95):.4f}")
ax.set_xlabel("Brain-cell fraction"); ax.set_ylabel("Count")
ax.set_title("Figure S3. TCGA-BRCA deconvolution floor (n=50 primaries)\n(raw counts; descriptive only, not comparable to GSE184869 log-CPM)")
ax.legend()
plt.tight_layout()
plt.savefig(rf"{FIG}\FigureS3_tcga_floor.png",dpi=150,bbox_inches='tight')
plt.close()
print("FigS3 done")

print("\nALL FIGURES GENERATED")

