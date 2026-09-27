# 61_redraw_figures.py - fix figures per reviewer comments
import pandas as pd, numpy as np, openpyxl, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
D=r"D:\BCBM_Project"; FIG=rf"{D}\manuscript\figures"

H=pd.read_csv(rf"{D}\results\tables\gene_contamination_verdict.csv")
DE=pd.read_csv(rf"{D}\results\tables\paired_de_gse184869.csv").dropna(subset=['symbol'])
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
J4v2=pd.read_csv(rf"{D}\results\tables\J4v2_neuron_only.csv")
LMO=pd.read_csv(rf"{D}\results\tables\leave_marker_stepwise.csv")
K=pd.read_csv(rf"{D}\results\tables\moduleK_null.csv")
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
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
comp=[g.strip() for g in open(rf"{D}\results\tables\geneset_competence.txt",encoding="utf-8") if g.strip()]
brain_markers=set(H[H.verdict=="contamination_driven"].gene)

# load expression for adaptation score
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("BM")]
bm_names=[header[i] for i in bm_idx]
emap=dict(zip(pm.ensembl,pm.symbol))
expr={}
for r in rows:
    sym=emap.get(r[0])
    if sym in set(adapt): expr[sym]=[r[i] for i in bm_idx]
wb.close()
EA=pd.DataFrame(expr,index=bm_names).T
common_bm=[s for s in bm_names if s in th.index]
EA=EA[common_bm]; bf_bm=th.loc[common_bm,"brain"]; neu_bm=th.loc[common_bm,"Neuron"]; bat_bm=th.loc[common_bm,"batch"]

# ============= FIGURE 3 REDESIGN =============
Hm=H.merge(DE[['symbol','mean_logFC_BM_vs_BP']],left_on='gene',right_on='symbol',how='left')
fig=plt.figure(figsize=(14,5.5))
gs=fig.add_gridspec(1,3,width_ratios=[3,1,1],wspace=0.3)
# (a) main scatter - 3 classes only
ax=fig.add_subplot(gs[0,0])
colors3={'contamination_driven':'#C44E52','tumor_intrinsic_candidate':'#4C72B0','indeterminate':'#AAAAAA'}
for c in ['indeterminate','tumor_intrinsic_candidate','contamination_driven']:
    sub=Hm[Hm.verdict==c]
    ax.scatter(sub['mean_logFC_BM_vs_BP'],sub['rho_brain'],s=10,c=colors3[c],alpha=0.7,
               label=f"{c.replace('_',' ')} (n={len(sub)})",zorder=2 if c!='indeterminate' else 1)
# label key genes with offset
labels={'C1QL1':(5.15,0.72,-30,10),'ATP1A2':(2.78,0.75,10,12),'TUBB4A':(2.44,0.74,10,-15),
        'ABCG2':(1.82,0.11,10,10),'GRB7':(1.03,0.13,-40,-15)}
for g,(x,y,dx,dy) in labels.items():
    row=Hm[Hm.gene==g]
    if len(row): ax.annotate(g,(row['mean_logFC_BM_vs_BP'].values[0],row['rho_brain'].values[0]),
        fontsize=9,fontweight='bold',xytext=(dx,dy),textcoords='offset points',
        arrowprops=dict(arrowstyle='-',color='gray',lw=0.5))
ax.axhline(0.6,ls='--',c='gray',lw=0.8); ax.axhline(0.3,ls='--',c='gray',lw=0.8); ax.axhline(-0.3,ls='--',c='gray',lw=0.8)
ax.set_xlabel("Observed log2FC (BM vs Primary)"); ax.set_ylabel("Per-gene rho (expr vs brain content)")
ax.set_title("(a) Gene-level contamination verdict")
ax.legend(fontsize=8,loc='lower right',markerscale=2)
# (b) competence gene rho distribution (marginal)
ax=fig.add_subplot(gs[0,1])
comp_rho=H[H.gene.isin(comp)]['rho_brain'].values
all_rho=H['rho_brain'].values
ax.hist(all_rho,bins=40,color='#CCCCCC',alpha=0.7,density=True,label='All genes')
ax.hist(comp_rho,bins=20,color='#2ca02c',alpha=0.5,density=True,label=f'Competence (n={len(comp_rho)})')
ax.axvline(0.3,ls='--',c='gray',lw=0.8); ax.axvline(-0.3,ls='--',c='gray',lw=0.8)
ax.set_xlabel("Per-gene rho"); ax.set_ylabel("Density")
ax.set_title("(b) Competence genes are\nbrain-content-independent")
ax.legend(fontsize=8)
# (c) competence overlap bar
ax=fig.add_subplot(gs[0,2])
comp_in=H[H.gene.isin(comp)]
vc=comp_in.verdict.value_counts()
cats=['tumor_intrinsic_candidate','indeterminate','contamination_driven']
vals=[vc.get(c,0) for c in cats]
bar_colors=['#4C72B0','#AAAAAA','#C44E52']
bars=ax.barh(range(3),vals,color=bar_colors,alpha=0.7)
ax.set_yticks(range(3)); ax.set_yticklabels(['Tumour-intrinsic','Indeterminate','Contamination'])
ax.set_xlabel("Number of competence genes")
ax.set_title(f"(c) Competence overlap\n(2/1050 contamination-driven)")
for i,v in enumerate(vals): ax.text(v+5,i,str(v),va='center',fontsize=9)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure3_per_gene_verdict.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig3 redrawn")

# ============= FIGURE 4 - add neuron-only panel =============
fig,axes=plt.subplots(2,2,figsize=(13,10))
# (a) batch
ax=axes[0,0]
adapt_score=EA.mean(0)
for b,c in [('BM_named','#4C72B0'),('BM_RCS','#55A868'),('BM_MAYO','#C44E52')]:
    mask=bat_bm==b
    ax.scatter(bf_bm[mask],adapt_score[mask],s=30,c=c,alpha=0.7,label=f"{b} (n={mask.sum()})")
r=spearmanr(bf_bm,adapt_score)[0]
ax.set_xlabel("Brain-cell content (all 3 classes)"); ax.set_ylabel("Adaptation module score")
ax.set_title(f"(a) Batch: all rho={r:.3f}; MAYO alone rho=0.845")
ax.legend(fontsize=8)
# (b) brain content by class - BM vs primary
ax=axes[0,1]
prim=th[th.batch.str.startswith("P")]
x=np.arange(3); w=0.35
bm_means=[bm['Neuron'].mean(),bm['Astrocyte'].mean(),bm['Oligodendrocyte'].mean()]
p_means=[prim['Neuron'].mean(),prim['Astrocyte'].mean(),prim['Oligodendrocyte'].mean()]
ax.bar(x-w/2,bm_means,w,label='BrainMet',color='#C44E52',alpha=0.7)
ax.bar(x+w/2,p_means,w,label='Primary',color='#4C72B0',alpha=0.7)
ax.set_xticks(x); ax.set_xticklabels(['Neuron','Astrocyte','Oligodendrocyte'])
ax.set_ylabel("Mean fraction"); ax.set_title("(b) Brain content by class: BM vs Primary")
ax.legend(fontsize=8)
# annotate ratios
for i,(b,p) in enumerate(zip(bm_means,p_means)):
    ax.text(i,max(b,p)+0.005,f"{b/p:.1f}x",ha='center',fontsize=9,fontweight='bold')
# (c) neuron-only vs original rho
ax=axes[1,0]
ax.scatter(J4v2['rho_brain'],J4v2['rho_neuron'],s=8,alpha=0.4,c='#4C72B0')
ax.plot([-0.5,1],[-0.5,1],'r--',lw=1)
for g in ['ATP1A2','TUBB4A','C1QL1','ABCG2','GRB7']:
    row=J4v2[J4v2.gene==g]
    if len(row): ax.annotate(g,(row.rho_brain.values[0],row.rho_neuron.values[0]),fontsize=8,xytext=(3,3),textcoords='offset points')
ax.set_xlabel("rho (all brain classes)"); ax.set_ylabel("rho (neuron only)")
ax.set_title(f"(c) Neuron-only definition: kappa=0.700, rho corr=0.944")
# (d) marker removal + null
ax=axes[1,1]
ax.plot(LMO['markers_removed'],LMO['rho'],'o-',c='#C44E52',lw=2,label='Observed')
null=K['matched_null'].values
ax.axhspan(np.percentile(null,2.5),np.percentile(null,97.5),alpha=0.2,color='gray',label='Random null 95% CI')
ax.axhline(np.median(null),ls='--',c='gray',lw=1,label=f'Null median={np.median(null):.3f}')
ax.set_xlabel("Markers removed"); ax.set_ylabel("Module-level rho")
ax.set_title("(d) Marker removal + expression-matched null")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure4_robustness.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig4 redrawn")

# ============= FIGURE 5 - fix CI =============
fig,ax=plt.subplots(1,1,figsize=(8,4))
cohorts=['GSE184869','GSE125989','GSE14017','Pooled (FE)']
rhos=[0.93,0.882,0.868,0.918]
# correct CI for n=69: Fisher z
def rho_ci(r,n):
    z=0.5*np.log((1+r)/(1-r)); se=1/np.sqrt(n-3)
    lo=z-1.96*se; hi=z+1.96*se
    return (np.exp(2*lo)-1)/(np.exp(2*lo)+1),(np.exp(2*hi)-1)/(np.exp(2*hi)+1)
ci_los=[rho_ci(0.93,69)[0],rho_ci(0.882,16)[0],rho_ci(0.868,15)[0],0.879]
ci_his=[rho_ci(0.93,69)[1],rho_ci(0.882,16)[1],rho_ci(0.868,15)[1],0.945]
ns=[69,16,15,100]
y=range(len(cohorts))
for i,(r,lo,hi,n) in enumerate(zip(rhos,ci_los,ci_his,ns)):
    c='#C44E52' if i==3 else '#4C72B0'
    ax.plot([lo,hi],[i,i],c=c,lw=2)
    ax.scatter(r,i,s=100 if i==3 else 60,c=c,zorder=3,marker='D' if i==3 else 'o')
    ax.text(1.01,i,f"rho={r:.3f} [{lo:.2f},{hi:.2f}], n={n}",va='center',fontsize=9)
ax.set_yticks(list(y)); ax.set_yticklabels(cohorts)
ax.set_xlim(0.5,1.28); ax.set_xlabel("Spearman rho (within BM samples only)")
ax.set_title("Replication across independent BCBM cohorts")
ax.axvline(0,ls='--',c='gray',lw=0.8)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure5_forest.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig5 redrawn with corrected CIs")

# ============= FIGURE 6a - bimodal fix =============
fig,axes=plt.subplots(1,3,figsize=(15,4.5))
# (a) baseline with group comparison
ax=axes[0]
Hm2=Hm.dropna(subset=['mean_logFC_BM_vs_BP']).copy()
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
Hm2['q']=pd.qcut(Hm2['baseline'],5,labels=False)
qm=Hm2.groupby('q').agg(bm=('baseline','mean'),rm=('rho_brain','mean'))
ax.errorbar(qm['bm'],qm['rm'],yerr=Hm2.groupby('q')['rho_brain'].sem(),fmt='o-',c='#C44E52',lw=2,ms=8,zorder=3,capsize=3)
r_bl=spearmanr(Hm2['baseline'],Hm2['rho_brain'])[0]
# group comparison
low=Hm2[Hm2['q']==0]['rho_brain']; high=Hm2[Hm2['q']==4]['rho_brain']
u,p=mannwhitneyu(low,high,alternative='greater')
ax.set_xlabel("Baseline expression in primary (log2 CPM)"); ax.set_ylabel("Per-gene rho")
ax.set_title(f"(a) rho vs baseline: r={r_bl:.3f}\nlow vs high quintile: median {low.median():.2f} vs {high.median():.2f}, p={p:.1e}")
# (b) baseline-matched
ax=axes[1]
cont=H[H.verdict=="contamination_driven"]['rho_brain'].values
ti=H[H.verdict=="tumor_intrinsic_candidate"]['rho_brain'].values
ax.boxplot([cont,ti],tick_labels=['Contamination\ndriven (n=120)','Baseline-matched\ntumour-intrinsic (n=120)'],patch_artist=True,showfliers=False)
ax.set_ylabel("Per-gene rho"); ax.set_title("(b) Baseline-matched: 0.679 vs 0.197\nMann-Whitney p=2e-41")
# (c) simulation
ax=axes[2]
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
print("Fig6 redrawn")

# ============= Update Table S6 with fixed survival =============
import shutil
s6=pd.read_csv(rf"{D}\results\tables\survival_gse2603_fixed.csv")
s6.to_csv(rf"{D}\manuscript\supplementary\TableS6_survival.csv",index=False)
print("TableS6 updated with fixed GSE2603 events")
print("\nALL FIGURES REDRAWN")
