# 64_fix_figures.py - fix Fig4b ratio, Fig4c labels, Fig5 x-axis
import pandas as pd, numpy as np, re, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
D=r"D:\BCBM_Project"; FIG=rf"{D}\manuscript\figures"

th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
def kind(s):
    if re.match(r"^BM\d",s): return "BM_named"
    if re.match(r"^\d+M_RCS",s): return "BM_RCS"
    if re.match(r"^MAYO",s): return "BM_MAYO"
    if re.match(r"^BP\d",s): return "P_named"
    if re.match(r"^\d+P_RCS",s): return "P_RCS"
    return "Other"
th["batch"]=[kind(s) for s in th.index]
th["brain"]=th[["Neuron","Astrocyte","Oligodendrocyte"]].sum(1)
bm=th[th.batch.str.startswith("BM")]; prim=th[th.batch.str.startswith("P")]

# ===== Fig4b: use MEAN, label clearly =====
fig,ax=plt.subplots(1,1,figsize=(6,5))
x=np.arange(3); w=0.35
bm_m=[bm['Neuron'].mean(),bm['Astrocyte'].mean(),bm['Oligodendrocyte'].mean()]
p_m=[prim['Neuron'].mean(),prim['Astrocyte'].mean(),prim['Oligodendrocyte'].mean()]
ax.bar(x-w/2,bm_m,w,label='BrainMet (mean)',color='#C44E52',alpha=0.7)
ax.bar(x+w/2,p_m,w,label='Primary (mean)',color='#4C72B0',alpha=0.7)
ax.set_xticks(x); ax.set_xticklabels(['Neuron','Astrocyte','Oligodendrocyte'])
ax.set_ylabel("Mean fraction"); ax.set_title("Brain content by class: BM vs Primary")
ax.legend(fontsize=9)
for i,(b,p) in enumerate(zip(bm_m,p_m)):
    ax.text(i,max(b,p)+0.004,f"{b/p:.1f}x\n(mean)",ha='center',fontsize=9,fontweight='bold')
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure4b_panel.png",dpi=200,bbox_inches='tight')
plt.close()
print(f"Fig4b: neuron mean ratio={bm_m[0]/p_m[0]:.2f}, astro={bm_m[1]/p_m[1]:.2f}, odg={bm_m[2]/p_m[2]:.2f}")
print(f"  neuron median: BM={bm['Neuron'].median():.4f}, P={prim['Neuron'].median():.4f}, ratio={bm['Neuron'].median()/prim['Neuron'].median():.2f}")

# ===== Fig4c: fix label overlap =====
J4v2=pd.read_csv(rf"{D}\results\tables\J4v2_neuron_only.csv")
fig,ax=plt.subplots(1,1,figsize=(6,5))
ax.scatter(J4v2['rho_brain'],J4v2['rho_neuron'],s=8,alpha=0.4,c='#4C72B0')
ax.plot([-0.5,1],[-0.5,1],'r--',lw=1)
# labels with manual offsets
labels={'ATP1A2':(0.75,0.69,10,-15),'TUBB4A':(0.74,0.71,-50,10),
        'C1QL1':(0.72,0.59,10,10),'ABCG2':(0.11,0.01,-40,-15),'GRB7':(0.13,0.28,10,8)}
for g,(x0,y0,dx,dy) in labels.items():
    row=J4v2[J4v2.gene==g]
    if len(row): ax.annotate(g,(row.rho_brain.values[0],row.rho_neuron.values[0]),
        fontsize=9,fontweight='bold',xytext=(dx,dy),textcoords='offset points',
        arrowprops=dict(arrowstyle='-',color='gray',lw=0.5))
ax.set_xlabel("rho (all brain classes)"); ax.set_ylabel("rho (neuron only)")
ax.set_title("Neuron-only definition: kappa=0.700, rho corr=0.944")
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure4c_panel.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig4c labels fixed")

# ===== Fig5: x-axis capped at 1.0 =====
fig,ax=plt.subplots(1,1,figsize=(8,4))
cohorts=['GSE184869','GSE125989','GSE14017','Pooled (FE)']
rhos=[0.93,0.882,0.868,0.918]
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
ax.set_xlim(0.5,1.02); ax.set_xlabel("Spearman rho (within BM samples only)")
ax.set_title("Replication across independent BCBM cohorts")
ax.axvline(0,ls='--',c='gray',lw=0.8)
plt.tight_layout()
plt.savefig(rf"{FIG}\Figure5_forest.png",dpi=200,bbox_inches='tight')
plt.close()
print("Fig5 x-axis capped at 1.0")

# Now rebuild full Figure 4 with fixed panels
import openpyxl
adapt=[g.strip() for g in open(rf"{D}\results\tables\geneset_adaptation.txt",encoding="utf-8") if g.strip()]
pm=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv"); emap=dict(zip(pm.ensembl,pm.symbol))
wb=openpyxl.load_workbook(rf"{D}\data\raw\GSE184869_expression.xlsx",read_only=True)
ws=wb["log2TMMCPM"]; rows=ws.iter_rows(values_only=True); header=list(next(rows))
bm_idx=[i for i,h in enumerate(header) if h and kind(str(h)).startswith("BM")]
bm_names=[header[i] for i in bm_idx]
expr={}
for r in rows:
    sym=emap.get(r[0])
    if sym in set(adapt): expr[sym]=[r[i] for i in bm_idx]
wb.close()
EA=pd.DataFrame(expr,index=bm_names).T
common=[s for s in bm_names if s in th.index]
EA=EA[common]; bf=th.loc[common,"brain"]; bat=th.loc[common,"batch"]
adapt_score=EA.mean(0)
LMO=pd.read_csv(rf"{D}\results\tables\leave_marker_stepwise.csv")
K=pd.read_csv(rf"{D}\results\tables\moduleK_null.csv")

fig,axes=plt.subplots(2,2,figsize=(13,10))
# (a) batch
ax=axes[0,0]
for b,c in [('BM_named','#4C72B0'),('BM_RCS','#55A868'),('BM_MAYO','#C44E52')]:
    mask=bat==b
    ax.scatter(bf[mask],adapt_score[mask],s=30,c=c,alpha=0.7,label=f"{b} (n={mask.sum()})")
r_all=np.corrcoef(bf.rank(),adapt_score.rank())[0,1]
ax.set_xlabel("Brain-cell content (all 3 classes)"); ax.set_ylabel("Adaptation module score")
ax.set_title(f"(a) Batch: all rho={r_all:.3f}; MAYO alone rho=0.845")
ax.legend(fontsize=8)
# (b) fixed
ax=axes[0,1]
ax.bar(x-w/2,bm_m,w,label='BrainMet',color='#C44E52',alpha=0.7)
ax.bar(x+w/2,p_m,w,label='Primary',color='#4C72B0',alpha=0.7)
ax.set_xticks(x); ax.set_xticklabels(['Neuron','Astrocyte','Oligodendrocyte'])
ax.set_ylabel("Mean fraction"); ax.set_title("(b) Brain content by class (mean)")
ax.legend(fontsize=8)
for i,(b,p) in enumerate(zip(bm_m,p_m)):
    ax.text(i,max(b,p)+0.004,f"{b/p:.1f}x",ha='center',fontsize=9,fontweight='bold')
# (c) fixed labels
ax=axes[1,0]
ax.scatter(J4v2['rho_brain'],J4v2['rho_neuron'],s=8,alpha=0.4,c='#4C72B0')
ax.plot([-0.5,1],[-0.5,1],'r--',lw=1)
for g,(x0,y0,dx,dy) in labels.items():
    row=J4v2[J4v2.gene==g]
    if len(row): ax.annotate(g,(row.rho_brain.values[0],row.rho_neuron.values[0]),
        fontsize=9,fontweight='bold',xytext=(dx,dy),textcoords='offset points',
        arrowprops=dict(arrowstyle='-',color='gray',lw=0.5))
ax.set_xlabel("rho (all brain classes)"); ax.set_ylabel("rho (neuron only)")
ax.set_title("(c) Neuron-only: kappa=0.700, rho corr=0.944")
# (d)
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
print("Full Figure 4 rebuilt")
print("\nALL FIGURE FIXES DONE")
