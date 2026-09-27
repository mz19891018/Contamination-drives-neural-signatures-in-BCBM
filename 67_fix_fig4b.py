import pandas as pd, numpy as np, re
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
D=r"D:\BCBM_Project"

th=pd.read_csv(rf"{D}\results\tables\bp_combined_fractions.csv",index_col=0)
def kind(s):
    if re.match(r"^BM\d",s): return "BM"
    if re.match(r"^\d+M_RCS",s): return "BM"
    if re.match(r"^MAYO",s): return "BM"
    if re.match(r"^BP\d",s): return "P"
    if re.match(r"^\d+P_RCS",s): return "P"
    return "Other"
th["batch"]=[kind(s) for s in th.index]
bm=th[th.batch=="BM"]; prim=th[th.batch=="P"]

fig,ax=plt.subplots(1,1,figsize=(6,5))
x=np.arange(3); w=0.35
bm_m=[bm['Neuron'].mean(),bm['Astrocyte'].mean(),bm['Oligodendrocyte'].mean()]
p_m=[prim['Neuron'].mean(),prim['Astrocyte'].mean(),prim['Oligodendrocyte'].mean()]
ax.bar(x-w/2,bm_m,w,label='BrainMet',color='#C44E52',alpha=0.7)
ax.bar(x+w/2,p_m,w,label='Primary',color='#4C72B0',alpha=0.7)
ax.set_xticks(x); ax.set_xticklabels(['Neuron','Astrocyte','Oligodendrocyte'])
ax.set_ylabel("Mean fraction")
ax.set_title("Brain content by class: BM vs Primary",pad=15)
ax.legend(fontsize=10,loc='upper left')
for i,(b,p) in enumerate(zip(bm_m,p_m)):
    ax.text(i,max(b,p)+0.004,f"{b/p:.1f}x",ha='center',fontsize=11,fontweight='bold')
ax.set_ylim(0,max(bm_m)*1.2)
plt.tight_layout()
plt.savefig(rf"{D}\manuscript\figures\Figure4b_panel.png",dpi=200,bbox_inches='tight')
plt.close()
print("Figure4b_panel.png re-exported")
print(f"Ratios: neuron={bm_m[0]/p_m[0]:.2f}, astro={bm_m[1]/p_m[1]:.2f}, odg={bm_m[2]/p_m[2]:.2f}")
