"""11_figures.py - generate main figures."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from common import RESULTS

FIG = os.path.join(RESULTS,"figures"); os.makedirs(FIG, exist_ok=True)
de = pd.read_csv(os.path.join(RESULTS,"tables","paired_de_gse184869.csv"))
prog = pd.read_csv(os.path.join(RESULTS,"tables","program_decomposition.csv"))

# Fig: volcano
fig, ax = plt.subplots(figsize=(7,6))
sig_up = de[(de.FDR<0.05)&(de.mean_logFC_BM_vs_BP>0.3)]
sig_dn = de[(de.FDR<0.05)&(de.mean_logFC_BM_vs_BP<-0.3)]
ax.scatter(de.mean_logFC_BM_vs_BP, -np.log10(de.paired_t_p+1e-300), s=3, c="lightgrey", alpha=0.5)
ax.scatter(sig_up.mean_logFC_BM_vs_BP, -np.log10(sig_up.paired_t_p), s=4, c="#d6604d", alpha=0.6, label=f"up in BM (n={len(sig_up)})")
ax.scatter(sig_dn.mean_logFC_BM_vs_BP, -np.log10(sig_dn.paired_t_p), s=4, c="#4393c3", alpha=0.6, label=f"down in BM (n={len(sig_dn)})")
for _,r in de.sort_values("paired_t_p").head(12).iterrows():
    if pd.notna(r.symbol):
        ax.annotate(r.symbol,(r.mean_logFC_BM_vs_BP,-np.log10(r.paired_t_p)),fontsize=7)
ax.set_xlabel("paired log2 fold-change (BM vs primary), 45 patients")
ax.set_ylabel("-log10 p (paired t)")
ax.set_title("Fig 2. Paired evolutionary landscape (GSE184869, n=45 pairs)")
ax.legend()
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig2_volcano.png"), dpi=130); plt.close()

# Fig: 2x2 decomposition - x=mean_BP, y=logFC, color by program
fig, ax = plt.subplots(figsize=(7.5,6))
colors = {"Adaptation_de_novo":"#c51b7d","Competence_preexisting":"#2166ac","Intermediate":"#bdbdbd"}
for name,sub in prog.groupby("program"):
    ax.scatter(sub.mean_BP, sub.mean_logFC_BM_vs_BP, s=6, alpha=0.5, c=colors.get(name,"grey"), label=f"{name} (n={len(sub)})")
ax.axvline(prog.mean_BP.median(), ls="--", c="k", lw=0.7)
ax.axvline(prog.mean_BP.quantile(0.25), ls=":", c="k", lw=0.7)
ax.axhline(0.3, ls="--", c="k", lw=0.7)
ax.set_xlabel("mean expression in primary tumor (BP)")
ax.set_ylabel("paired log2 FC (BM vs BP)")
ax.set_title("Fig 3. Competence vs Adaptation decomposition")
ax.legend(markerscale=2)
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig3_decomposition.png"), dpi=130); plt.close()

# Fig: top adaptation & competence genes
adapt = prog[prog.program=="Adaptation_de_novo"].sort_values("mean_logFC_BM_vs_BP",ascending=False).head(15)
comp  = prog[prog.program=="Competence_preexisting"].sort_values("mean_BP",ascending=False).head(15)
fig, axes = plt.subplots(1,2, figsize=(13,6))
axes[0].barh(adapt.symbol[::-1], adapt.mean_logFC_BM_vs_BP[::-1], color="#c51b7d")
axes[0].set_title("Top de-novo ADAPTATION genes"); axes[0].set_xlabel("log2 FC BM/BP")
axes[1].barh(comp.symbol[::-1], comp.mean_logFC_BM_vs_BP[::-1], color="#2166ac")
axes[1].set_title("Top pre-existing COMPETENCE genes (high in primary)"); axes[1].set_xlabel("log2 FC BM/BP")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig3b_topgenes.png"), dpi=130); plt.close()

# Fig: organ specificity
org = pd.read_csv(os.path.join(RESULTS,"tables","organ_specificity_scores.csv"), index_col=0)
fig, axes = plt.subplots(1,2, figsize=(11,4.5), sharey=True)
for ax,(dset,d) in zip(axes, org.groupby(level=0)):
    means = d.groupby("organ")[["adapt_score","comp_score"]].mean()
    means.plot(kind="bar", ax=ax, color=["#c51b7d","#2166ac"])
    ax.set_title(dset); ax.set_xlabel(""); ax.axhline(0,color="k",lw=0.5)
axes[0].set_ylabel("mean z-score")
plt.suptitle("Fig 4. Brain specificity vs other-organ metastases")
plt.tight_layout(); plt.savefig(os.path.join(FIG,"fig4_organ.png"), dpi=130); plt.close()
print("figures saved:", os.listdir(FIG))
