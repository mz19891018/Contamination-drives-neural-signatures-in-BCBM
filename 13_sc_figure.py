"""13_sc_figure.py - barplot of program scores by cell type."""
import pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
agg = pd.read_csv(r"D:\BCBM_Project\results\tables\singlecell_by_celltype.csv")
agg = agg.sort_values("adapt", ascending=True)
fig, ax = plt.subplots(figsize=(8,5))
y = range(len(agg))
ax.barh([i+0.2 for i in y], agg.adapt, height=0.4, color="#c51b7d", label="Adaptation (de novo)")
ax.barh([i-0.2 for i in y], agg.comp, height=0.4, color="#2166ac", label="Competence (pre-existing)")
ax.set_yticks(list(y)); ax.set_yticklabels(agg.common_cell_clusters)
ax.set_xlabel("mean log-normalized module score")
ax.set_title("Fig 5. Program localization in GSE324453 snRNA-seq (n=56,268 cells)")
ax.legend(loc="lower right")
plt.tight_layout(); plt.savefig(r"D:\BCBM_Project\results\figures\fig5_singlecell.png", dpi=130); plt.close()
print("saved fig5")
