# 38_Anew_metrics.py
import pandas as pd, numpy as np
D=r"D:\BCBM_Project"
# read sim full logFC (R rds -> use pyreadr? not available. Use feather? Simpler: R exported rds.
# We'll re-export via R to csv. Do that in a tiny R call.)
import subprocess
subprocess.run([r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe","-e",
  'x<-readRDS("D:/BCBM_Project/results/tables/sim_full_logFC.rds"); write.csv(x,"D:/BCBM_Project/results/tables/sim_full_logFC.csv")'],
  capture_output=True)
sim=pd.read_csv(rf"{D}\results\tables\sim_full_logFC.csv",index_col=0)
print("sim shape:",sim.shape)

prog=pd.read_csv(rf"{D}\results\tables\program_decomposition.csv")
obs=dict(zip(prog.symbol, prog.mean_logFC_BM_vs_BP))
up=prog[(prog.program!="Intermediate")&(prog.mean_logFC_BM_vs_BP>0.3)]
up_genes=[g for g in up.symbol.tolist() if g in sim.index]
print("up genes in sim:",len(up_genes))

fs=[float(c) for c in sim.columns]
o=np.array([obs[g] for g in up_genes])
S=sim.loc[up_genes].values  # genes x f

# Metric 1: RMSE
rmse=np.sqrt(np.mean((S-o[:,None])**2,axis=0))
fstar_rmse=fs[np.argmin(rmse)]
print("\n=== RMSE per f ===")
for f,r in zip(fs,rmse): print(f"  f={f:.3f} RMSE={r:.3f}")
print(f"f* (RMSE min) = {fstar_rmse}")

# Metric 2: Deming slope (approx via total least squares)
def deming_slope(x,y):
    xm,ym=x.mean(),y.mean()
    sxx=np.sum((x-xm)**2); syy=np.sum((y-ym)**2); sxy=np.sum((x-xm)*(y-ym))
    # assume equal variance -> slope = (syy-sxx+sqrt((syy-sxx)**2+4*sxy**2))/(2*sxy)
    return (syy-sxx+np.sqrt((syy-sxx)**2+4*sxy**2))/(2*sxy+1e-12)
print("\n=== Deming slope per f ===")
slopes=[]
for j,f in enumerate(fs):
    s=deming_slope(S[:,j],o)
    slopes.append(s)
    print(f"  f={f:.3f} slope={s:.3f}")
# f where slope=1
fstar_slope=np.interp(1.0, slopes, fs)  # slopes increasing with f? check
print(f"f* (slope=1) ~ {fstar_slope:.3f}")

# Metric 3: per-gene inverse f for marker genes
markers=["TUBB4A","ATP1A2","MAG","APLP1","PAX6","C1QL1","CHGB","ZIC2"]
print("\n=== per-gene f_g (logFC_sim=logFC_obs) ===")
fg=[]
for g in markers:
    if g not in sim.index: continue
    target=obs.get(g,np.nan)
    if np.isnan(target): continue
    curve=sim.loc[g].values
    # interpolate
    if target>curve.max():
        print(f"  {g}: obs={target:.2f} > max sim={curve.max():.2f} -> f>{fs[-1]}")
        fg.append(np.nan)
    else:
        f_g=np.interp(target, curve, fs)
        fg.append(f_g)
        print(f"  {g}: obs={target:.2f} -> f_g={f_g:.3f}")
fg_arr=np.array(fg)
print(f"median f_g = {np.nanmedian(fg_arr):.3f}, IQR = {np.nanpercentile(fg_arr,25):.3f}-{np.nanpercentile(fg_arr,75):.3f}")

# Stratified CAF at f=0.22 (combined-ref BM brain fraction)
f_use=0.22; j_use=fs.index(f_use) if f_use in fs else np.argmin(np.abs(np.array(fs)-f_use))
print(f"\n=== stratified CAF at f={fs[j_use]} ===")
for lo,hi,label in [(3,100,"strong >3"),(1,3,"medium 1-3"),(0.3,1,"weak 0.3-1")]:
    mask=(o>=lo)&(o<hi)
    caf=np.clip(S[mask,j_use]/o[mask],0,1)
    print(f"  {label}: n={mask.sum()} median CAF={np.median(caf):.3f}")

# Top-N recovery at f=0.22: need simulated significance; approximate by sim logFC>0.3
top200=up.sort_values("mean_logFC_BM_vs_BP",ascending=False).head(200).symbol.tolist()
top200=[g for g in top200 if g in sim.index]
rec200=np.mean(sim.loc[top200].iloc[:,j_use]>0.3)
print(f"\nTop-200 recovery (sim logFC>0.3) at f={fs[j_use]}: {rec200:.1%}")

res=pd.DataFrame({"f":fs,"RMSE":rmse,"Deming_slope":slopes})
res.to_csv(rf"{D}\results\tables\Anew_metrics.csv",index=False)
print("\nsaved Anew_metrics.csv")
