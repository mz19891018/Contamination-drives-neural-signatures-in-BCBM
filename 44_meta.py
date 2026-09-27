# 44_meta.py - Fisher z fixed-effect meta of BM-internal rhos
import numpy as np
from scipy.stats import norm
studies=[
 ("GSE184869",69,0.93),
 ("GSE125989",16,0.882),
 ("GSE14017",15,0.868),
 # GSE14018 n=7 excluded (n<10)
]
zs=[]; ws=[]
for name,n,r in studies:
    z=0.5*np.log((1+r)/(1-r)); w=n-3
    zs.append(z*w); ws.append(w)
    print(f"{name}: n={n} r={r} z={z:.3f}")
zbar=sum(zs)/sum(ws)
rbar=np.tanh(zbar)
se=1/np.sqrt(sum(ws))
lo=np.tanh(zbar-1.96*se); hi=np.tanh(zbar+1.96*se)
p=2*(1-norm.cdf(abs(zbar/se)))
print(f"\nFisher-z fixed-effect meta (k={len(studies)}, total n={sum(s[1] for s in studies)}):")
print(f"  pooled rho = {rbar:.3f}, 95% CI [{lo:.3f}, {hi:.3f}], p={p:.2e}")
