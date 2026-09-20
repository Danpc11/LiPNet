#!/usr/bin/env bash
# Reproduces every table and figure. Run from the repository root.
# Stage 0 is fast (< 1 min). Stages 1–3 are the physics campaigns (≈ 1–2 h in total).
set -e
export PYTHONPATH=src
mkdir -p results && cp data/* results/          # published-series table and cached states

# ---- Stage 0: clinical analysis (the manuscript) -------------------------------------
python src/series_from_pdfs.py                    # -> results/series_pdfs_z.tsv (values transcribed from the 12 PDFs)
python - <<'PY'                                   # pooled logistic on z_P (5 SFSS groups) with bootstrap CI
import numpy as np, pandas as pd, json
from scipy.optimize import minimize
df=pd.read_csv('results/series_pdfs_z.tsv',sep='\t')
s=df[df.desenlace.str.contains('SFSS|disfunción') & df.z_P.notna()]
z,ev,n=s.z_P.values,s.eventos.values,s.n.values
def fit(z,ev,n):
    nll=lambda b: -np.sum(ev*np.log(np.clip(1/(1+np.exp(-(b[0]+b[1]*z))),1e-9,1-1e-9))+(n-ev)*np.log(np.clip(1-1/(1+np.exp(-(b[0]+b[1]*z))),1e-9,1-1e-9)))
    return minimize(nll,[-4,1.5]).x
b=fit(z,ev,n); rng=np.random.default_rng(0)
B=np.array([fit(z,rng.binomial(n,ev/n),n) for _ in range(2000)]); lo,hi=np.percentile(B,[2.5,97.5],axis=0)
zz=np.linspace(0.8,4,33); P=1/(1+np.exp(-(b[0]+b[1]*zz))); PB=1/(1+np.exp(-(B[:,[0]]+B[:,[1]]*zz))); ci=np.percentile(PB,[2.5,97.5],axis=0)
json.dump(dict(b0=float(b[0]),b1=float(b[1]),b0_ci=[float(lo[0]),float(hi[0])],b1_ci=[float(lo[1]),float(hi[1])],groups=int(len(s)),events=int(ev.sum()),patients=int(n.sum()),z=zz.tolist(),p=P.tolist(),lo=ci[0].tolist(),hi=ci[1].tolist()),open('results/logit_zP.json','w'))
print('logit =',b)
PY
python src/make_figs_jhep.py                      # -> results/JHEP_Fig1-4.pdf/.png  (needs results/partB_states.pkl and liver_*.png from data/)
python src/make_supp_figs.py                      # -> results/JHEP_FigS1-S3.pdf/.png
cp results/JHEP_Fig*.p* figures/
(cd manuscript && node make_manuscript.js && node make_supplement.js)   # -> JHEP_manuscript.docx, JHEP_supplement.docx

# ---- Stage 1: network model on the 2D liver (Figs. S-series; ≈ 25 min) ---------------
# python src/run_liver.py results                 # optimum, adaptation, fluctuating demand; writes partB_states.pkl
# python src/make_figs.py

# ---- Stage 2: scaling with size, 2D and 3D (≈ 30 min) ---------------------------------
# python src/scaling.py results                   # closed-form exponents and allometric closure
# for R in 3.0 4.5 6.5 9.0 12.0; do python src/graph_scaling.py $R; done
# for R in 3.0 4.0 5.0 6.5 8.5; do python src/liver3d.py $R results; done
# python src/scaling_figs.py

# ---- Stage 3: controls (≈ 1 min) --------------------------------------------------------
# python src/certify.py                           # exhaustive optimum on 7 lobules vs heuristic
