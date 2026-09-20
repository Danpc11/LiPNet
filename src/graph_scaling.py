import numpy as np, liver_canopy as lc, warnings, time, sys, pandas as pd; warnings.filterwarnings('ignore')
R = float(sys.argv[1]); rows=[]
G = lc.build_liver_graph(R=R); gin, gout = lc.split_liver(G)
for b in (1.0, 0.75, 0.5):
    t=time.time(); tree = lc.best_sinusoid_tree(gin, b, seed=0, sweeps=8)
    out = lc.tree_allocation(gout, np.arange(gout.m), b)
    rows.append(dict(R=R, N=gin.N, b=b, D_in=tree['D'], D_out=out['D'], D_in_perF2=tree['D']/gin.N**2, t=time.time()-t))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f'results/graph_scaling_R{R}.tsv', sep='\t', index=False)
