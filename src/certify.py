"""Exact certificate on small domains: enumerate every spanning tree of the inlet problem
(inlet tree fixed + one feeding sinusoid per lobule -> 6^N candidates) and compare with the
coordinate-descent optimum. Also: dispersion of the heuristic vs N."""
import numpy as np, itertools, time, pandas as pd, warnings; warnings.filterwarnings('ignore')
import liver_canopy as lc, liver3d as l3
rows=[]
for R in (1.2,):                      # 7 and 13-19 lobules in 2D
    G = lc.build_liver_graph(R=R); gin,_ = lc.split_liver(G)
    for b in (1.0, 0.75, 0.5):
        ft = l3.FastTree(gin, b); cs = [ft.by_centre[c] for c in ft.centres]
        t=time.time(); best=np.inf; n=0
        tree=ft.tree
        for combo in itertools.product(*cs):
            k=np.zeros(gin.m)
            for e in combo: k[ft.path[gin.E[e][0]]]+=1; k[e]=1
            act=np.flatnonzero(k>0); S=ft.term(act,k[act]).sum()
            if S<best: best=S
            n+=1
        Dex=best**(1+1/ft.alpha)
        heur=[ft.solve(seed=s)['D'] for s in range(10)]
        rows.append(dict(R=R,N=gin.N,b=b,candidates=n,D_exact=Dex,D_heuristic_best=min(heur),
                         rel_gap=min(heur)/Dex-1,frac_seeds_at_opt=np.mean(np.abs(np.array(heur)/Dex-1)<1e-9),t=time.time()-t))
        print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv('results/certificate_small.tsv',sep='\t',index=False)
