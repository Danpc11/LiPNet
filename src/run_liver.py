"""Liver campaign: analytic generation model (Part A) and explicit canopy graph (Part B)."""
import numpy as np, pandas as pd, warnings, time, json, sys
import liver_canopy as lc
warnings.filterwarnings('ignore')
OUT = sys.argv[1] if len(sys.argv) > 1 else 'results'
import os; os.makedirs(OUT, exist_ok=True)
RHO = 3 ** (-1 / 3)
B_LIST = (1.0, 0.75, 2 / 3, 0.5)

# ------------------------------------------------------------------ Part A
# Lorente Table 1 (Debbaut): (d ratio, sd, L ratio, sd, splitting, sd, #generations)
TAB = dict(HA=(0.74, 0.10, 0.66, 0.14, 2.76, 1.01, 11), PV=(0.70, 0.10, 0.72, 0.18, 2.80, 0.61, 11),
           HV=(0.59, 0.12, 0.64, 0.21, 3.22, 1.46, 8))
rowsA = []
for name, (rd, sd, rl, sl, n, sn, ng) in TAB.items():
    lo, med, hi = lc.infer_b_mc(rd, sd / np.sqrt(ng), rl, sl / np.sqrt(ng), n, sn / np.sqrt(ng))
    lo3, med3, hi3 = lc.infer_b_mc(rd, sd / np.sqrt(ng), rl, sl / np.sqrt(ng), 3.0, 0.0)
    rowsA.append(dict(tree=name, rho_d=rd, rho_L=rl, n=n, b_hat=lc.infer_b(rd, rl, n), b_lo=lo, b_hi=hi,
                      b_hat_n3=lc.infer_b(rd, rl, 3.0), b_lo_n3=lo3, b_hi_n3=hi3,
                      rho_d_pred_b1=lc.optimal_diameter_ratio(1.0, rl, n),
                      rho_d_pred_b075=lc.optimal_diameter_ratio(0.75, rl, n)))
dfA = pd.DataFrame(rowsA); dfA.to_csv(f'{OUT}/partA_infer_b.tsv', sep='\t', index=False, float_format='%.3f')
print(dfA.to_string(index=False))

# efficiency of the measured generation geometry against its own-budget optimum,
# for a symmetric canopy: g=10 inlet (PV) + lobule layer + g=8 outlet (HV) (Debbaut generations)
def canopy(rd_in, rl_in, rd_out, rl_out, g_in=10, g_out=8, n=3.0):
    inl = lc.gen_tree(g_in, rd_in, rl_in, n); out = lc.gen_tree(g_out, rd_out, rl_out, n)
    N = int(round(n ** (g_in - 1)))
    # lobule layer: one equivalent channel per lobule of length L_h ~ last inlet length, radius
    # from Murray continuity with the last inlet generation (r^3 ∝ flow per channel)
    d_last, l_last = inl['d'][-1], inl['l'][-1]
    lob = dict(N=N, d=d_last * (inl['f'][-1] * n ** (g_in - 1) / 1.0) ** 0 * (1 / 1.0), l=l_last)
    lob['d'] = d_last * ((1.0 / N) / inl['f'][-1]) ** (1 / 3)
    return lc.canopy_series_parallel(inl, out, lob)
rowsE = []
geoms = {'Lorente theory (3^-1/3)': (RHO, RHO, RHO, RHO),
         'Debbaut PV in / HV out': (0.70, 0.72, 0.59, 0.64),
         'Debbaut HA in / HV out': (0.74, 0.66, 0.59, 0.64),
         'Ma d=0.79, L=0.69': (0.79, RHO, 0.80, RHO)}
for name, (a1, a2, a3, a4) in geoms.items():
    net = canopy(a1, a2, a3, a4)
    for b in B_LIST:
        opt = lc.sp_optimum(net, b, lc.sp_cost(net, b))
        rowsE.append(dict(geometry=name, b=b, eta=lc.efficiency(net, b),
                          tau_ratio_opt=opt['tau'].max() / opt['tau'].min()))
dfE = pd.DataFrame(rowsE); dfE.to_csv(f'{OUT}/partA_efficiency.tsv', sep='\t', index=False, float_format='%.4f')
print(dfE.pivot(index='geometry', columns='b', values='eta').round(3).to_string())

# ------------------------------------------------------------------ Part B
t0 = time.time()
G = lc.build_liver_graph(R=4.5); gin, gout = lc.split_liver(G)
meta = dict(nodes=G['n'], edges=G['m'], lobules=G['n_lob'], triads=G['n_tri'], inlet_edges=int(gin.m),
            sinusoids=int((gin.kind == 'sin').sum()), outlet_edges=int(gout.m), R=4.5)
print(meta); json.dump(meta, open(f'{OUT}/partB_meta.json', 'w'))
np.savez(f'{OUT}/liver_graph.npz', xy=G['xy'], E=np.array(G['E']), kind=G['kind'], depth=G['depth'], L=G['L'],
         src=G['src'], sink=G['sink'])
SIG = (0.0, 1.0, 2.0, 4.0)
rowsB = []; states = {}
for b in B_LIST:
    alpha = b / 2
    # --- outlet tree: pure tree -> closed form is exact
    out_keep = np.arange(gout.m); out_sol = lc.tree_allocation(gout, out_keep, b)
    A_out = out_sol['D']                                     # at C=1
    for sigma in SIG:
        gin.set_sigma(sigma)
        # reference: best sinusoid-assignment tree (closed form, exact on its support)
        tree = min((lc.best_sinusoid_tree(gin, b, sigma=sigma, seed=s) for s in range(3)), key=lambda s: s['D'])
        A_tree = tree['D']
        # Lorente baseline: generation ratios 3^-1/3 in the tree, Murray sinusoids, all 6 sinusoids
        r_lor = lc.murray_radii(gin, rho=RHO); w_lor = lc.to_budget(gin, r_lor, b)
        D_lor_norm = gin.d_norm(w_lor, b)
        # local adaptation from the Lorente state on the fixed support
        ad = lc.adapt(gin, b, r_lor, sigma=sigma, max_steps=3000)
        D_ad_norm = gin.d_norm(ad['w'], b, ad['keep'])
        # whole-organ numbers at total budget 1 (inlet+outlet share)
        D_tot_opt, C_in, C_out = lc.combine_budget(A_tree, A_out, b, 1.0)
        D_tot_ad, _, _ = lc.combine_budget(D_ad_norm, A_out, b, 1.0)
        D_tot_lor, _, _ = lc.combine_budget(D_lor_norm, A_out, b, 1.0)
        sh = lc.shear_fit(gin, tree['w'], tree['keep'], b)
        row = dict(b=b, sigma=sigma, D_tree=A_tree, D_adapted=D_ad_norm, D_lorente=D_lor_norm,
                   eta_lorente=A_tree / D_lor_norm, eta_adapted=A_tree / D_ad_norm,
                   excess_adapted_pct=100 * (D_ad_norm / A_tree - 1),
                   active_adapted=ad['n_active'], sin_adapted=ad['n_sin'], sin_per_lobule=ad['n_sin'] / gin.N,
                   beta_adapted=ad['beta'], steps=ad['steps'], converged=ad['converged'],
                   active_tree=tree['n_active'], C_in_share=C_in,
                   eta_organ_lorente=D_tot_opt / D_tot_lor, eta_organ_adapted=D_tot_opt / D_tot_ad,
                   shear_b_r=sh['b_r'], shear_b_l=sh['b_l'], tau_ratio=sh['tau_ratio'])
        rowsB.append(row); states[(b, sigma)] = dict(tree=tree, adapted=ad, lorente=w_lor)
        print(f"b={b:.3f} s={sigma}: eta_lor={row['eta_lorente']:.3f} eta_ad={row['eta_adapted']:.3f} "
              f"sin/lob={row['sin_per_lobule']:.2f} beta={ad['beta']} steps={ad['steps']} conv={ad['converged']} "
              f"b_r={sh['b_r']:.3f} b_l={sh['b_l']:.3f}  [{time.time()-t0:.0f}s]", flush=True)
dfB = pd.DataFrame(rowsB); dfB.to_csv(f'{OUT}/partB_campaign.tsv', sep='\t', index=False, float_format='%.5g')
import pickle; pickle.dump({f'{b:.4f}_{s}': dict(tree_keep=v['tree']['keep'], tree_w=v['tree']['w'],
                                               ad_keep=v['adapted']['keep'], ad_w=v['adapted']['w'], lor_w=v['lorente'])
                            for (b, s), v in states.items()}, open(f'{OUT}/partB_states.pkl', 'wb'))
print('done', time.time() - t0)
