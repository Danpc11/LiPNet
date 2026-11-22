"""The network model applied to the human liver, with nothing added to it.

Everything here is a consequence of the cost-constrained optimum of src/model/scaling.py and network.py. It answers
five questions that the manuscript asserts or leaves open, and that decide how far the model reaches into a liver:

1. What is the wall shear stress of a vessel?  On every edge of the optimum D_e / c_e = tau0^2 (single exponent b),
   but the Poiseuille wall shear is tau_e = f_e / r_e^3, and the two differ unless b = 1:

       tau_e = tau0 * m_e^((b - 1) / 2),      m_e = r_e^2 l_e the vascular mass of the edge.

   tau0 is the scale of the design, not the shear of any vessel. `canopy` returns the shear of every generation.

2. Which exponent b does a human liver have?  The optimum fixes the diameter ratio between successive generations,

       r_{i+1} / r_i = n^(-1/(2+b)) * lambda^((1-b) / (2 (2+b))),

   with n the branching number and lambda the length ratio. `exponent_from_cast` inverts it for the casts already
   quoted in network.py (Debbaut et al.). The portal vein and the hepatic vein do not give the same b. A model with
   one exponent per tree is not closed: the costs m^b of the two trees then have different units, so their joint
   optimum depends on one more number (the budget of one tree relative to the other) that the model does not
   supply. The exponent of the hepatic vein is reported and not used: the same radii give, by Poiseuille, a drop
   larger than the gradient of the whole liver (src/model/human.py), so they are not the calibres of a perfused
   liver. Results are reported across b, and where they depend on the split between the trees, over every split.

3. How is the portocaval pressure drop shared between the portal tree, the sinusoids and the hepatic venous tree
   at the size of a human liver (about 10^6 lobules), and how much of it sits in the trunks?  `canopy`.

4. What does the model predict for a partial graft?  A graft is the sub-network below one first-order branch. Its
   vessels keep their conductance, so with z_F times the donor flow in every lobule the pressure drop across it is
   z_F times the donor drop minus the trunks it left behind:

       z_R(graft) = 1 - (share of the donor pressure drop in the generations removed).

   `lobar_graft` returns it, and `graft_on_graph` repeats the calculation on the explicit 3D network, where the
   graft keeps the hepatic venous edges that drain it at the calibre they had in the whole liver.

5. What does the second inflow change?  The hepatic artery joins the portal blood at the sinusoid, so the portal
   tree carries portal flow and the sinusoids and hepatic veins carry portal plus arterial flow. The clinical flow
   load is a portal load. With a the arterial-to-portal flow ratio of the donor and phi the arterial flow per lobule
   of the graft relative to the donor,

       z_P = s_in z_F + (1 - s_in) (z_F + phi a) / (1 + a),

   where s_in is the share of the portal tree in the donor pressure drop. `two_inflows` evaluates it.

    python src/model/liver.py      -> results/model_liver.json
"""
import json, os, sys
import numpy as np
from scipy.optimize import brentq

# Human hepatic casts as quoted in src/model/network.py and plots.shear_profile (Debbaut et al.):
# diameter ratio, length ratio and branching number between successive generations.
CASTS = {'portal vein': dict(d_ratio=0.70, l_ratio=0.72, n=2.80),
         'hepatic vein': dict(d_ratio=0.59, l_ratio=0.64, n=3.22),
         'hepatic artery': dict(d_ratio=0.74, l_ratio=0.66, n=2.76)}
HUMAN_GENERATIONS = 14                      # 3^13 = 1.6e6 lobules, the order of a human liver
ARTERIAL_FRACTION = 0.25                    # hepatic artery share of total hepatic inflow in the donor
B_CANOPY = (1.0, 0.75, 0.5)                 # maintenance exponents at which the symmetric canopy enters the rigid range, b in [0.5, 1]


def diameter_ratio(b, n=3.0, l_ratio=3.0 ** (-1 / 3)):
    """r_{i+1} / r_i at the optimum, from w_e ~ (f_e^2 / a_e)^(1/(1+alpha)) with w = r^4 / l and a = l^(3b/2)."""
    return n ** (-1.0 / (2 + b)) * l_ratio ** ((1 - b) / (2 * (2 + b)))


def exponent_from_cast(d_ratio, l_ratio, n):
    """The maintenance exponent b whose optimum reproduces a measured diameter ratio."""
    return float(brentq(lambda b: diameter_ratio(b, n, l_ratio) - d_ratio, 0.05, 3.0))


def canopy(g, b, d=3, F=1.0, C=1.0, Lh=1.0, n=3.0):
    """The optimum of the symmetric canopy, generation by generation (same optimum as scaling.canopy_optimum).

    w_e = lam (f_e^2 / a_e)^(1/(1+alpha)) with lam fixed by the budget sum_e a_e w_e^alpha = C.
    """
    alpha = b / 2; i = np.arange(g); cnt = n ** i; f = F / cnt; l = Lh * n ** ((g - 1 - i) / d); N = n ** (g - 1)
    cnt = np.concatenate([cnt, [N], cnt[::-1]]); f = np.concatenate([f, [F / N], f[::-1]]); l = np.concatenate([l, [Lh], l[::-1]])
    a = l ** (1.5 * b); u = (f ** 2 / a) ** (1 / (1 + alpha)); w = u * (C / np.sum(cnt * a * u ** alpha)) ** (1 / alpha)
    r = (w * l) ** 0.25; m = r ** 2 * l; dP = f / w; tau = f / r ** 3; k = g; tot = dP.sum()
    return dict(N=float(N), D=float(np.sum(cnt * f ** 2 / w)), dP_total=float(tot), cnt=cnt, f=f, l=l, w=w, r=r, m=m, dP=dP, tau=tau, k=k,
                V=float(np.sum(cnt * m)),
                share_inlet=float(dP[:k].sum() / tot), share_sinusoid=float(dP[k] / tot), share_outlet=float(dP[k + 1:].sum() / tot),
                share_trunks=float((dP[0] + dP[-1]) / tot), trunk_within_tree=float(dP[0] / dP[:k].sum()),
                tau_sinusoid=float(tau[k]), tau_trunk=float(tau[0]))


def lobar_graft(g=HUMAN_GENERATIONS, b=1.0, level=1, **kw):
    """z_R of the sub-network below one branch of generation `level` (1 = a lobe), relative to the whole liver."""
    c = canopy(g, b, **kw); dP = c['dP']
    removed = (dP[:level].sum() + dP[len(dP) - level:].sum()) / dP.sum()
    return dict(level=level, fraction_of_lobules=float(3.0 ** -level), zR=float(1 - removed), removed_share=float(removed),
                share_inlet=c['share_inlet'], share_sinusoid=c['share_sinusoid'], share_outlet=c['share_outlet'])


def two_inflows(zF, share_inlet, zR_network=1.0, arterial_fraction=ARTERIAL_FRACTION, arterial_load=1.0):
    """Pressure and resistance loads when the flow load is the portal one and the artery joins at the sinusoid.

    arterial_load is the arterial flow per lobule of the graft relative to the donor (phi): 0 if the artery
    contributes nothing, 1 if each lobule keeps its donor arterial flow, zF if the artery scales with the portal vein.
    zR_network multiplies the whole drop (the trunks a graft leaves behind).
    """
    a = arterial_fraction / (1 - arterial_fraction)
    zP = zR_network * (share_inlet * zF + (1 - share_inlet) * (zF + arterial_load * a) / (1 + a))
    return dict(zF=float(zF), zP=float(zP), zR=float(zP / zF), arterial_load=float(arterial_load))


def allometry(b, flow_exp=0.77, lobule_exp=0.10, n_exp=0.59, gs=range(8, 19)):
    """Closure across body mass M at invariant portocaval pressure drop, solved numerically.

    N ~ M^n_exp, F ~ M^flow_exp, lobule size ~ M^lobule_exp; the budget C(M) is whatever keeps D/F constant.
    Returns the exponents with M of the physical wall shear in the sinusoid and in the trunk, of the vascular
    mass, and of the scale tau0 = sqrt(D/C) the manuscript reports.
    """
    gs = list(gs); ref = 3.0 ** (gs[0] - 1); rows = []; target = None
    for g in gs:
        M = (3.0 ** (g - 1) / ref) ** (1 / n_exp); F, Lh = M ** flow_exp, M ** lobule_exp
        drop = lambda logC: np.log(canopy(g, b, F=F, C=np.exp(logC), Lh=Lh)['D'] / F)
        if target is None:
            logC, target = 0.0, drop(0.0)
        else:
            lo, hi = -5.0, 5.0
            while drop(lo) < target: lo -= 20
            while drop(hi) > target: hi += 20
            logC = brentq(lambda x: drop(x) - target, lo, hi, xtol=1e-12)
        c = canopy(g, b, F=F, C=np.exp(logC), Lh=Lh)
        rows.append((M, c['tau_sinusoid'], c['tau_trunk'], c['V'], np.sqrt(c['D'] / np.exp(logC))))
    A = np.log(np.array(rows)); ex = lambda j: float(np.polyfit(A[:, 0], A[:, j], 1)[0])
    return dict(b=float(b), flow_exp=flow_exp, lobule_exp=lobule_exp, shear_sinusoid=ex(1), shear_trunk=ex(2), vascular_mass=ex(3), tau0=ex(4))


def graft_on_graph(R=5.0, b=1.0, seed=0, starts=3):
    """The same graft on the explicit 3D network: one first-order portal branch, the sinusoids it feeds and the
    hepatic venous edges that drain them, every edge at the conductance it has in the whole liver."""
    import networkx as nx
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    from model import network as lc, network3d as l3
    G = l3.build_liver3d(R=R, seed=seed); gin, gout = lc.split_liver(G)
    ft = l3.FastTree(gin, b); best = min((ft.solve(seed=s) for s in range(starts)), key=lambda s: s['D'])
    keep = np.asarray(best['keep']); alpha = b / 2
    Tin = nx.Graph(); Tin.add_edges_from((gin.E[k][0], gin.E[k][1], {'k': int(k)}) for k in keep)
    Tout = nx.Graph(); Tout.add_edges_from((gout.E[k][0], gout.E[k][1], {'k': int(k)}) for k in range(gout.m))
    pin, pout = nx.shortest_path(Tin, gin.src), nx.shortest_path(Tout, gout.src)
    cen = list(gin.dem); N = len(cen)
    e_in = {c: [Tin[u][v]['k'] for u, v in zip(pin[c][:-1], pin[c][1:])] for c in cen}
    e_out = {c: [Tout[u][v]['k'] for u, v in zip(pout[c][:-1], pout[c][1:])] for c in cen}

    def counts(sel):
        ki, ko = np.zeros(gin.m), np.zeros(gout.m)
        for c in sel: ki[e_in[c]] += 1; ko[e_out[c]] += 1
        return ki, ko
    ki, ko = counts(cen); mi, mo = ki > 0, ko > 0; sin = gin.kind == 'sin'
    wi = np.zeros(gin.m); wi[mi] = (ki[mi] ** 2 / gin.a(b)[mi]) ** (1 / (1 + alpha))
    wo = np.zeros(gout.m); wo[mo] = (ko[mo] ** 2 / gout.a(b)[mo]) ** (1 / (1 + alpha))
    diss = lambda a_, b_: np.array([np.sum(a_[mi & ~sin] ** 2 / wi[mi & ~sin]), np.sum(a_[mi & sin] ** 2 / wi[mi & sin]), np.sum(b_[mo] ** 2 / wo[mo])])
    Dd = diss(ki, ko); grafts = []
    for br in sorted({e_in[c][0] for c in cen}):
        sel = [c for c in cen if e_in[c][0] == br]; kgi, kgo = counts(sel); Dg = diss(kgi, kgo)
        rho = (Dg / len(sel)) / (Dd / N)
        grafts.append(dict(fraction_of_lobules=len(sel) / N, zR=float((Dg.sum() / len(sel)) / (Dd.sum() / N)),
                           inlet=float(rho[0]), sinusoid=float(rho[1]), outlet=float(rho[2])))
    return dict(R=float(R), b=float(b), seed=int(seed), n_lobules=N, share_inlet=float(Dd[0] / Dd.sum()),
                share_sinusoid=float(Dd[1] / Dd.sum()), share_outlet=float(Dd[2] / Dd.sum()), grafts=grafts)


def summary(graphs=True):
    out = {}
    out['casts'] = {k: dict(v, b=exponent_from_cast(**v), model_ratio_b_half=float(diameter_ratio(0.5, v['n'], v['l_ratio'])),
                            model_ratio_b_one=float(diameter_ratio(1.0, v['n'], v['l_ratio'])),
                            shear_per_generation=float((1 / v['n']) / v['d_ratio'] ** 3)) for k, v in CASTS.items()}
    b_pv = out['casts']['portal vein']['b']
    settings = [('b = 1', 1.0), ('portal vein cast', b_pv), ('b = 0.75', 0.75), ('b = 2/3', 2 / 3), ('b = 0.5', 0.5)]
    out['human_liver'] = []
    for name, b in settings:
        c = canopy(HUMAN_GENERATIONS, b); tau0 = np.sqrt(c['D'])
        gs = np.arange(HUMAN_GENERATIONS - 2, HUMAN_GENERATIONS + 2); R = [canopy(int(g), b) for g in gs]
        out['human_liver'].append(dict(
            setting=name, b=float(b), n_lobules=c['N'], share_inlet=c['share_inlet'], share_sinusoid=c['share_sinusoid'],
            share_outlet=c['share_outlet'], share_trunks=c['share_trunks'], trunk_within_tree=c['trunk_within_tree'],
            zR_lobar_graft=lobar_graft(HUMAN_GENERATIONS, b, 1)['zR'], zR_sectoral_graft=lobar_graft(HUMAN_GENERATIONS, b, 2)['zR'],
            local_exponent_D_N=float(np.polyfit(np.log([r['N'] for r in R]), np.log([r['D'] for r in R]), 1)[0]),
            asymptotic_exponent=float((2 - b) / b), sinusoid_shear_over_tau0=float(c['tau_sinusoid'] / tau0),
            trunk_shear_over_tau0=float(c['tau_trunk'] / tau0), sinusoid_over_trunk_shear=float(c['tau_sinusoid'] / c['tau_trunk']),
            shear_identity_holds=bool(np.allclose(c['tau'], tau0 * c['m'] ** ((b - 1) / 2)))))
    out['two_inflows'] = []
    for name, b in settings:
        lob = lobar_graft(HUMAN_GENERATIONS, b, 1)
        for zF in (3.3, 4.1):
            for phi in (0.0, 1.0, zF):
                out['two_inflows'].append(dict(setting=name, share_inlet=lob['share_inlet'], zF=zF, arterial_load=float(phi),
                                               zR=two_inflows(zF, lob['share_inlet'], lob['zR'], arterial_load=phi)['zR']))
    out['two_inflows_any_split'] = [dict(zF=zF, arterial_load=float(phi), zR_interval=[two_inflows(zF, 0.0, arterial_load=phi)['zR'], 1.0])
                                    for zF in (3.3, 4.1) for phi in (0.0, 1.0)]
    out['allometry'] = []
    for flow_exp, flab in ((0.77, 'portal flow'), (0.88, 'total hepatic flow')):
        for lob_exp, n_exp in ((0.0, 0.89), (0.10, 0.59)):
            for b in (1.0, 0.9, 0.85, 0.8, 0.75, 2 / 3, 0.5):
                r = allometry(b, flow_exp, lob_exp, n_exp); r.update(flow=flab, inside_observed_blood_volume=bool(0.80 <= r['vascular_mass'] <= 0.92)); out['allometry'].append(r)
    rng = lambda k: [float(min(r[k] for r in out['allometry'])), float(max(r[k] for r in out['allometry']))]
    out['allometry_ranges'] = dict({k: rng(k) for k in ('tau0', 'shear_sinusoid', 'shear_trunk', 'vascular_mass')},
                                   observed_hepatic_blood_volume=[0.80, 0.92],
                                   inside_observed=[dict(flow=r['flow'], lobule_exp=r['lobule_exp'], b=r['b']) for r in out['allometry'] if r['inside_observed_blood_volume']],
                                   combinations=len(out['allometry']))
    if graphs:
        out['graft_on_explicit_graph'] = [graft_on_graph(R, b, seed) for b in (1.0, 0.75) for R in (4.0, 5.0) for seed in (0, 1)]
        z = [g['zR'] for G in out['graft_on_explicit_graph'] for g in G['grafts']]
        out['graft_on_explicit_graph_range'] = [float(min(z)), float(max(z))]
    return out


if __name__ == '__main__':
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'); os.makedirs(f'{root}/results', exist_ok=True)
    s = summary(); json.dump(s, open(f'{root}/results/model_liver.json', 'w'), indent=1)
    print('exponent that reproduces each human cast')
    for k, v in s['casts'].items(): print(f"   {k:15s} diameter ratio {v['d_ratio']:.2f} (model {v['model_ratio_b_half']:.3f} at b = 0.5, {v['model_ratio_b_one']:.3f} at b = 1) -> b = {v['b']:.2f}")
    print('human liver (1.6e6 lobules): share of the pressure drop, shear, lobar graft, exponent of D* with N')
    for r in s['human_liver']:
        print(f"   {r['setting']:18s} b {r['b']:.2f} | inlet {r['share_inlet']:.3f} sinusoid {r['share_sinusoid']:.3f} outlet {r['share_outlet']:.3f} trunks {r['share_trunks']:.3f} | "
              f"shear sinusoid/trunk {r['sinusoid_over_trunk_shear']:.3g} | zR lobar {r['zR_lobar_graft']:.3f}, sectoral {r['zR_sectoral_graft']:.3f} | exponent {r['local_exponent_D_N']:.2f} (asymptote {r['asymptotic_exponent']:.2f})")
    if 'graft_on_explicit_graph_range' in s: print(f"   explicit 3D graphs, one first-order portal branch: zR {s['graft_on_explicit_graph_range'][0]:.2f} to {s['graft_on_explicit_graph_range'][1]:.2f}")
    print('portal flow load with the artery joining at the sinusoid: zR of a lobar graft')
    for r in s['two_inflows']: print(f"   {r['setting']:18s} zF {r['zF']}: arterial flow per lobule x{r['arterial_load']:.1f} -> zR {r['zR']:.3f}")
    for r in s['two_inflows_any_split']: print(f"   any split         zF {r['zF']}: arterial flow per lobule x{r['arterial_load']:.1f} -> zR {r['zR_interval'][0]:.3f} to 1")
    a = s['allometry_ranges']
    print(f"allometric closure over b in [0.5, 1] and the four scaling inputs: tau0 {a['tau0'][0]:+.3f} to {a['tau0'][1]:+.3f}; sinusoidal shear {a['shear_sinusoid'][0]:+.3f} to {a['shear_sinusoid'][1]:+.3f}; "
          f"trunk shear {a['shear_trunk'][0]:+.3f} to {a['shear_trunk'][1]:+.3f}; vascular mass {a['vascular_mass'][0]:.2f} to {a['vascular_mass'][1]:.2f}, {len(a['inside_observed'])} of {a['combinations']} inside 0.80-0.92:")
    for r in a['inside_observed']: print(f"      {r['flow']}, lobule ~ M^{r['lobule_exp']:.2f}, b = {r['b']:.3g}")
