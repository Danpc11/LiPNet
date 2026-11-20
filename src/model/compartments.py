"""Distension site by site: which vessels distend, and at which pressure.

distension.py applies each published relation to the whole liver: one factor d(G) multiplies the resistance of the
rigid network. That is the extreme case in which every vessel distends as the whole organ did in the cat. This module
replaces the global factor by the vessels of the network itself.

Sites. Lautt, Greenway and Legare (Am J Physiol 1991;260:G325) validated the cube law separately for the
presinusoidal and the postsinusoidal resistance sites, each with its own distending pressure (the mean of the
pressures upstream and downstream of that site), and showed that the indices of the two sites add in series. The same
form is applied here to every generation of the network: generation i keeps its own index, R_i = IC_i / Pd_i^3, with
Pd_i the mean pressure in that generation, which the network's own pressure profile gives. For the compliance
relation each generation changes its volume with its own change of pressure, and R_i goes as V_i^-2.

Which sites distend. In basal conditions the resistance of the cat liver lies mainly at postsinusoidal sites in the
hepatic veins, and presinusoidal resistance can reach half of the total under constriction (Lautt et al. 1991). The
network places 42 to 48% of the donor drop in the portal tree, 3 to 20% in the sinusoids and 37 to 48% in the hepatic
venous tree. Four distributions are compared, each with the share s of the resistance that can distend fixed by the
anatomy of the network, not by the six groups:
    every vessel                      s = 1
    pre- and postsinusoidal sites     s = 1 - share of the sinusoids (the two sites of Lautt et al.)
    postsinusoidal site               s = share of the hepatic venous tree (the basal resistance site of the cat)
    presinusoidal site                s = share of the portal tree
The vessels that do not distend keep their rigid resistance: R / R0 = (1 - s) + s d, with d now a weighted mean over
the distending vessels at their own pressures.

Distending pressure. The global relations of distension.py take Pd = CVP + G / 2 for the cube law and the portal
pressure CVP + G for the compliance relation. `pressure_variants` replaces them by Pd = CVP - Pext + lambda G with
lambda from 0.25 to 1 for both relations, and an external pressure Pext of 0 or 2 mmHg where the outflow pressure stays
above 1 mmHg. The site-by-site model needs no lambda: each vessel sits at the pressure the network gives it.

    python src/model/compartments.py     -> printed; the results enter results/model_distension.json
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import indices
from model import liver, human

A_RATIO = liver.ARTERIAL_FRACTION / (1 - liver.ARTERIAL_FRACTION)
SITES = ('every vessel', 'pre- and postsinusoidal sites', 'postsinusoidal site', 'presinusoidal site')
LAMBDAS = (0.25, 0.5, 0.75, 1.0)
EXTERNAL = (0.0, 2.0)


def geometries():
    """Donor pressure-drop shares by generation of the networks that make up the rigid range: the canopy at b = 1,
    0.75, 0.5 and the calibrated human liver at b = 0.6 to 0.93. Returns (name, shares, k, trunk share removed)."""
    out = []
    for b in (1.0, 0.75, 0.5):
        c = liver.canopy(liver.HUMAN_GENERATIONS, b); dP = c['dP'] / c['dP'].sum()
        out.append((f'canopy b={b:g}', dP, int(c['k']), float(dP[0] + dP[-1])))
    for b in (0.6, 2 / 3, 0.75, 0.85, 0.93):
        h = human.liver(b); dP = h['dP_mmHg'] / h['dP_mmHg'].sum()
        out.append((f'human b={b:.2f}', dP, int(h['g']), float(h['trunks_removed'])))
    return out


def mask(site, n, k):
    m = np.zeros(n)
    if site == 'every vessel': m[:] = 1
    elif site == 'pre- and postsinusoidal sites': m[:] = 1; m[k] = 0
    elif site == 'postsinusoidal site': m[k + 1:] = 1
    elif site == 'presinusoidal site': m[:k] = 1
    return m


def share(site, w, k):
    return float((mask(site, len(w), k) * w).sum())


def _profile(drop, cvp):
    """Mean pressure in each generation, from the drops of all generations downstream of it."""
    down = np.concatenate([np.cumsum(drop[::-1])[::-1][1:], [0.0]])
    return cvp + down + drop / 2


def solve(zF, geom, phi, cvp, law, p=None, site='every vessel', g0=indices.NORMAL_GRADIENT, external=0.0, iters=600):
    """zR of a graft whose vessels in `site` follow `law` at their own mean pressure ('cube law' or 'compliance', with
    doubling pressure p). The rigid limit reproduces liver.two_inflows. Returns zR, zP and the share s that distends."""
    name, w, k, rem = geom; n = len(w); m = mask(site, n, k)
    load = np.where(np.arange(n) < k, zF, (zF + phi * A_RATIO) / (1 + A_RATIO))
    P0 = _profile(g0 * w, cvp) - external
    d = np.ones(n)
    for _ in range(iters):
        drop = g0 * (1 - rem) * w * load * d; P = _profile(drop, cvp) - external
        if law == 'cube law': new = (P0 / np.maximum(P, 1e-6)) ** 3
        elif law == 'compliance': new = np.maximum(1.0 + (P - P0) / p, 1e-6) ** -2.0
        else: raise ValueError(law)
        new = m * new + (1 - m)
        step = np.exp(0.5 * np.log(new) + 0.5 * np.log(d))
        if np.max(np.abs(step - d)) < 1e-10: d = step; break
        d = step
    zP = float((1 - rem) * (w * load * d).sum()); return dict(zR=zP / zF, zP=zP, s=float((m * w).sum()))


def rigid(zF, geom, phi):
    name, w, k, rem = geom; n = len(w); load = np.where(np.arange(n) < k, zF, (zF + phi * A_RATIO) / (1 + A_RATIO))
    return float((1 - rem) * (w * load).sum() / zF)


def global_variant(zF, zR_rigid, cvp, law, p, lam, external=0.0, g0=indices.NORMAL_GRADIENT):
    """distension.distensible with the distending pressure CVP - Pext + lambda G."""
    out = cvp - external
    if law == 'cube law': d = lambda G: ((out + lam * g0) / (out + lam * G)) ** 3
    else: d = lambda G: np.maximum(1.0 + lam * (G - g0) / p, 1e-9) ** -2.0
    lo, hi = 0.0, g0 * zF * zR_rigid * 2 + 1; f = lambda G: g0 * zF * zR_rigid * d(G) - G
    for _ in range(80):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi) / g0 / zF


def site_ranges(groups, cvps_default, doubling):
    """For each group and each distribution of distensible sites: the range of zR over the networks, the arterial
    load (none, the donor's, scaling with the portal one), the central venous pressure and the doubling pressure."""
    G = geometries(); out = []
    for g in groups:
        cvps = (g['cvp'],) if g['cvp'] is not None else cvps_default; row = dict(study=g['study'], group=g['group'], by_site={})
        for site in SITES:
            res = {}
            for law, ps in (('cube law', (None,)), ('compliance', doubling)):
                z = [solve(g['zF'], geom, phi, c, law, p, site)['zR'] for geom in G for phi in (0.0, 1.0, g['zF']) for c in cvps for p in ps]
                res[law] = [round(min(z), 3), round(max(z), 3)]
            row['by_site'][site] = dict(share_that_distends=[round(min(share(site, x[1], x[2]) for x in G), 3), round(max(share(site, x[1], x[2]) for x in G), 3)], **res)
        out.append(row)
    return out


def pressure_variants(groups, cvps_default, doubling, rigid_range):
    """The global relations with the distending pressure CVP - Pext + lambda G: the range of zR for each group over
    lambda, Pext, the rigid range, the central venous pressure and the doubling pressure."""
    out = []
    for g in groups:
        lo, hi = rigid_range(g['zF']); cvps = (g['cvp'],) if g['cvp'] is not None else cvps_default; row = dict(study=g['study'], group=g['group'], by_lambda={})
        for lam in LAMBDAS:
            res = {}
            for law, ps in (('cube law', (None,)), ('compliance', doubling)):
                z = [global_variant(g['zF'], zr, c, law, p, lam, e) for zr in (lo, hi) for c in cvps for e in EXTERNAL if c - e >= 1.0 for p in ps]
                res[law] = [round(min(z), 3), round(max(z), 3)]
            row['by_lambda'][f'{lam:g}'] = res
        out.append(row)
    return out


if __name__ == '__main__':
    import pandas as pd
    from model import distension as ds
    d = indices.compute(pd.read_csv(os.path.join(ds.ROOT, 'data/series.tsv'), sep='\t')); gs = ds.groups(d)
    geom = geometries()
    for x in geom: print(x[0], 'pre', round(share('presinusoidal site', x[1], x[2]), 3), 'post', round(share('postsinusoidal site', x[1], x[2]), 3))
    for g in gs:
        lo, hi = ds.rigid_range(g['zF']); chk = [rigid(g['zF'], x, phi) for x in geom for phi in (0.0, 1.0, g['zF'])]
        print(g['study'], round(lo, 3), round(hi, 3), 'recomputed', round(min(chk), 3), round(max(chk), 3))
    for r in site_ranges(gs, ds.CVPS, ds.DOUBLING_MMHG): print(r)
