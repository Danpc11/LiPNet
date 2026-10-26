"""What the lobule does to the clinical indices, once the network is allowed to dilute it.

The FEM says the lobule's own resistance falls when its flow rises: z_R_lobule(4) = 0.56.
The clinical index z_R is measured across the whole graft, and the sinusoidal bed is only
one resistance in a series that also contains the portal tree, the hepatic venous tree and
two anastomoses. So the question this module asks is not whether the effect exists but how
much of it survives.

Two steps.

1. `sinusoidal_share` takes the dissipation-optimal network of src/model/network.py and
   measures the fraction of its dissipation that sits in the sinusoidal edges. That is the
   only part of the series that the lobule model governs.

2. `graft` repeats the series calculation of src/model/whole_graft.py with the lobular
   share made flow dependent, and solves the fixed point, since a lower resistance admits
   more flow and more flow dilates further.

The answer is in `sweep`, and it is mostly a negative one. That is the useful thing to
know before anyone attributes the published z_R values to sinusoidal autoregulation.
"""
import numpy as np


# ------------------------------------------------------------------ network share
def sinusoidal_share(R=4.5, b=0.75, s=1.0, sweeps=4, seed=0):
    """Fraction of the inlet dissipation carried by the sinusoidal edges at the optimum.

    This is a property of the cost-constrained optimum, not a measurement: the optimiser
    assigns the conductance of every edge under a budget, so the sinusoidal share is what
    the model thinks it should be, not what a cast of a liver would show. It is reported
    for the inlet problem (portal tree plus sinusoids); adding the hepatic venous tree can
    only lower the share, so this is an upper bound.
    """
    from model import network as nw
    G = nw.build_liver_graph(R=R, s=s, seed=seed)
    inlet, _ = nw.split_liver(G)
    sol = nw.best_sinusoid_tree(inlet, b=b, sweeps=sweeps, seed=seed)
    keep = sol['keep']
    D, _, f = inlet.solve(sol['w'], keep)
    w = sol['w'][keep]
    diss = f[keep] ** 2 / np.maximum(w, 1e-300)
    kind = inlet.kind[keep]
    return dict(R=float(R), b=float(b), n_lobules=int(G['n_lob']),
                sinusoidal_share=float(diss[kind == 'sin'].sum() / diss.sum()),
                dissipation=float(D))


def share_scaling(radii=(3.0, 4.5, 6.5), b=0.75):
    """How the sinusoidal share moves with the size of the liver. It falls slowly, which
    matters because a real liver has of order a million lobules and these graphs have
    hundreds: extrapolating the trend is the only way to guess the real share, and four
    decades of extrapolation is not a measurement."""
    rows = [sinusoidal_share(R=r, b=b) for r in radii]
    n = np.array([r['n_lobules'] for r in rows], dtype=float)
    sh = np.array([r['sinusoidal_share'] for r in rows])
    slope = float(np.polyfit(np.log(n), np.log(sh), 1)[0])
    return dict(rows=rows, log_slope=slope,
                extrapolated_at_1e6=float(sh[0] * (1e6 / n[0]) ** slope))


# ------------------------------------------------------------------ the graft
def zr_interpolator(load):
    """z_R_lobule(z_F) from coupling.resistance_load, log-log, flat outside the range."""
    zf = np.asarray(load['z_F'], dtype=float)
    zr = np.asarray(load['z_R'], dtype=float)
    order = np.argsort(zf)
    lx, ly = np.log(zf[order]), np.log(zr[order])

    def curve(x):
        return float(np.exp(np.interp(np.log(max(float(x), 1e-6)), lx, ly)))
    return curve


def graft(h, share, curve=None, r_in=0.05, r_out=0.05, r_coll=12.0,
          d_in=1.0, d_out=1.0, tol=1e-10, max_iter=100):
    """The series model of whole_graft.py with the sinusoidal part made flow dependent.

    h      recipient splanchnic inflow relative to the donor's
    share  fraction of the lobular resistance that is sinusoidal, from sinusoidal_share
    curve  z_R_lobule(z_F); None means the rigid Poiseuille bed, i.e. whole_graft.py

    Returns the same indices as whole_graft, plus the lobule's own z_R and the part of it
    that reaches the measurement.
    """
    r_lob = 1.0 - r_in - r_out
    r_sin = share * r_lob
    r_tree = r_lob - r_sin
    R_in, R_out = r_in * d_in ** -4, r_out * d_out ** -4
    par = lambda RR: 1.0 / (1.0 / RR + 1.0 / r_coll)
    f0 = h * par(1.0)

    m, f = 1.0, h
    for _ in range(max_iter):
        R_lob = r_tree + r_sin * m
        R = R_in + R_lob + R_out
        f_abs = h * par(R) / R
        f_new = h * f_abs / f0
        m_new = 1.0 if curve is None else curve(f_new)
        if abs(f_new - f) < tol and abs(m_new - m) < tol:
            f, m = f_new, m_new
            break
        f, m = f_new, 0.5 * (m + m_new)
    R_lob = r_tree + r_sin * m
    zP = f * (R_lob + R_out) / (r_lob + r_out)
    return dict(zF=float(f), zP=float(zP), zP_over_zF=float(zP / f) if f > 0 else np.nan,
                z_R_lobule=float(m), z_R_graft=float(zP / f) if f > 0 else np.nan,
                sinusoidal_share=float(share),
                R_total=float(R_in + R_lob + R_out))


def sweep(load, share, inflow=(1.0, 1.5, 2.0, 3.0, 4.0), **kw):
    """Rigid against adaptive, side by side."""
    curve = zr_interpolator(load)
    out = dict(inflow=[], zF_rigid=[], zF_adaptive=[], zR_lobule=[], zR_graft=[],
               zP_rigid=[], zP_adaptive=[], share=float(share))
    for h in inflow:
        a = graft(h, share, curve, **kw)
        r = graft(h, share, None, **kw)
        out['inflow'].append(float(h))
        out['zF_rigid'].append(r['zF'])
        out['zF_adaptive'].append(a['zF'])
        out['zR_lobule'].append(a['z_R_lobule'])
        out['zR_graft'].append(a['z_R_graft'])
        out['zP_rigid'].append(r['zP'])
        out['zP_adaptive'].append(a['zP'])
    zl = np.array(out['zR_lobule'])
    zg = np.array(out['zR_graft'])
    with np.errstate(divide='ignore', invalid='ignore'):
        transmitted = np.where(zl < 1, (1 - zg) / (1 - zl), np.nan)
    out['fraction_transmitted'] = [None if not np.isfinite(t) else float(t) for t in transmitted]
    return out


def share_vs_cost_exponent(radii=(3.0, 4.5, 6.5, 9.0), exponents=(0.5, 0.75, 1.0)):
    """The sinusoidal share against liver size, for each cost exponent b.

    This is the test that decides whether the share's drift with N is a problem or a
    constraint. It is slow: the optimisation at R = 9 takes about half a minute per b.
    """
    out = []
    for b in exponents:
        rows = [sinusoidal_share(R=r, b=b) for r in radii]
        n = np.array([r['n_lobules'] for r in rows], dtype=float)
        sh = np.array([r['sinusoidal_share'] for r in rows])
        slope = float(np.polyfit(np.log(n), np.log(sh), 1)[0])
        out.append(dict(b=float(b), n_lobules=n.tolist(), share=sh.tolist(),
                        log_slope=slope,
                        extrapolated_1e6=float(sh[0] * (1e6 / n[0]) ** slope)))
    return out


def predict_series(load, share, root=None, **kw):
    """The prediction against the five published groups that report both indices.

    The flow load barely varies across them (z_F from 3.3 to 4.1), so the lobule model
    predicts a z_R that is the same to three decimals for all five. Whatever separates
    0.31 from 0.74 in the data is not flow through the lobules.
    """
    import os
    import pandas as pd
    import sys
    here = os.path.dirname(os.path.abspath(__file__))
    root = root or os.path.join(here, '..', '..')
    sys.path.insert(0, os.path.join(root, 'src'))
    import indices
    d = indices.compute(pd.read_csv(os.path.join(root, 'data', 'series.tsv'), sep='\t'))
    d = d.dropna(subset=['zF', 'zP'])
    curve = zr_interpolator(load)
    rows = []
    for _, r in d.iterrows():
        g = graft(float(r.zF), share, curve, **kw)
        rows.append(dict(study=str(r.study), outflow=str(r.outflow), zF=float(r.zF),
                         zR_observed=float(r.zR), zR_lobule=g['z_R_lobule'],
                         zR_predicted=g['z_R_graft'],
                         unexplained=float(r.zR) / g['z_R_graft']))
    obs = np.array([x['zR_observed'] for x in rows])
    pred = np.array([x['zR_predicted'] for x in rows])
    return dict(groups=rows,
                observed_spread=float(obs.max() - obs.min()),
                predicted_spread=float(pred.max() - pred.min()),
                spread_ratio=float((obs.max() - obs.min()) / (pred.max() - pred.min())))
