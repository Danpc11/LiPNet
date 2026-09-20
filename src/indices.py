"""Single source of the two indices, the correlations and the pooled logistic fit.

    python src/indices.py            -> results/series_with_indices.tsv, results/logit_zP.json, results/sensitivity.tsv

Definitions (the calculator and plots.py use these same functions)
    zF = graft PVF per 100 g / donor PVF per 100 g   (donor reference from the same series when measured, else 90)
    zP = (PVP - CVP) / 5 mmHg                          (CVP = 5 when not reported; provenance in the *_basis columns)
The stored table holds only raw values and their provenance; the indices are always recomputed here.
Sensitivity: correlations and the fit are repeated excluding groups with any imputed pressure value.
Bootstrap: 2000 binomial resamples of the group event counts, seed 0. Curve grid 0.5-8; the range of the
groups that entered the fit is stored so that the calculator can flag extrapolation.
"""
import json, os
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr

ROOT = os.path.join(os.path.dirname(__file__), '..')
DATA = f'{ROOT}/data/series.tsv'; OUT = f'{ROOT}/results'
NORMAL_GRADIENT = 5.0; DEFAULT_CVP = 5.0; DEFAULT_DONOR_REF = 90.0

def zF(pvf_per_100g, donor_ref=DEFAULT_DONOR_REF):
    return pvf_per_100g / donor_ref

def zP(pvp, cvp=DEFAULT_CVP):
    return (pvp - cvp) / NORMAL_GRADIENT

def compute(d):
    """add zF, zP and an 'imputed' flag to a series table (raw columns only)."""
    d = d.copy()
    d['zF'] = zF(d.PVF_per_100g, d.donor_PVF_per_100g_ref.fillna(DEFAULT_DONOR_REF))
    d['zP'] = zP(d.PVP, d.CVP.fillna(DEFAULT_CVP))
    d['any_imputed'] = d[['PVP_basis', 'CVP_basis']].isin(['imputed']).any(axis=1)
    return d

def fit_logistic(z, events, n):
    def nll(b):
        p = np.clip(1 / (1 + np.exp(-(b[0] + b[1] * z))), 1e-9, 1 - 1e-9)
        return -np.sum(events * np.log(p) + (n - events) * np.log(1 - p))
    return minimize(nll, [-4.0, 1.5]).x

def pooled_fit(s, boot=2000, seed=0):
    z, ev, n = s.zP.values, s.events.values.astype(float), s.n.values.astype(float)
    b = fit_logistic(z, ev, n); rng = np.random.default_rng(seed)
    B = np.array([fit_logistic(z, rng.binomial(n.astype(int), ev / n), n) for _ in range(boot)])
    lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
    zz = np.round(np.arange(0.5, 8.01, 0.1), 2)
    P = 1 / (1 + np.exp(-(b[0] + b[1] * zz))); PB = 1 / (1 + np.exp(-(B[:, [0]] + B[:, [1]] * zz))); ci = np.percentile(PB, [2.5, 97.5], axis=0)
    return dict(b0=float(b[0]), b1=float(b[1]), b0_ci=[float(lo[0]), float(hi[0])], b1_ci=[float(lo[1]), float(hi[1])],
                groups=int(len(s)), events=int(ev.sum()), patients=int(n.sum()), zP_min=float(z.min()), zP_max=float(z.max()),
                normal_gradient_mmHg=NORMAL_GRADIENT, default_cvp_mmHg=DEFAULT_CVP,
                z=zz.tolist(), p=P.round(5).tolist(), lo=ci[0].round(5).tolist(), hi=ci[1].round(5).tolist(),
                groups_used=[dict(study=r.study, group=r.group, n=int(r.n), events=int(r.events), zP=round(float(r.zP), 3),
                                  CVP_basis=r.CVP_basis, PVP_basis=r.PVP_basis) for _, r in s.iterrows()])

def main(boot=2000, seed=0):
    os.makedirs(OUT, exist_ok=True)
    d = compute(pd.read_csv(DATA, sep='\t'))
    d.to_csv(f'{OUT}/series_with_indices.tsv', sep='\t', index=False, float_format='%.3f')
    rows = []
    for label, sub in (('all groups', d), ('no imputed pressure values', d[~d.any_imputed])):
        for col in ('zP', 'zF'):
            s = sub.dropna(subset=[col])
            if len(s) >= 3:
                rho, p = spearmanr(s[col], s.pct); rows.append(dict(analysis=label, index=col, groups=len(s), spearman_rho=rho, p=p))
                print(f'{label:28s} {col}: {len(s):2d} groups, Spearman rho = {rho:.2f}, p = {p:.3f}')
    sf = d[(d.outcome_type == 'SFSS_or_dysfunction') & d.zP.notna()]
    res = pooled_fit(sf, boot, seed); json.dump(res, open(f'{OUT}/logit_zP.json', 'w'), indent=1)
    print(f"logit(SFSS) = {res['b0']:.2f} + {res['b1']:.2f} zP (slope 95% CI {res['b1_ci'][0]:.2f}-{res['b1_ci'][1]:.2f}); "
          f"{res['groups']} groups, {res['events']} events / {res['patients']}; zP range {res['zP_min']:.2f}-{res['zP_max']:.2f}")
    for r in (0.05, 0.10, 0.20): print(f'  zP at {int(r*100)}% risk: {(np.log(r/(1-r)) - res["b0"]) / res["b1"]:.2f}')
    sf2 = sf[~sf.any_imputed]
    if len(sf2) >= 2 and sf2.events.sum() > 0:
        r2 = pooled_fit(sf2, boot, seed); rows.append(dict(analysis='fit without imputed pressure values', index='zP', groups=len(sf2), spearman_rho=np.nan, p=np.nan, b0=r2['b0'], b1=r2['b1']))
    else:
        print(f'Sensitivity fit without imputed pressure values: not estimable ({len(sf2)} SFSS group(s) with fully reported pressures, {int(sf2.events.sum())} events)')
        rows.append(dict(analysis='fit without imputed pressure values', index='zP', groups=len(sf2), spearman_rho=np.nan, p=np.nan, b0=np.nan, b1=np.nan))
    pd.DataFrame(rows).to_csv(f'{OUT}/sensitivity.tsv', sep='\t', index=False, float_format='%.4f')

if __name__ == '__main__':
    main()
