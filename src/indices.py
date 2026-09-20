"""Compute the two indices and the pooled logistic fit from data/series.tsv.

    python src/indices.py            -> results/series_with_indices.tsv, results/logit_zP.json

Definitions
    zF = graft PVF per 100 g / donor PVF per 100 g   (donor reference from the same series when measured, else 90)
    zP = (PVP - CVP) / 5 mmHg                          (CVP = 5 when not reported; flagged in CVP_measured)
The logistic fit uses the groups reporting SFSS or early graft dysfunction counts together with zP.
Bootstrap: 2000 resamples of the group event counts (binomial), seed 0.
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr

DATA = os.path.join(os.path.dirname(__file__), '..', 'data', 'series.tsv')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results')

def zF(pvf_per_100g, donor_ref=90.0):
    return pvf_per_100g / donor_ref

def zP(pvp, cvp=5.0, normal_gradient=5.0):
    return (pvp - cvp) / normal_gradient

def fit_logistic(z, events, n):
    def nll(b):
        p = np.clip(1 / (1 + np.exp(-(b[0] + b[1] * z))), 1e-9, 1 - 1e-9)
        return -np.sum(events * np.log(p) + (n - events) * np.log(1 - p))
    return minimize(nll, [-4.0, 1.5]).x

def main(boot=2000, seed=0):
    os.makedirs(OUT, exist_ok=True)
    d = pd.read_csv(DATA, sep='\t')
    # recompute the indices from the raw columns (they are stored, but this is the reference implementation)
    d['zF_calc'] = zF(d.PVF_per_100g, d.donor_PVF_per_100g_ref.fillna(90))
    d['zP_calc'] = zP(d.PVP, d.CVP.fillna(5))
    # groups whose zP came from a reported gradient keep the stored value
    d['zP_calc'] = d.zP_calc.where(d.zP_calc.notna(), d.zP)
    d.to_csv(f'{OUT}/series_with_indices.tsv', sep='\t', index=False, float_format='%.3f')

    for col in ('zP', 'zF'):
        s = d.dropna(subset=[col])
        rho, p = spearmanr(s[col], s.pct)
        print(f'{col}: {len(s)} groups, Spearman rho = {rho:.2f}, p = {p:.3f}')

    s = d[(d.outcome_type == 'SFSS_or_dysfunction') & d.zP.notna()]
    z, ev, n = s.zP.values, s.events.values.astype(float), s.n.values.astype(float)
    b = fit_logistic(z, ev, n)
    rng = np.random.default_rng(seed)
    B = np.array([fit_logistic(z, rng.binomial(n.astype(int), ev / n), n) for _ in range(boot)])
    lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
    zz = np.linspace(0.8, 4.0, 33)
    P = 1 / (1 + np.exp(-(b[0] + b[1] * zz)))
    PB = 1 / (1 + np.exp(-(B[:, [0]] + B[:, [1]] * zz)))
    ci = np.percentile(PB, [2.5, 97.5], axis=0)
    res = dict(b0=float(b[0]), b1=float(b[1]), b0_ci=[float(lo[0]), float(hi[0])], b1_ci=[float(lo[1]), float(hi[1])],
               groups=int(len(s)), events=int(ev.sum()), patients=int(n.sum()),
               z=zz.round(4).tolist(), p=P.round(5).tolist(), lo=ci[0].round(5).tolist(), hi=ci[1].round(5).tolist(),
               groups_used=[dict(study=r.study, group=r.group, n=int(r.n), events=int(r.events), zP=float(r.zP)) for _, r in s.iterrows()])
    json.dump(res, open(f'{OUT}/logit_zP.json', 'w'), indent=1)
    for r in (0.05, 0.10, 0.20):
        print(f'zP at {int(r*100)}% risk: {(np.log(r/(1-r)) - b[0]) / b[1]:.2f}')
    print(f'logit(SFSS) = {b[0]:.2f} + {b[1]:.2f} zP   (slope 95% CI {lo[1]:.2f}-{hi[1]:.2f}); {len(s)} groups, {int(ev.sum())} events / {int(n.sum())}')

if __name__ == '__main__':
    main()
