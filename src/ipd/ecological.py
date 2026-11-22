"""What the group means hide: a sensitivity analysis for the spread of the pressure load inside a group.

The primary analysis in src/bayes.py regresses the event count of a group on that group's MEAN pressure load. It is an
association between groups. Every series that reports a dispersion says its recipients did not all sit at the mean:

    Ogura 2010, final PVP <15      gradient 6.2 mmHg   +/- 0.3  (n = 86)
    Ogura 2010, final PVP >=15     gradient 9.9 mmHg   +/- 0.4  (n = 43)
    Uemura 2016, GRWR <0.7         gradient 5.1 mmHg   SD 4.1
    Uemura 2016, GRWR 0.7-0.8      gradient 4.75 mmHg  SD 3.8
    Uemura 2016, GRWR >=0.8        gradient 4.73 mmHg  SD 3.1
    Ishizaki 2012                  gradient 12.4 mmHg  SD 4.4
    Botha 2010                     gradient 6.2 mmHg   SD 3.9  (computed from Table 5)
    Chan 2011                      gradient 6 mmHg     SD 2.7  (from the range, n = 46)
    Wu 2011                        gradient 10.6 mmHg  SD 4.5

Two of those numbers are not standard deviations. Ogura's Results section reports the same cohort's MELD as
23.2 +/- 1.1 (n = 86) and 22.9 +/- 1.7 (n = 43), while its own Table 1 gives the MELD of that cohort as 22.9 +/- 10.6.
Multiplying by sqrt(n) recovers 10.2 and 11.1: the Results-section spreads are standard errors of the mean, whatever
the methods section says. Applying the same conversion to the gradient gives within-group standard deviations of 2.8
and 2.6 mmHg, not 0.3 and 0.4. `sem_to_sd` does this and `WITHIN_GROUP_SD` records it.

Once corrected, the nine dispersions run from 2.6 to 4.5 mmHg, which is comparable to the whole between-group range of
the meta-analysis.

What this module does. It replaces the group-mean likelihood with one that averages an ASSUMED individual curve over
that spread. Patient i of group g is given an unobserved load x_i ~ Normal(mu_g, sigma_g^2), and the event probability
of the group is the average over that distribution,

    pbar_g = Integral logistic(alpha_s + beta_s (x - xbar_s)) Normal(x; mu_g, sigma_g^2) dx,

evaluated by Gauss-Hermite quadrature. events_g ~ Binomial(n_g, pbar_g) as before. With sigma = 0 the integral
collapses to the original model, which is how the sampler is checked against the published fit.

What the slope then is, and what it is not. The data are the same group counts as in src/bayes.py: no outcome is
observed with the load of the recipient who had it, so the slope of an individual relation is not identified by them.
beta is the slope the individual curve would need to have for the group counts to arise IF
    (1) one logistic curve holds for every recipient of a stratum,
    (2) the load is normally distributed inside each group with the spread given, and
    (3) inside a group the load is unrelated to anything else that moves the outcome.
None of the three can be checked with group counts, and the second fails by construction where a group is defined by
a cut-off of the pressure itself (Osman 2017, Ogura 2010), whose spread is bounded by the cut-off. The module therefore
answers one narrow question: how much does the between-group slope depend on having entered the mean of a group where
its recipients are spread? It measures the effect of aggregating a nonlinear curve. It does not remove ecological
confounding, and its slope is not an association observed in patients.

The direction is not fixed either. For small spreads the averaged curve is flatter than the individual one, so the
slope that reproduces the same counts is larger; for spreads comparable to the distance between the groups the counts
stop distinguishing the groups and the estimate returns towards the prior. `sweep` shows both.

Imputed values. Only some groups report a dispersion. `fit` gives each group its own where there is one and a single
value where there is none (the pooled median by default); `imputation` repeats that fit over what is given to the
series that report none, each varied on its own, from no spread to the largest any series reports. An assumed central
venous pressure needs no such analysis: it shifts every group of a study by the same amount and the load is centred
within the study (tests/test_ipd.py).

The sampler is a plain adaptive random-walk Metropolis in NumPy, so this module has no dependency beyond NumPy, SciPy,
pandas and indices.py. At sigma = 0 it reproduces the NUTS fit of src/bayes.py to within Monte Carlo error; that
agreement is the test in tests/test_ipd.py.
"""
import itertools
import numpy as np
import pandas as pd
import indices

GH_NODES = 21                      # Gauss-Hermite nodes; 21 is exact enough for a logistic
MMHG_PER_ZP = indices.NORMAL_GRADIENT   # z_P = (PVP - CVP) / 5 mmHg


def sem_to_sd(sem, n):
    """A standard error reported as if it were a standard deviation."""
    return float(sem) * np.sqrt(float(n))


# Within-group SD of the portocaval gradient in mmHg, by study, after the corrections
# above. Ogura: standard errors converted with sqrt(n). Botha: computed from the published
# patient-level table (patient_level.group_row). Chan: estimated from the published range
# (audit.sd_from_range). Others: as reported.
WITHIN_GROUP_SD = {'Ogura 2010': (2.78, 2.62), 'Uemura 2016': (4.1, 3.8, 3.1),
                   'Ishizaki 2012': (4.4,), 'Botha 2010': (3.91,), 'Chan 2011': (2.71,), 'Wu 2011': (4.5,)}


def corrected_sigma_mmHg():
    """Every within-group SD now available, in mmHg, as a flat list."""
    out = []
    for v in WITHIN_GROUP_SD.values():
        out.extend(v)
    return sorted(out)


def reported_sigma(d):
    """Within-group SD of the pressure load of each group, in z_P units: the value of the table for that group
    (`gradient_sd`) and, for a series whose dispersion is computed here and not printed in its paper (Botha, Chan),
    the value of WITHIN_GROUP_SD. NaN where there is neither."""
    s = []
    for _, row in d.iterrows():
        vals = WITHIN_GROUP_SD.get(str(row.get('study')))
        if not pd.isna(row.get('gradient_sd')):
            s.append(float(row['gradient_sd']) / MMHG_PER_ZP)
        elif vals is not None and len(vals) == 1:
            s.append(float(vals[0]) / MMHG_PER_ZP)
        else:
            s.append(np.nan)
    return np.array(s)


def pooled_sigma(kind='median'):
    """A single within-group SD to stand in where none is reported, pooled over every
    series that reports one, after the standard-error correction above."""
    s = np.array(corrected_sigma_mmHg()) / MMHG_PER_ZP
    return float(np.median(s) if kind == 'median' else np.mean(s))


def _loglik(theta, S, mu_g, sigma_g, xbar_s, n, ev, n_s, gh_x, gh_w,
            alpha_scale, mu_scale, tau_scale):
    alpha = theta[:n_s]
    raw = theta[n_s:2 * n_s]
    mu_b = theta[2 * n_s]
    log_tau = theta[2 * n_s + 1]
    tau = np.exp(log_tau)
    beta = mu_b + tau * raw

    a = alpha[S]
    b = beta[S]
    # x = mu_g + sqrt(2) sigma_g t,  weights gh_w / sqrt(pi)
    x = mu_g[:, None] + np.sqrt(2.0) * sigma_g[:, None] * gh_x[None, :]
    eta = a[:, None] + b[:, None] * (x - xbar_s[S][:, None])
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -500, 500)))
    pbar = np.clip((p * gh_w[None, :]).sum(axis=1) / np.sqrt(np.pi), 1e-12, 1 - 1e-12)

    ll = float(np.sum(ev * np.log(pbar) + (n - ev) * np.log1p(-pbar)))
    lp = (-0.5 * np.sum((alpha / alpha_scale) ** 2)
          - 0.5 * np.sum(raw ** 2)
          - 0.5 * (mu_b / mu_scale) ** 2
          - 0.5 * (tau / tau_scale) ** 2 + log_tau)        # HalfNormal + Jacobian
    return ll + lp


def usable(d, by_outcome=True):
    """The rows a fit uses: strata (a study and, by default, an outcome definition) with at least two groups at
    different pressure loads."""
    d = d.copy()
    d['stratum'] = d.study.astype(str) + ((' / ' + d.outcome_type.astype(str)) if by_outcome else '')
    return d.groupby('stratum').filter(lambda x: len(x) >= 2 and x.zP.nunique() > 1)


def fit(d, sigma=None, fill=None, chains=4, warmup=4000, samples=8000, seed=0,
        alpha_scale=2.5, mu_scale=1.5, tau_scale=0.5, by_outcome=True):
    """Fit the integrated model.

    sigma : None uses each group's own within-group SD where the series gives one and `fill` where it does not; a float
            applies that SD (in z_P units) to every group; 0 reproduces the group-mean model of src/bayes.py.
    fill  : with sigma None, the SD (in z_P units) given to the groups with none: one value, or one per study in a
            dict (studies not named take the pooled median). None is the pooled median for all of them.
    """
    d = usable(d, by_outcome)
    strata = sorted(d.stratum.unique())
    S = np.array([strata.index(s) for s in d.stratum])
    n_s = len(strata)
    mu_g = d.zP.values.astype(float)
    xbar_s = np.array([d.loc[d.stratum == s, 'zP'].mean() for s in strata])
    n = d.n.values.astype(float)
    ev = d.events.values.astype(float)

    if sigma is None:
        sg = reported_sigma(d)
        missing = np.isnan(sg)
        pooled = pooled_sigma()
        if isinstance(fill, dict):
            given = np.array([float(fill.get(str(st), pooled)) for st in d.study])
            sigma_label = 'reported; where missing, ' + ', '.join('%s %.2f' % (k, v) for k, v in sorted(fill.items()))
        else:
            given = np.full(len(d), pooled if fill is None else float(fill))
            sigma_label = ('reported, pooled median %.2f where missing' % pooled) if fill is None else 'reported, %.2f where missing' % float(fill)
        sg = np.where(missing, given, sg)
    else:
        sg = np.full(len(d), float(sigma))
        missing = np.zeros(len(d), bool)
        sigma_label = 'fixed at %.2f' % float(sigma)

    gh_x, gh_w = np.polynomial.hermite.hermgauss(GH_NODES)
    args = (S, mu_g, sg, xbar_s, n, ev, n_s, gh_x, gh_w, alpha_scale, mu_scale, tau_scale)
    dim = 2 * n_s + 2

    draws = []
    acc_rates = []
    for c in range(chains):
        rng = np.random.default_rng(seed * 100 + c)
        th = np.concatenate([np.zeros(n_s), np.zeros(n_s), [0.0], [np.log(0.3)]])
        th += 0.1 * rng.standard_normal(dim)
        lp = _loglik(th, *args)
        step = np.full(dim, 0.25)
        chain = np.empty((samples, dim))
        acc = 0
        for it in range(warmup + samples):
            prop = th + step * rng.standard_normal(dim)
            lp_p = _loglik(prop, *args)
            if np.log(rng.random()) < lp_p - lp:
                th, lp = prop, lp_p
                acc += 1
            if it < warmup:
                if (it + 1) % 100 == 0:
                    rate = acc / (it + 1)
                    step *= np.exp((rate - 0.234) * 1.5)
            else:
                chain[it - warmup] = th
        draws.append(chain)
        acc_rates.append(acc / (warmup + samples))
    draws = np.array(draws)                                   # (chains, samples, dim)

    mu_draws = draws[:, :, 2 * n_s].reshape(-1)
    tau_draws = np.exp(draws[:, :, 2 * n_s + 1]).reshape(-1)
    q = lambda x: [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]
    return dict(sigma=sigma_label,
                sigma_by_group=[float(x) for x in sg],
                strata=strata, groups=int(len(d)), groups_imputed=int(missing.sum()), patients_imputed=int(n[missing].sum()),
                patients=int(n.sum()), events=int(ev.sum()),
                mu_beta=float(mu_draws.mean()), mu_beta_ci=q(mu_draws),
                tau_beta=float(tau_draws.mean()), tau_beta_ci=q(tau_draws),
                prob_mu_positive=float((mu_draws > 0).mean()),
                OR_per_zP=float(np.exp(mu_draws.mean())),
                OR_per_mmHg=float(np.exp(mu_draws.mean() / MMHG_PER_ZP)),
                rhat=float(_rhat(draws[:, :, 2 * n_s])),
                ess=float(_ess(draws[:, :, 2 * n_s])),
                acceptance=float(np.mean(acc_rates)))


def _rhat(x):
    """Split R-hat on a (chains, samples) array."""
    c, s = x.shape
    h = s // 2
    y = np.concatenate([x[:, :h], x[:, h:2 * h]], axis=0)
    m, nn = y.shape
    means = y.mean(axis=1)
    W = y.var(axis=1, ddof=1).mean()
    B = nn * means.var(ddof=1)
    var = ((nn - 1) * W + B) / nn
    return np.sqrt(var / W) if W > 0 else np.nan


def _ess(x):
    """Crude effective sample size from the lag-1 autocorrelation of the pooled chains."""
    f = x.reshape(-1)
    f = f - f.mean()
    r1 = np.sum(f[1:] * f[:-1]) / np.sum(f ** 2)
    return len(f) * (1 - r1) / (1 + r1)


def sweep(d, sigmas=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), **kw):
    """The slope of the integrated model as a function of one within-group spread assumed for every group.

    sigma is in z_P units, so 0.8 corresponds to a within-group SD of 4 mmHg, which is what Uemura and Ishizaki report.
    `ratio_to_sigma0` is the slope over the group-mean slope (sigma = 0). It is above one for spreads below the distance
    between the groups and falls back, towards the prior, beyond it: the sweep states the direction for each spread and
    does not assume one.
    """
    out = []
    for s in sigmas:
        r = fit(d, sigma=s, **kw)
        out.append(dict(sigma_zP=float(s), sigma_mmHg=float(s) * MMHG_PER_ZP,
                        mu_beta=r['mu_beta'], mu_beta_ci=r['mu_beta_ci'],
                        OR_per_zP=r['OR_per_zP'], OR_per_mmHg=r['OR_per_mmHg'],
                        prob_positive=r['prob_mu_positive'], rhat=r['rhat']))
    base = out[0]['mu_beta']
    for o in out:
        o['ratio_to_sigma0'] = o['mu_beta'] / base if base != 0 else float('nan')
    return out


def imputation(d, values_mmHg=None, by_outcome=True, **kw):
    """The fit with the reported dispersions, repeated over what is given to the series that report none.

    Each series without a dispersion takes, on its own, every value of `values_mmHg` (by default no spread, the smallest
    and the largest dispersion any series reports, and their pooled median), so two such series give sixteen fits. The
    groups that report a dispersion keep it throughout. Returns the fits and the range of the slope over them."""
    f = usable(d, by_outcome); missing = np.isnan(reported_sigma(f)); series = sorted(set(f.study[missing].astype(str)))
    sds = corrected_sigma_mmHg()
    values = [float(v) for v in (values_mmHg if values_mmHg is not None else (0.0, sds[0], pooled_sigma() * MMHG_PER_ZP, sds[-1]))]
    fits = []
    for combo in itertools.product(values, repeat=len(series)):
        r = fit(d, sigma=None, fill={st: v / MMHG_PER_ZP for st, v in zip(series, combo)}, by_outcome=by_outcome, **kw)
        fits.append(dict(sd_mmHg=dict(zip(series, combo)), mu_beta=r['mu_beta'], mu_beta_ci=r['mu_beta_ci'], prob_positive=r['prob_mu_positive'], rhat=r['rhat']))
    mu = [x['mu_beta'] for x in fits]
    return dict(series_without_dispersion=series, values_mmHg=values, groups=int(len(f)), groups_imputed=int(missing.sum()),
                patients=int(f.n.sum()), patients_imputed=int(f.n[missing].sum()), fits=fits,
                mu_beta_range=[float(min(mu)), float(max(mu))], lowest_lower_bound=float(min(x['mu_beta_ci'][0] for x in fits)),
                highest_upper_bound=float(max(x['mu_beta_ci'][1] for x in fits)), prob_positive_min=float(min(x['prob_positive'] for x in fits)),
                rhat_max=float(max(x['rhat'] for x in fits)))
