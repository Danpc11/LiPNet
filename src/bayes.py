"""Primary analysis: hierarchical binomial model of the outcome on the within-study centred pressure index.

    E_gs ~ Binomial(n_gs, p_gs)
    logit(p_gs) = alpha_s + beta_s (zP_gs - mean zP_s)
    beta_s ~ Normal(mu_beta, tau_beta)

Centring zP within each study separates the study's own level (alpha_s) from the association (beta_s), so the
two are not traded off against each other. mu_beta is the quantity of interest; tau_beta is the variation of the
association between studies. Fitted with NUTS (numpyro), four chains, with R-hat, effective sample size,
divergences, posterior predictive checks and sensitivity to the prior on tau_beta.

Outcome order, fixed in advance:
  1. primary      SFSS or early graft dysfunction
  2. secondary    mortality or graft loss
  3. exploratory  both combined (one stratum per study and outcome definition)
"""
import os
os.environ.setdefault('XLA_FLAGS', '--xla_force_host_platform_device_count=4')
import numpy as np, pandas as pd
import jax, jax.numpy as jnp, numpyro, numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS, Predictive
numpyro.set_host_device_count(4)

PRIORS = (0.25, 0.5, 1.0)           # half-normal scales for tau_beta, reported side by side


def _model(S, z_c, n, events=None, n_studies=1, tau_scale=0.5, alpha_scale=2.5, mu_scale=1.5):
    alpha = numpyro.sample('alpha', dist.Normal(0., alpha_scale).expand([n_studies]))
    mu = numpyro.sample('mu_beta', dist.Normal(0., mu_scale))
    tau = numpyro.sample('tau_beta', dist.HalfNormal(tau_scale))
    with numpyro.plate('studies', n_studies):
        beta_raw = numpyro.sample('beta_raw', dist.Normal(0., 1.))
    beta = numpyro.deterministic('beta', mu + tau * beta_raw)
    logit = alpha[S] + beta[S] * z_c
    numpyro.sample('obs', dist.BinomialLogits(logit, total_count=n), obs=events)


def _prep(d, by_outcome=True, by_centre=False):
    """Return the design: one stratum per study (or per centre), optionally split by outcome definition, with zP
    centred within stratum. by_centre=True pools publications from the same centre, which is the right unit when
    their recruitment periods overlap (Kyoto)."""
    d = d.copy()
    unit = d.centre_id.astype(str) if by_centre else d.study.astype(str)
    d['stratum'] = unit + ((' / ' + d.outcome_type.astype(str)) if by_outcome else '')
    d = d.groupby('stratum').filter(lambda x: len(x) >= 2 and x.zP.nunique() > 1)
    strata = sorted(d.stratum.unique())
    S = np.array([strata.index(s) for s in d.stratum])
    z_c = (d.zP - d.groupby('stratum').zP.transform('mean')).values
    return d, strata, S, z_c


def fit(d, tau_scale=0.5, chains=4, warmup=1500, samples=2000, seed=0, by_outcome=True, by_centre=False, quiet=True, keep_draws=0):
    d, strata, S, z_c = _prep(d, by_outcome, by_centre)
    n = d.n.values.astype(int); ev = d.events.values.astype(int)
    kernel = NUTS(_model, target_accept_prob=0.95)
    mcmc = MCMC(kernel, num_warmup=warmup, num_samples=samples, num_chains=chains, progress_bar=not quiet)
    mcmc.run(jax.random.PRNGKey(seed), jnp.array(S), jnp.array(z_c), jnp.array(n), events=jnp.array(ev),
             n_studies=len(strata), tau_scale=tau_scale)
    post = mcmc.get_samples(group_by_chain=True)
    flat = {k: np.array(v).reshape((-1,) + np.array(v).shape[2:]) for k, v in post.items()}
    summ = numpyro.diagnostics.summary(post, prob=0.95)
    div = int(np.sum(mcmc.get_extra_fields(group_by_chain=False).get('diverging', np.zeros(1)))) if 'diverging' in mcmc.get_extra_fields() else 0
    mu = flat['mu_beta']; tau = flat['tau_beta']
    # posterior predictive check: replicated event counts
    pred = Predictive(_model, posterior_samples={k: v for k, v in flat.items() if k in ('alpha', 'beta_raw', 'mu_beta', 'tau_beta')})
    rep = np.array(pred(jax.random.PRNGKey(seed + 1), jnp.array(S), jnp.array(z_c), jnp.array(n),
                        n_studies=len(strata), tau_scale=tau_scale)['obs'])
    ppc = [dict(stratum=d.stratum.iloc[i], group=d.group.iloc[i], observed=int(ev[i]),
                predicted_median=float(np.median(rep[:, i])),
                predicted_ci=[float(np.percentile(rep[:, i], 2.5)), float(np.percentile(rep[:, i], 97.5))],
                inside=bool(np.percentile(rep[:, i], 2.5) <= ev[i] <= np.percentile(rep[:, i], 97.5))) for i in range(len(ev))]
    q = lambda x: [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]
    rng = np.random.default_rng(seed)
    pred_new = mu + tau * rng.standard_normal(len(mu))              # a new study's slope
    return dict(tau_prior=tau_scale, strata=strata, groups=int(len(d)), events=int(ev.sum()), patients=int(n.sum()),
                mu_beta=float(mu.mean()), mu_beta_ci=q(mu), tau_beta=float(tau.mean()), tau_beta_ci=q(tau),
                prob_mu_positive=float((mu > 0).mean()),
                OR_per_zP=float(np.exp(mu.mean())), OR_per_zP_ci=[float(np.exp(x)) for x in q(mu)],
                OR_per_mmHg=float(np.exp(mu.mean() / 5.0)),
                prediction_interval=q(pred_new), prob_new_study_positive=float((pred_new > 0).mean()),
                beta_by_stratum={s: dict(mean=float(flat['beta'][:, i].mean()), ci=q(flat['beta'][:, i])) for i, s in enumerate(strata)},
                mu_draws=(np.asarray(mu)[np.linspace(0, len(mu) - 1, keep_draws).astype(int)].tolist() if keep_draws else None),
                rhat=dict(mu_beta=float(summ['mu_beta']['r_hat']), tau_beta=float(summ['tau_beta']['r_hat']),
                          max=float(max(np.max(summ[k]['r_hat']) for k in summ))),
                ess=dict(mu_beta=float(summ['mu_beta']['n_eff']), tau_beta=float(summ['tau_beta']['n_eff']),
                         min=float(min(np.min(summ[k]['n_eff']) for k in summ))),
                divergences=div, ppc_inside=int(sum(p['inside'] for p in ppc)), ppc_total=len(ppc), ppc=ppc)


def prior_sensitivity(d, **kw):
    return [fit(d, tau_scale=t, **kw) for t in PRIORS]


def overlap_sets(d, **kw):
    """Kyoto series overlap in time. Refit keeping one Kyoto publication at a time, and one cohort per overlap set."""
    out = []
    ky = sorted(d[d.centre_id == 'Kyoto'].study.unique())
    for keep in ky:
        sub = d[(d.centre_id != 'Kyoto') | (d.study == keep)]
        r = fit(sub, **kw); r['kept_from_Kyoto'] = keep; out.append(r)
    for keep in ky:                                            # one cohort per overlap_set: also collapse Ghent/Kaohsiung
        sub = d[(d.centre_id != 'Kyoto') | (d.study == keep)]
        for cset in ('Ghent 1999-2004', 'Kaohsiung 2007-2012'):
            g = sub[sub.overlap_set == cset]
            if g.study.nunique() > 1: sub = sub[(sub.overlap_set != cset) | (sub.study == sorted(g.study.unique())[0])]
        r = fit(sub, **kw); r['kept_from_Kyoto'] = keep + ' (one cohort per overlap set)'; out.append(r)
    return out


def recovery(d, mu_true=(0.0, 0.5, 1.0, 2.0), tau_true=0.3, reps=40, seed=0, **kw):
    """Can this design recover the effect? Keep the real sizes, gradients and strata, simulate outcomes with a
    known mu_beta and tau_beta, refit, and report bias, coverage and the power to exclude zero."""
    d0, strata, S, z_c = _prep(d, kw.get('by_outcome', True), kw.get('by_centre', False))
    n = d0.n.values.astype(int); rng = np.random.default_rng(seed)
    p_obs = np.clip(d0.events.values / n, 1e-3, 1 - 1e-3)
    base = float(np.log(p_obs / (1 - p_obs)).mean())          # logit of the observed rates, not log
    rows = []
    for mt in mu_true:
        est, cov, pos = [], [], []
        for r in range(reps):
            b = mt + tau_true * rng.standard_normal(len(strata))
            p = 1 / (1 + np.exp(-(base + b[S] * z_c)))
            ev = rng.binomial(n, p)
            sim = d0.copy(); sim['events'] = ev
            f = fit(sim, warmup=600, samples=800, chains=2, seed=r, **kw)
            est.append(f['mu_beta']); cov.append(f['mu_beta_ci'][0] <= mt <= f['mu_beta_ci'][1]); pos.append(f['mu_beta_ci'][0] > 0)
        rows.append(dict(mu_true=mt, mean_estimate=float(np.mean(est)), bias=float(np.mean(est) - mt),
                         coverage=float(np.mean(cov)), power_to_exclude_zero=float(np.mean(pos)), reps=reps))
    return rows


def _conditional_loglik(ev, n, dz, beta):
    """Log of the conditional likelihood of how the events split between two groups, given their total.
    The study's own level cancels exactly (Fisher non-central hypergeometric with log odds ratio beta*dz),
    so this tests the slope without using the held-out centre's outcomes to fit an intercept."""
    a, b = int(ev[0]), int(ev[1]); m1, m2 = int(n[0]), int(n[1]); t = a + b
    psi = np.exp(beta * dz)
    ks = np.arange(max(0, t - m2), min(m1, t) + 1)
    from scipy.special import gammaln
    logw = (gammaln(m1 + 1) - gammaln(ks + 1) - gammaln(m1 - ks + 1)
            + gammaln(m2 + 1) - gammaln(t - ks + 1) - gammaln(m2 - t + ks + 1) + ks * np.log(max(psi, 1e-12)))
    return float(logw[list(ks).index(a)] - (np.max(logw) + np.log(np.sum(np.exp(logw - np.max(logw))))))


def conditional_iecv(d, by='centre_id', **kw):
    """Internal-external cross-validation of the slope alone. Leave out one centre, estimate mu_beta on the rest,
    and score the within-centre contrast of the held-out centre with the conditional likelihood, which removes its
    intercept. Compared with beta = 0, so a positive difference means the slope learned elsewhere predicts the
    contrast better than no association at all. No parameter is fitted on the held-out outcomes."""
    out = []
    for held in sorted(d[by].dropna().unique()):
        train, test = d[d[by] != held], d[d[by] == held]
        tr, strata, S, z_c = _prep(train, kw.get('by_outcome', True), kw.get('by_centre', False))
        te, strata_t, _, _ = _prep(test, kw.get('by_outcome', True), kw.get('by_centre', False))
        if len(tr) < 4 or len(te) < 2: continue
        f = fit(train, warmup=800, samples=1200, chains=2, **kw); mu = f['mu_beta']
        rows = []
        for st, g in te.groupby('stratum'):
            if len(g) != 2: continue
            g = g.sort_values('zP'); dz = float(g.zP.iloc[1] - g.zP.iloc[0])
            ev, n = g.events.values[::-1], g.n.values[::-1]      # events of the higher-zP group first
            ll_m = _conditional_loglik(ev, n, dz, mu); ll_0 = _conditional_loglik(ev, n, dz, 0.0)
            own = _prep(g, kw.get('by_outcome', True))
            rows.append(dict(stratum=st, dzP=dz, loglik_model=ll_m, loglik_null=ll_0, gain=ll_m - ll_0,
                             observed_high_zP_rate=float(g.pct.iloc[1]), observed_low_zP_rate=float(g.pct.iloc[0]),
                             direction_correct=bool(g.pct.iloc[1] >= g.pct.iloc[0])))
        if rows:
            out.append(dict(held_out=held, mu_from_others=mu, mu_ci=f['mu_beta_ci'], strata=rows,
                            total_gain=float(sum(r['gain'] for r in rows)),
                            all_directions_correct=bool(all(r['direction_correct'] for r in rows))))
    return out


def monte_carlo_measurement(d, draws=200, seed=0, **kw):
    """Propagate the uncertainty of the indices. Each group's zP is redrawn from what the paper actually reports:
    the standard error of a reported mean where it exists, a plausible spread otherwise, an imputed CVP sampled
    over 3-9 mmHg, and a cut-off group's value sampled over its side of the cut-off. The model is refitted on each
    draw and the posteriors are pooled, so measurement and statistical uncertainty appear together."""
    rng = np.random.default_rng(seed); mus = []
    for _ in range(draws):
        e = d.copy()
        sd = e.gradient_sd.fillna(3.5).values                     # reported SD where available
        se = sd / np.sqrt(e.n.values.astype(float))               # error of each group's mean
        shift = rng.normal(0, se)
        # The assumed CVP is a property of the study, not of the group: draw one value per study over 3-9 mmHg.
        # zP = (PVP - CVP)/5, so a CVP of 3 raises the gradient by 2 and a CVP of 9 lowers it by 4.
        cvp_imp = e.CVP_basis.eq('imputed').values
        per_study = {st: 5.0 - rng.uniform(3.0, 9.0) for st in e.study.unique()}
        shift = shift + np.where(cvp_imp, e.study.map(per_study).values, 0.0)
        e['zP'] = e.zP + shift / 5.0
        f = fit(e, warmup=400, samples=600, chains=2, quiet=True, keep_draws=200, **kw)
        mus.extend(f['mu_draws'])                 # pool the posteriors, not their means, so the reported
    mus = np.array(mus)                           # interval carries measurement AND statistical uncertainty
    return dict(draws=int(draws), mu_beta=float(mus.mean()),
                mu_beta_ci=[float(np.percentile(mus, 2.5)), float(np.percentile(mus, 97.5))],
                prob_positive=float((mus > 0).mean()))


def spec_curve(d, quick=True):
    """One table instead of many scattered sensitivities: the same quantity, mu_beta, under every reasonable
    specification. The message should not depend on any single one of them."""
    kw = dict(warmup=800, samples=1200, chains=2) if quick else {}
    prim = d[d.outcome_type == 'SFSS_or_dysfunction']; sec = d[d.outcome_type != 'SFSS_or_dysfunction']
    specs = [('primary outcome (SFSS or early dysfunction)', prim, {}),
             ('secondary outcome (mortality or graft loss)', sec, {}),
             ('all outcomes', d, {}),
             ('gradients reported or derived only', d[~d.zP_imputed], {}),
             ('no imputed CVP', d[~d.CVP_basis.eq('imputed')], {}),
             ('events reported only (no counts from percentages)', d[d.events_basis.eq('reported')], {}),
             ('one Kyoto publication (Ogura 2010)', d[(d.centre_id != 'Kyoto') | (d.study == 'Ogura 2010')], {}),
             ('one Kyoto publication (Uemura 2016)', d[(d.centre_id != 'Kyoto') | (d.study == 'Uemura 2016')], {}),
             ('one Kyoto publication (Yagi 2005)', d[(d.centre_id != 'Kyoto') | (d.study == 'Yagi 2005')], {}),
             ('tau prior HalfNormal(0.25)', d, dict(tau_scale=0.25)),
             ('tau prior HalfNormal(1.0)', d, dict(tau_scale=1.0)),
             ('one stratum per centre, not per publication', d, dict(by_centre=True))]
    for cvp in (3.0, 7.0, 9.0):
        e = d.copy(); m = e.CVP_basis.eq('imputed')
        e.loc[m, 'zP'] = e.loc[m, 'zP'] + (5.0 - cvp) / 5.0
        specs.append((f'assumed CVP {cvp:.0f} mmHg', e, {}))
    out = []
    for name, sub, extra in specs:
        try:
            s2 = _prep(sub, True, extra.get('by_centre', False))[0]
            if len(s2) < 4: out.append(dict(spec=name, groups=int(len(s2)), mu_beta=None)); continue
            f = fit(sub, **{**kw, **extra})
            out.append(dict(spec=name, groups=f['groups'], strata=len(f['strata']), events=f['events'],
                            mu_beta=f['mu_beta'], mu_beta_ci=f['mu_beta_ci'], prob_positive=f['prob_mu_positive'],
                            OR_per_zP=f['OR_per_zP']))
        except Exception as e:                                   # a specification that leaves too little data
            out.append(dict(spec=name, groups=None, mu_beta=None, error=str(e)[:80]))
    return out


def main():
    """Fit and cache everything the primary analysis needs (a few minutes)."""
    import json, os, indices
    root = os.path.join(os.path.dirname(__file__), '..'); out = f'{root}/results'; os.makedirs(out, exist_ok=True)
    d = indices.compute(pd.read_csv(f'{root}/data/series.tsv', sep='\t'))
    c = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]
    pri = c[c.outcome_type == 'SFSS_or_dysfunction']; sec = c[c.outcome_type != 'SFSS_or_dysfunction']
    main_ = {'primary (SFSS or early dysfunction)': fit(pri), 'secondary (mortality or graft loss)': fit(sec),
             'exploratory (all outcomes)': fit(c)}
    json.dump(main_, open(f'{out}/bayes_main.json', 'w'), indent=1)
    for name, res in (('bayes_priors', prior_sensitivity(pri)), ('bayes_overlap', overlap_sets(c)),
                      ('bayes_iecv', conditional_iecv(c)), ('spec_curve', spec_curve(c))):
        json.dump(res, open(f'{out}/{name}.json', 'w'), indent=1); print('wrote', name)
    json.dump(monte_carlo_measurement(c, draws=40), open(f'{out}/bayes_measurement.json', 'w'), indent=1); print('wrote bayes_measurement')
    json.dump(recovery(c, mu_true=(0.0, 1.0, 2.0), reps=10), open(f'{out}/bayes_recovery.json', 'w'), indent=1); print('wrote bayes_recovery')
    r = main_['primary (SFSS or early dysfunction)']
    print(f"primary: mu {r['mu_beta']:+.2f} ({r['mu_beta_ci'][0]:+.2f} to {r['mu_beta_ci'][1]:+.2f}), "
          f"P(mu>0) = {r['prob_mu_positive']:.3f}, R-hat {r['rhat']['max']:.3f}, divergences {r['divergences']}")


if __name__ == '__main__':
    main()
