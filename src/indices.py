"""LiPNet indices: the three loads a partial liver graft carries, and the models fitted to them.

    zF = normalised portal flow load      graft PVF per 100 g / donor PVF per 100 g
    zP = normalised portal pressure load  (PVP - CVP) / 5 mmHg
    zR = zP / zF = R_graft / R_donor      normalised resistance load

so that

    zP = zF * zR        pressure load = flow load x resistance load

All three equal 1 in a healthy donor. The identity is why flow and pressure stop being interchangeable: they
coincide only while zR = 1, and any manoeuvre that enlarges the venous outflow lowers zR, which is how a graft can
carry three or four times the donor's flow without the pressure that would normally come with it.

    python src/indices.py     -> results/series_with_indices.tsv, results/within_study_fit.json,
                                 results/meta_slope.json, results/sensitivity.tsv, results/pooled_fit.json

Indices (the calculator and plots.py use these same functions)
    zF = graft PVF per 100 g / donor PVF per 100 g   (donor reference from the same series when measured, else 90)
    zP = (PVP - CVP) / 5 mmHg                          (CVP = 5 when not reported; provenance in the *_basis columns)
The stored table holds raw values and their provenance only; the indices are always recomputed here.

Risk model. A single pooled logistic model with a common intercept does not hold across centres: at the same
gradient the reported incidence ranges from 0% (Ishizaki 2012) to 17% (Uemura 2016), because the baseline level
depends on recipient severity, outcome definition and technique. The estimable quantity is the WITHIN-STUDY
slope, fitted with one intercept per study and a common slope:

    logit(p) = alpha_study + beta * zP

The primary estimate of the effect is the hierarchical model in src/bayes.py (NUTS, zP centred within study);
this stratified fit provides the per-series levels used by the calculator and the descriptive views; meta_slope() repeats it as a random-effects meta-analysis of the
per-stratum slopes (sensitivity, and a measure of the heterogeneity), overdispersion() checks for extra-binomial
spread and centred() is a descriptive view on a common scale. beta is the transferable parameter (the change in odds per unit of zP, i.e. per 5 mmHg of gradient); alpha is
not transferable and must be supplied by the user as their own baseline rate. The pooled common-intercept fit is
still computed, for the record, and reported as not robust.

Main analysis: groups whose haemodynamic value is measured or derived (in_main_analysis = True). Excluded and
reported in sensitivity: groups defined only by a cut-off (Vasavada 2014, Yao 2018, Kanetkar 2017) and the second
partition of a cohort already represented (the Wang 2014 pressure groups are the same patients as its splenectomy
groups). Intercepts are per study AND outcome definition, so a study reporting two outcomes does not share one.
"""
import json, os
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr

ROOT = os.path.join(os.path.dirname(__file__), '..')
DATA = f'{ROOT}/data/series.tsv'; OUT = f'{ROOT}/results'
NORMAL_GRADIENT = 5.0; DEFAULT_CVP = 5.0; DEFAULT_DONOR_REF = 90.0
IMPUTED_BASES = ('imputed', 'threshold', 'group mean')

def zF(pvf_per_100g, donor_ref=DEFAULT_DONOR_REF):
    return pvf_per_100g / donor_ref

def zR(zP_value, zF_value):
    """Normalised resistance load, zR = zP / zF = R_graft / R_donor, the third term of zP = zF * zR."""
    return zP_value / zF_value


def zP(pvp=None, cvp=DEFAULT_CVP, gradient=None):
    """zP from the portocaval gradient. Pass `gradient` when the paper reports it directly; otherwise it is
    PVP - CVP, with CVP defaulting to 5 mmHg."""
    g = (pvp - cvp) if gradient is None else gradient
    return g / NORMAL_GRADIENT

def compute(d):
    """add zF, zP and provenance flags to a series table (raw columns only).
    Any value filled in here (missing CVP -> 5, missing donor reference -> 90) is recorded in the *_filled columns
    and counts as imputed, whatever the *_basis label says; a label that contradicts the data raises an error."""
    d = d.copy()
    has_g = d.gradient_mmHg.notna()
    cvp_filled = d.PVP.notna() & d.CVP.isna() & ~has_g          # CVP defaulted to 5 to build the gradient
    ref_filled = d.PVF_per_100g.notna() & d.donor_PVF_per_100g_ref.isna()
    d['CVP_filled'] = cvp_filled; d['donor_ref_filled'] = ref_filled
    bad = d[cvp_filled & d.CVP_basis.eq('reported')]
    if len(bad): raise ValueError('CVP missing but labelled as reported: ' + ', '.join(bad.study + ' / ' + bad.group))
    bad = d[d.CVP.notna() & d.CVP_basis.eq('not applicable')]
    if len(bad): raise ValueError('CVP present but labelled not applicable: ' + ', '.join(bad.study + ' / ' + bad.group))
    bad = d[has_g & d.gradient_basis.eq('not applicable')]
    if len(bad): raise ValueError('gradient present but labelled not applicable: ' + ', '.join(bad.study + ' / ' + bad.group))
    d['zF'] = zF(d.PVF_per_100g, d.donor_PVF_per_100g_ref.fillna(DEFAULT_DONOR_REF))
    d['zP'] = np.where(has_g, zP(gradient=d.gradient_mmHg), zP(d.PVP, d.CVP.fillna(DEFAULT_CVP)))
    d['zR'] = d.zP / d.zF                                      # resistance load, where both indices exist
    d['zP_imputed'] = np.where(has_g, d.gradient_basis.isin(IMPUTED_BASES),
                               d[['PVP_basis', 'CVP_basis']].isin(IMPUTED_BASES).any(axis=1) | cvp_filled)
    d['zF_imputed'] = d.PVF_basis.isin(IMPUTED_BASES) | ref_filled
    d['any_imputed'] = d.zP_imputed | d.zF_imputed
    if 'in_main_analysis' not in d: d['in_main_analysis'] = ~d.PVF_basis.isin(['threshold', 'group mean'])
    return d

# ---------------------------------------------------------------- risk models
def _nll_fixed(par, S, z, ev, n, k):
    p = np.clip(1 / (1 + np.exp(-(par[:k][S] + par[-1] * z))), 1e-9, 1 - 1e-9)
    return -np.sum(ev * np.log(p) + (n - ev) * np.log(1 - p))

def strata(s):
    """One stratum per study AND outcome definition: a study reporting two outcomes does not share an intercept."""
    return s.study.astype(str) + ' / ' + s.outcome_type.astype(str)

def contributing(s):
    """Keep the strata that contribute a within-stratum contrast (>=2 groups with different zP)."""
    return s.groupby(strata(s).values).filter(lambda x: len(x) >= 2 and x.zP.nunique() > 1)

def within_study_fit(s, boot=2000, seed=0, _inner=False):
    """One intercept per stratum (study x outcome), common slope. s must have study, outcome_type, n, events, zP."""
    st = strata(s); studies = sorted(st.unique()); k = len(studies)
    S = np.array([studies.index(x) for x in st]); z = s.zP.values
    ev = s.events.values.astype(float); n = s.n.values.astype(float)
    x0 = np.r_[np.full(k, -2.0), 1.0]
    fit = minimize(_nll_fixed, x0, args=(S, z, ev, n, k), method='BFGS')
    beta = float(fit.x[-1]); alphas = {st: float(a) for st, a in zip(studies, fit.x[:k])}
    n_failed = 0; zgrid = np.round(np.arange(0.6, 4.001, 0.05), 3); preds = []
    if boot:
        rng = np.random.default_rng(seed); B = []
        for _ in range(boot):
            evb = rng.binomial(n.astype(int), ev / n)
            r = minimize(_nll_fixed, fit.x, args=(S, z, evb, n, k), method='BFGS')
            if not r.success:                                    # retry once with a derivative-free method
                r = minimize(_nll_fixed, fit.x, args=(S, z, evb, n, k), method='Nelder-Mead',
                             options=dict(maxiter=20000, fatol=1e-10, xatol=1e-8))
            if r.success:
                B.append(r.x[-1])
                preds.append(1 / (1 + np.exp(-(r.x[:k][:, None] + r.x[-1] * zgrid[None, :]))))
            else: n_failed += 1
        B = np.array(B); lo, hi = np.percentile(B, [2.5, 97.5])
    else:
        lo = hi = beta
    # profile-likelihood interval as an independent check on the bootstrap
    prof_ci = None
    if boot:
        l0 = fit.fun
        f = lambda b: minimize(lambda a: _nll_fixed(np.r_[a, b], S, z, ev, n, k), fit.x[:k]).fun - l0 - 1.92
        try:
            from scipy.optimize import brentq
            prof_ci = [float(brentq(f, beta - 3, beta - 1e-3)), float(brentq(f, beta + 1e-3, beta + 3))]
        except Exception:
            prof_ci = None
    loo, by_outcome = {}, {}
    if not _inner:
        for one in sorted(s.study.unique()):
            sub = contributing(s[s.study != one])
            if len(sub) >= 2: loo[one] = float(within_study_fit(sub, boot=0, _inner=True)['beta'])
        for oc, g in s.groupby('outcome_type'):
            g = contributing(g)
            if len(g) >= 2: by_outcome[oc] = float(within_study_fit(g, boot=0, _inner=True)['beta'])
    return dict(beta=beta, beta_ci=[float(lo), float(hi)], OR_per_zP=float(np.exp(beta)),
                OR_per_zP_ci=[float(np.exp(lo)), float(np.exp(hi))],
                OR_per_mmHg=float(np.exp(beta / NORMAL_GRADIENT)),
                OR_per_mmHg_ci=[float(np.exp(lo / NORMAL_GRADIENT)), float(np.exp(hi / NORMAL_GRADIENT))],
                profile_ci=prof_ci, bootstrap_failures=int(n_failed), by_outcome=by_outcome,
                studies=studies, alphas=alphas, groups=int(len(s)), patients=int(n.sum()), events=int(ev.sum()),
                zgrid=zgrid.tolist(),
                pred_band={st: dict(lo=np.percentile(np.array(preds)[:, i, :], 2.5, axis=0).round(5).tolist(),
                                    hi=np.percentile(np.array(preds)[:, i, :], 97.5, axis=0).round(5).tolist())
                           for i, st in enumerate(studies)} if preds else {},
                zP_min=float(z.min()), zP_max=float(z.max()), normal_gradient_mmHg=NORMAL_GRADIENT,
                default_cvp_mmHg=DEFAULT_CVP, leave_one_study_out=loo,
                study_groups=[dict(study=r.study, group=r.group, n=int(r.n), events=int(r.events),
                                   zP=round(float(r.zP), 3), pct=float(r.pct), outcome=r.outcome, stratum=f'{r.study} / {r.outcome_type}',
                                   outcome_type=r.outcome_type, zP_imputed=bool(r.zP_imputed)) for _, r in s.iterrows()])

def meta_slope(s):
    """Exploratory sensitivity analysis: random-effects meta-analysis of the within-stratum slopes.

    DerSimonian-Laird for tau2 with a Hartung-Knapp interval, which is the one reported: with five strata a normal
    interval on the random-effects standard error is too narrow. The primary estimate remains within_study_fit();
    I2 and tau2 here describe consistency among five estimates and have little power to detect heterogeneity.

    Each stratum contributes its own logistic slope of the outcome on zP with its standard error; the slopes are
    pooled with an additive between-study variance tau2. This is the standard way of combining grouped-exposure
    data across studies, and unlike a single fit with stratum intercepts it *measures* the heterogeneity (tau2, I2)
    instead of assuming it away, and does not spend degrees of freedom estimating the intercepts.
    """
    rows = []
    for st, g in s.groupby(strata(s).values):
        z = g.zP.values; ev = g.events.values.astype(float); n = g.n.values.astype(float)
        f = lambda par: _nll_fixed(par, np.zeros(len(z), int), z, ev, n, 1)
        r = minimize(f, [-3.0, 1.0], method='BFGS')
        se = float(np.sqrt(np.diag(r.hess_inv))[1])
        rows.append(dict(stratum=st, groups=int(len(g)), patients=int(n.sum()), events=int(ev.sum()),
                         beta=float(r.x[1]), se=se, zP_span=float(z.max() - z.min())))
    p = pd.DataFrame(rows); w = 1 / p.se ** 2
    mu_fe = float((w * p.beta).sum() / w.sum()); Q = float((w * (p.beta - mu_fe) ** 2).sum()); k = len(p)
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
    w2 = 1 / (p.se ** 2 + tau2); mu = float((w2 * p.beta).sum() / w2.sum()); se = float(np.sqrt(1 / w2.sum()))
    # Hartung-Knapp: scale the variance by the observed dispersion of the estimates and use a t quantile.
    # With five strata the normal interval is too narrow, so this is the one reported.
    from scipy.stats import t as _t
    q = float((w2 * (p.beta - mu) ** 2).sum() / (k - 1)) if k > 1 else 1.0
    q = max(q, 1.0)                                   # ad hoc truncation: never narrower than the normal interval
    se_hk = float(np.sqrt(q / w2.sum())); tcrit = float(_t.ppf(0.975, k - 1)) if k > 1 else 1.96
    lo, hi = mu - tcrit * se_hk, mu + tcrit * se_hk
    I2 = float(max(0.0, 100 * (Q - (k - 1)) / Q)) if Q > 0 else 0.0
    from scipy.stats import chi2 as _chi2
    return dict(beta=mu, se=se, se_hartung_knapp=se_hk, beta_ci=[lo, hi], beta_ci_normal=[mu - 1.96 * se, mu + 1.96 * se], beta_fixed_effect=mu_fe, tau2=float(tau2), Q=Q, Q_df=k - 1,
                Q_p=float(1 - _chi2.cdf(Q, k - 1)) if k > 1 else float('nan'), I2=I2,
                OR_per_zP=float(np.exp(mu)), OR_per_zP_ci=[float(np.exp(lo)), float(np.exp(hi))],
                OR_per_mmHg=float(np.exp(mu / NORMAL_GRADIENT)),
                OR_per_mmHg_ci=[float(np.exp(lo / NORMAL_GRADIENT)), float(np.exp(hi / NORMAL_GRADIENT))],
                strata=rows)

def centred(s):
    """Descriptive view with the cohort level removed (not the primary estimator).

    Each group's empirical log-odds (Haldane-corrected) and zP are centred on the mean of its own stratum, so every
    group lands on one scale: outcome odds and gradient, each relative to that cohort's own average. The weighted
    regression through the origin gives a slope of the same sign and order as the stratified logistic fit, but it is
    NOT the same estimator: it works on transformed proportions, and its interval ignores the dependence induced by
    centring observations of one study. Use it to show the data on a common scale; take the effect from
    within_study_fit(), with meta_slope() as sensitivity.
    """
    d = s.copy()
    d['log_odds'] = np.log((d.events + 0.5) / (d.n - d.events + 0.5))
    d['var_log_odds'] = 1 / (d.events + 0.5) + 1 / (d.n - d.events + 0.5)
    g = strata(d).values
    d['log_odds_c'] = d.log_odds - d.groupby(g).log_odds.transform('mean')
    d['zP_c'] = d.zP - d.groupby(g).zP.transform('mean')
    d['OR_vs_own_cohort'] = np.exp(d.log_odds_c)
    w = 1 / d.var_log_odds
    beta = float((w * d.zP_c * d.log_odds_c).sum() / (w * d.zP_c ** 2).sum())
    se = float(np.sqrt(1 / (w * d.zP_c ** 2).sum()))
    r2 = float(1 - (w * (d.log_odds_c - beta * d.zP_c) ** 2).sum() / (w * d.log_odds_c ** 2).sum())
    lo, hi = beta - 1.96 * se, beta + 1.96 * se
    return dict(beta=beta, se=se, beta_ci=[lo, hi], r2_weighted=r2,
                OR_per_zP=float(np.exp(beta)), OR_per_zP_ci=[float(np.exp(lo)), float(np.exp(hi))],
                OR_per_mmHg=float(np.exp(beta / NORMAL_GRADIENT)),
                points=[dict(study=r.study, group=r.group, stratum=f'{r.study} / {r.outcome_type}', n=int(r.n),
                             events=int(r.events), zP_c=round(float(r.zP_c), 4),
                             log_odds_c=round(float(r.log_odds_c), 4), se_log_odds=round(float(np.sqrt(r.var_log_odds)), 4))
                        for _, r in d.iterrows()])

def overdispersion(s):
    """Pearson chi2 / df of the stratified fit: >1 means more spread than binomial sampling explains, and the
    standard error should be scaled by its square root (quasi-binomial)."""
    st = strata(s); studies = sorted(st.unique()); k = len(studies)
    S = np.array([studies.index(x) for x in st]); z = s.zP.values
    ev = s.events.values.astype(float); n = s.n.values.astype(float)
    fit = minimize(_nll_fixed, np.r_[np.full(k, -3.0), 1.0], args=(S, z, ev, n, k), method='BFGS')
    p = 1 / (1 + np.exp(-(fit.x[:k][S] + fit.x[-1] * z)))
    chi2 = float(np.sum((ev - n * p) ** 2 / (n * p * (1 - p)))); df = int(len(z) - k - 1)
    scale = float(np.sqrt(max(1.0, chi2 / df))) if df > 0 else 1.0
    se = float(np.sqrt(np.diag(fit.hess_inv))[-1])
    return dict(chi2=chi2, df=df, ratio=chi2 / df if df > 0 else float('nan'), se_scale=scale,
                beta=float(fit.x[-1]), se_model=se, se_quasi=se * scale,
                beta_ci_quasi=[float(fit.x[-1] - 1.96 * se * scale), float(fit.x[-1] + 1.96 * se * scale)])

def attenuation(s, within_group_sd_mmHg=3.5, draws=400, seed=1):
    """Scenario analysis, not an observed result: the within-group SD is assumed, not reported by the series.
    How much of the slope the sampling error of the group means could cost: each group's zP is perturbed by
    sd/sqrt(n) and the slope re-estimated. It says nothing about a per-patient slope: recovering an individual-level
    relation from group means is not possible here (ecological bias, between-study variation, error in the means)."""
    st = strata(s); studies = sorted(st.unique()); k = len(studies)
    S = np.array([studies.index(x) for x in st]); z = s.zP.values
    ev = s.events.values.astype(float); n = s.n.values.astype(float)
    fit = minimize(_nll_fixed, np.r_[np.full(k, -3.0), 1.0], args=(S, z, ev, n, k), method='BFGS'); b = float(fit.x[-1])
    rng = np.random.default_rng(seed); se_z = (within_group_sd_mmHg / np.sqrt(n)) / NORMAL_GRADIENT
    sims = [minimize(_nll_fixed, fit.x, args=(S, z + rng.normal(0, se_z), ev, n, k), method='BFGS').x[-1] for _ in range(draws)]
    return dict(beta=b, beta_with_group_mean_error=float(np.mean(sims)),
                loss_pct=float(100 * (1 - np.mean(sims) / b)), within_group_sd_mmHg=within_group_sd_mmHg)

def pooled_fit(s, boot=2000, seed=0):
    """Common-intercept logistic (kept for the record; not robust across centres)."""
    z, ev, n = s.zP.values, s.events.values.astype(float), s.n.values.astype(float)
    f = lambda b: _nll_fixed(np.r_[b[0], b[1]], np.zeros(len(z), int), z, ev, n, 1)
    b = minimize(f, [-4.0, 1.5]).x; rng = np.random.default_rng(seed)
    B = np.array([minimize(lambda p, e=rng.binomial(n.astype(int), ev / n): _nll_fixed(np.r_[p[0], p[1]], np.zeros(len(z), int), z, e, n, 1), b).x for _ in range(boot)]) if boot else np.empty((0, 2))
    lo, hi = (np.percentile(B, [2.5, 97.5], axis=0) if boot else (b, b))
    return dict(b0=float(b[0]), b1=float(b[1]), b1_ci=[float(lo[1]), float(hi[1])], groups=int(len(s)),
                patients=int(n.sum()), events=int(ev.sum()))

def recalibrate(outcomes, zP_values, beta, se=True):
    """Recalibration in the large: keep the slope, fit only the intercept for a new cohort.

    outcomes: one 0/1 per patient (grouped counts are not accepted).
    Returns the intercept alpha of logit(p) = alpha + beta*zP for that cohort, so that
    risk = 1 / (1 + exp(-(alpha + beta*zP))). This is the standard first step of prediction-model updating
    (Steyerberg; Vergouwe et al., Stat Med 2017): the coefficients travel, the level does not.
    With only an overall rate p and an average gradient, alpha = logit(p) - beta * mean(zP).
    """
    y = np.asarray(outcomes, float); z = np.asarray(zP_values, float)
    off = beta * z
    f = lambda a: -np.sum(y * np.log(np.clip(1 / (1 + np.exp(-(a[0] + off))), 1e-9, 1 - 1e-9))
                          + (1 - y) * np.log(np.clip(1 - 1 / (1 + np.exp(-(a[0] + off))), 1e-9, 1 - 1e-9)))
    r = minimize(f, [np.log(max(y.mean(), 1e-6) / max(1 - y.mean(), 1e-6))], method='BFGS')
    a = float(r.x[0])
    out = dict(alpha=a, n=int(len(y)), events=int(y.sum()), beta_used=float(beta))
    if se: out['se'] = float(np.sqrt(np.diag(r.hess_inv))[0])
    return out

def baseline_from_rate(rate, mean_zP, beta):
    """Calibration in the large from the two numbers a centre can always state: its overall outcome rate and
    its average zP. Returns the intercept alpha = logit(rate) - beta*mean(zP)."""
    p = np.clip(rate, 1e-6, 1 - 1e-6)
    return float(np.log(p / (1 - p)) - beta * mean_zP)

def risk_after(baseline_rate, delta_zP, beta):
    """Absolute risk implied by a change in zP, given the user's own baseline rate at the starting gradient."""
    p0 = np.clip(baseline_rate, 1e-6, 1 - 1e-6)
    odds = p0 / (1 - p0) * np.exp(beta * delta_zP)
    return odds / (1 + odds)

# ---------------------------------------------------------------- hierarchical model and validation
def cvp_scenarios(d, values=(3.0, 5.0, 7.0, 9.0)):
    """How the effect moves if the assumed CVP is 3, 5, 7 or 9 mmHg instead of 5, for the groups where it was imputed."""
    out = []
    for v in values:
        e = d.copy()
        m = e.CVP.isna() & e.PVP.notna() & e.gradient_mmHg.isna()
        e.loc[m, 'CVP'] = v; e.loc[m, 'CVP_basis'] = 'imputed'
        c = compute(e); c = contributing(c[c.in_main_analysis & c.zP.notna() & c.events.notna()])
        r = within_study_fit(c, boot=0, _inner=True)
        out.append(dict(assumed_CVP=v, beta=r['beta'], OR_per_zP=float(np.exp(r['beta'])), groups=r['groups']))
    return out


# ---------------------------------------------------------------- driver
def main(boot=2000, seed=0):
    os.makedirs(OUT, exist_ok=True)
    d = compute(pd.read_csv(DATA, sep='\t'))
    d.to_csv(f'{OUT}/series_with_indices.tsv', sep='\t', index=False, float_format='%.3f')

    rows = []
    for label, subs in (('main analysis (measured or derived haemodynamics)', {'zP': d[d.in_main_analysis], 'zF': d[d.in_main_analysis]}),
                        ('including groups defined by a cut-off', {'zP': d, 'zF': d}),
                        ('no imputed or filled value in the index itself', {'zP': d[~d.zP_imputed & d.in_main_analysis], 'zF': d[~d.zF_imputed & d.in_main_analysis]})):
        for col in ('zP', 'zF'):
            s = subs[col].dropna(subset=[col])
            if len(s) >= 3:
                rho, p = spearmanr(s[col], s.pct); rows.append(dict(analysis=label, index=col, groups=len(s), spearman_rho=rho, p=p))
                print(f'{label[:46]:46s} {col}: {len(s):2d} groups, Spearman rho = {rho:+.2f}, p = {p:.3f}')

    # within-study model: studies contributing at least two groups with different zP
    m = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]
    contrib = contributing(m)
    res = within_study_fit(contrib, boot, seed)
    json.dump(res, open(f'{OUT}/within_study_fit.json', 'w'), indent=1)
    print(f"\nWithin-study model: logit(p) = alpha_study + {res['beta']:.2f} zP   "
          f"(OR {res['OR_per_zP']:.2f} per unit zP, 95% CI {res['OR_per_zP_ci'][0]:.2f}-{res['OR_per_zP_ci'][1]:.2f}; "
          f"OR {res['OR_per_mmHg']:.2f} per mmHg of gradient)")
    print(f"  {res['groups']} groups from {len(res['studies'])} studies, {res['events']} events / {res['patients']} recipients")
    if res['profile_ci']: print(f"  profile-likelihood 95% CI for the slope: {res['profile_ci'][0]:.2f} to {res['profile_ci'][1]:.2f}"
                                f"   (bootstrap {res['beta_ci'][0]:.2f} to {res['beta_ci'][1]:.2f}; {res['bootstrap_failures']} replicates discarded)")
    print('  leave-one-study-out slopes: ' + ', '.join(f'{k} {v:.2f}' for k, v in res['leave_one_study_out'].items()))
    print('  by outcome: ' + ', '.join(f'{k} {v:.2f}' for k, v in res['by_outcome'].items()))
    for st, a in res['alphas'].items():
        p_at2 = 1 / (1 + np.exp(-(a + res['beta'] * 2)))
        print(f"    baseline {st:48s} alpha {a:+.2f}  -> risk at zP = 2: {100*p_at2:.1f}%")

    meta = meta_slope(contrib); od = overdispersion(contrib); att = attenuation(contrib); cen = centred(contrib)
    json.dump(dict(meta=meta, overdispersion=od, attenuation=att, centred=cen), open(f'{OUT}/meta_slope.json', 'w'), indent=1)
    print(f"\nCohort level removed by centring within stratum: slope {cen['beta']:.2f} +- {cen['se']:.2f}, "
          f"OR {cen['OR_per_zP']:.2f} per unit zP ({cen['OR_per_zP_ci'][0]:.2f}-{cen['OR_per_zP_ci'][1]:.2f}), "
          f"{cen['OR_per_mmHg']:.2f} per mmHg; weighted R2 of the common line = {cen['r2_weighted']:.2f}")
    print(f"\nSensitivity, random-effects meta-analysis (Hartung-Knapp): {meta['beta']:.2f} (95% CI {meta['beta_ci'][0]:.2f} to {meta['beta_ci'][1]:.2f}), "
          f"OR {meta['OR_per_zP']:.2f} per unit zP ({meta['OR_per_zP_ci'][0]:.2f}-{meta['OR_per_zP_ci'][1]:.2f}), {meta['OR_per_mmHg']:.2f} per mmHg")
    print(f"  heterogeneity: tau2 = {meta['tau2']:.3f}, Q = {meta['Q']:.2f} on {meta['Q_df']} df (p = {meta['Q_p']:.2f}), I2 = {meta['I2']:.0f}% "
          f"(five estimates: little power to detect heterogeneity)")
    for r in meta['strata']: print(f"    {r['stratum']:48s} slope {r['beta']:+.2f} (SE {r['se']:.2f}), zP span {r['zP_span']:.2f}, {r['events']}/{r['patients']}")
    print(f"  overdispersion: Pearson chi2/df = {od['ratio']:.2f} on {od['df']} df -> SE scale {od['se_scale']:.2f}; quasi-binomial CI {od['beta_ci_quasi'][0]:.2f} to {od['beta_ci_quasi'][1]:.2f}")
    print(f"  scenario analysis: assuming a within-group SD of {att['within_group_sd_mmHg']} mmHg, perturbing the reported "
          f"group means changes the slope by about {att['loss_pct']:.0f}%")

    cvps = cvp_scenarios(d)
    json.dump(dict(cvp_scenarios=cvps), open(f'{OUT}/cvp_scenarios.json', 'w'), indent=1)
    slopes = ", ".join("%.2f" % x["beta"] for x in cvps)
    print("\nAssumed CVP of 3, 5, 7 or 9 mmHg: slope " + slopes +
          "  (a constant shift within a study is absorbed by its intercept)")

    sf = m[(m.outcome_type == 'SFSS_or_dysfunction')]
    pf = pooled_fit(sf, boot, seed); json.dump(pf, open(f'{OUT}/pooled_fit.json', 'w'), indent=1)
    print(f"\nCommon-intercept pooled fit (not robust): logit = {pf['b0']:.2f} + {pf['b1']:.2f} zP "
          f"(slope 95% CI {pf['b1_ci'][0]:.2f} to {pf['b1_ci'][1]:.2f}), {pf['groups']} groups")
    rows.append(dict(analysis='common-intercept pooled fit (SFSS groups)', index='zP', groups=pf['groups'],
                     spearman_rho=np.nan, p=np.nan, b0=pf['b0'], b1=pf['b1']))
    rows.append(dict(analysis='within-study fit (common slope, stratum intercepts)', index='zP', groups=res['groups'],
                     spearman_rho=np.nan, p=np.nan, b1=res['beta'], b1_lo=res['beta_ci'][0], b1_hi=res['beta_ci'][1]))
    rows.append(dict(analysis='random-effects meta-analysis of stratum slopes', index='zP', groups=res['groups'],
                     spearman_rho=np.nan, p=np.nan, b1=meta['beta'], b1_lo=meta['beta_ci'][0], b1_hi=meta['beta_ci'][1],
                     tau2=meta['tau2'], I2=meta['I2'], Q=meta['Q'], Q_df=meta['Q_df'], Q_p=meta['Q_p']))
    rows.append(dict(analysis='overdispersion (Pearson chi2/df) and quasi-binomial CI', index='zP', groups=res['groups'],
                     spearman_rho=np.nan, p=np.nan, b1=od['beta'], b1_lo=od['beta_ci_quasi'][0], b1_hi=od['beta_ci_quasi'][1],
                     chi2_over_df=od['ratio']))
    rows.append(dict(analysis='centred within stratum (cohort level removed)', index='zP', groups=res['groups'],
                     spearman_rho=np.nan, p=np.nan, b1=cen['beta'], b1_lo=cen['beta_ci'][0], b1_hi=cen['beta_ci'][1],
                     r2_weighted=cen['r2_weighted']))
    rows.append(dict(analysis='regression dilution (sampling error of the group means)', index='zP',
                     groups=res['groups'], spearman_rho=np.nan, p=np.nan, b1=att['beta_with_group_mean_error']))
    pd.DataFrame(rows).to_csv(f'{OUT}/sensitivity.tsv', sep='\t', index=False, float_format='%.4f')

if __name__ == '__main__':
    main()
