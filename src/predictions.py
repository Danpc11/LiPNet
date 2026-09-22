"""Test the predictions of the perfusion-network model against independent data.

    python src/predictions.py        -> results/predictions.json, printed summary

The model (src/model/) describes the liver as a dissipation-optimal perfusion network at fixed vascular maintenance
cost. It fixes the wall shear set-point tau0 = sqrt(D*/C), the scaling of the optimum with the number of lobules,
and the local relation dP = f * R, from which the two indices follow. It says nothing about recipient severity,
donor age or steatosis. Four predictions follow, and each is tested here against data the model was not fitted to:

  P1  Comparative anatomy. With portal pressure invariant across mammals, the set-point shear must be invariant with
      body mass and the vascular mass must scale close to the observed exponent of hepatic blood volume (0.86).
  P2  Flow and pressure are the same variable only at donor outflow resistance. Grafts whose outflow was
      reconstructed must sit below the diagonal (zP/zF < 1) and tolerate high flow.
  P3  Within a cohort, the outcome must increase with zP, with a common slope and a centre-specific level.
  P4  The clinical thresholds of three different fields must coincide once normalised: HVPG 10 mmHg in cirrhosis,
      PVP 15 mmHg at CVP 5 in a graft, and the limit of safe resection all fall at zP = 2.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import indices
from model import scaling as sc
import json as _json

ROOT = os.path.join(os.path.dirname(__file__), '..'); OUT = f'{ROOT}/results'


def p1_allometry():
    """Set-point and vascular-mass exponents under an invariant portal pressure, for the measured flow and lobule
    scalings (Kruepunga 2019: portal flow ~ M^0.77, portocentral distance ~ M^0.10, portal pressure 6-11 mmHg)."""
    gs = np.arange(8, 19); rows = []
    for f, flab in ((0.77, 'portal flow ~ M^0.77'), (0.88, 'total hepatic flow ~ M^0.88')):
        for l, n in ((0.0, 0.89), (0.10, 0.59)):
            for b in (1.0, 0.9, 0.85, 0.8, 0.75, 2 / 3, 0.5):
                R = [sc.canopy_optimum(g, b, 3) for g in gs]; N = np.array([r['N'] for r in R])
                eD = np.polyfit(np.log(N), np.log([r['D'] for r in R]), 1)[0]
                eV = np.polyfit(np.log(N), np.log([r['V'] for r in R]), 1)[0]
                c = b * (f + eD * n + 3 * l) / 2
                rows.append(dict(flow=flab, lobule_exponent=l, b=b,
                                 tau0_exponent=(2 * f - (2 + b) * c / b + eD * n + 3 * l) / 2,
                                 vascular_mass_exponent=c / b + eV * n))
    t = pd.DataFrame(rows); t.to_csv(f'{OUT}/prediction_allometry.tsv', sep='\t', index=False, float_format='%.4f')
    # reported over the whole prespecified range of b, with no value selected for agreeing with the observation
    obs_lo, obs_hi = 0.80, 0.92          # hepatic blood volume exponent with its reported spread
    return dict(b_range=[float(t.b.min()), float(t.b.max())],
                tau0_exponent_range=[float(t.tau0_exponent.min()), float(t.tau0_exponent.max())],
                tau0_exponent_iqr=[float(t.tau0_exponent.quantile(.25)), float(t.tau0_exponent.quantile(.75))],
                tau0_invariant=bool(t.tau0_exponent.abs().max() < 0.1),
                vascular_mass_exponent_range=[float(t.vascular_mass_exponent.min()), float(t.vascular_mass_exponent.max())],
                observed_vascular_mass_exponent=0.86, observed_range=[obs_lo, obs_hi],
                fraction_compatible=float(((t.vascular_mass_exponent >= obs_lo) & (t.vascular_mass_exponent <= obs_hi)).mean()),
                note='no value of b was selected for fitting the observation; the whole prespecified range is reported')


def p2_discordance(d, draws=4000, seed=0):
    """Mechanistic contrast, kept apart from outcome. The ratio zP/zF is the graft's *relative effective
    resistance*: what the pressure per unit of flow is, compared with the donor. The model predicts it below one
    whenever the outflow was enlarged. Surgical configuration was taken from the methods of each paper, before
    looking at its results. Uncertainty comes from the reported spread of the pressures and flows; the direction is
    tested with an exact sign test, and outcome is reported separately, not as evidence for the mechanism."""
    from scipy.stats import binomtest
    b = d.dropna(subset=['zP', 'zF']).copy(); b['ratio'] = b.zP / b.zF
    rng = np.random.default_rng(seed); detail = []
    for _, r in b.iterrows():
        sd_p = r.gradient_sd if pd.notna(r.gradient_sd) else 3.5
        zp = (r.zP * 5 + rng.normal(0, sd_p / np.sqrt(r.n), draws)) / 5
        zf = r.zF * (1 + rng.normal(0, 0.25 / np.sqrt(r.n), draws))     # flow per gram, reported spread ~25%
        q = np.percentile(zp / zf, [2.5, 97.5])
        detail.append(dict(study=r.study, group=r.group, outflow=r.outflow, zP=float(r.zP), zF=round(float(r.zF), 2),
                           ratio=round(float(r.ratio), 3), ratio_ci=[round(float(q[0]), 3), round(float(q[1]), 3)],
                           prob_below_one=float((zp / zf < 1).mean()), outcome_pct=float(r.pct)))
    rec = [x for x in detail if x['outflow'].startswith('reconstructed')]
    k = sum(x['ratio'] < 1 for x in rec)
    return dict(groups_with_both=int(len(b)), reconstructed=len(rec), standard=len(detail) - len(rec),
                ratio_reconstructed=[x['ratio'] for x in rec], all_below_one=bool(k == len(rec)),
                sign_test_p=float(binomtest(k, len(rec), 0.5, alternative='greater').pvalue) if rec else None,
                min_prob_below_one=float(min(x['prob_below_one'] for x in rec)) if rec else None,
                outcome_reconstructed_pct=[x['outcome_pct'] for x in rec],
                note='configuration classified from the methods before seeing outcomes; outcome shown separately',
                detail=detail)


def p3_within_cohort(d):
    """Primary analysis: hierarchical binomial model with the pressure index centred within each study, fitted by
    NUTS. Reported in the prespecified order (primary outcome, secondary outcome, all outcomes combined), with the
    sensitivity to the prior on tau, the overlap analysis for the Kyoto series, and the conditional
    internal-external validation, which scores the slope without fitting anything on the held-out centre."""
    import bayes
    c = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]
    def read(name, make):
        f = f'{OUT}/{name}.json'
        if os.path.exists(f): return json.load(open(f))
        r = make(); json.dump(r, open(f, 'w'), indent=1); return r
    main = read('bayes_main', lambda: {k: bayes.fit(v) for k, v in
                (('primary (SFSS or early dysfunction)', c[c.outcome_type == 'SFSS_or_dysfunction']),
                 ('secondary (mortality or graft loss)', c[c.outcome_type != 'SFSS_or_dysfunction']),
                 ('exploratory (all outcomes)', c))})
    iecv_ = read('bayes_iecv', lambda: bayes.conditional_iecv(c))
    spec = read('spec_curve', lambda: bayes.spec_curve(c))
    prio = read('bayes_priors', lambda: bayes.prior_sensitivity(c[c.outcome_type == 'SFSS_or_dysfunction']))
    meas = json.load(open(f'{OUT}/bayes_measurement.json')) if os.path.exists(f'{OUT}/bayes_measurement.json') else None
    recov = json.load(open(f'{OUT}/bayes_recovery.json')) if os.path.exists(f'{OUT}/bayes_recovery.json') else None
    pairs = []
    for st, g in c.groupby(indices.strata(c).values):
        g = g.sort_values('zP')
        if len(g) >= 2: pairs.append(dict(stratum=st, higher_zP_worse=bool(g.pct.iloc[-1] >= g.pct.iloc[0])))
    mus = [x['mu_beta'] for x in spec if x.get('mu_beta') is not None]
    return dict(primary=main['primary (SFSS or early dysfunction)'], secondary=main['secondary (mortality or graft loss)'],
                exploratory=main['exploratory (all outcomes)'],
                every_pair_in_predicted_direction=bool(all(p['higher_zP_worse'] for p in pairs)), pairs=pairs,
                prior_sensitivity=[dict(tau_prior=r['tau_prior'], mu_beta=r['mu_beta'], ci=r['mu_beta_ci'],
                                        tau_beta=r['tau_beta'], prob_positive=r['prob_mu_positive']) for r in prio],
                conditional_iecv=[dict(centre=v['held_out'], mu_from_others=v['mu_from_others'],
                                       log_score_gain=v['total_gain'], directions_correct=v['all_directions_correct']) for v in iecv_],
                specification_curve=spec, mu_range_across_specifications=[float(min(mus)), float(max(mus))],
                all_specifications_positive=bool(all(x['mu_beta'] > 0 for x in spec if x.get('mu_beta') is not None)),
                measurement_monte_carlo=meas, design_recovery=recov)

def p4_thresholds():
    """A correspondence of scales, not a statistical result: independently developed thresholds occupy a similar
    region once normalised. The exact landing at 2 depends on the 5 mmHg reference, so the same thresholds are
    also shown against 4 and 6 mmHg."""
    z = indices.zP
    sens = {f'{ref:g} mmHg reference': {k: round(v / ref, 2) for k, v in
            dict(cirrhosis_HVPG_10=10, graft_PVP15_CVP5=10, resection_HVPG_10=10, variceal_HVPG_12=12).items()}
            for ref in (4, 5, 6)}
    return dict(reference_sensitivity=sens, clinically_significant_portal_hypertension=dict(source='HVPG >= 10 mmHg (Baveno)', zP=z(gradient=10)),
                graft_threshold=dict(source='PVP 15 mmHg at CVP 5 (ILTS 2023)', zP=z(15, 5)),
                variceal_bleeding=dict(source='HVPG >= 12 mmHg', zP=z(gradient=12)),
                resection_limit=dict(source='HVPG 10 mmHg, any resection puts the remnant above', zP=z(gradient=10)),
                all_at_two=bool(abs(z(gradient=10) - 2) < 1e-9 and abs(z(15, 5) - 2) < 1e-9))


def main():
    os.makedirs(OUT, exist_ok=True)
    d = indices.compute(pd.read_csv(f'{ROOT}/data/series.tsv', sep='\t'))
    res = dict(P1_allometry=p1_allometry(), P2_discordance=p2_discordance(d),
               P3_within_cohort=p3_within_cohort(d), P4_thresholds=p4_thresholds())
    json.dump(res, open(f'{OUT}/predictions.json', 'w'), indent=1)
    a, b2, c3, d4 = res['P1_allometry'], res['P2_discordance'], res['P3_within_cohort'], res['P4_thresholds']
    print('P1  comparative anatomy')
    print(f"    over b = {a['b_range'][0]:.2f}-{a['b_range'][1]:.2f}: exponent of the shear set-point stays in "
          f"[{a['tau0_exponent_range'][0]:+.3f}, {a['tau0_exponent_range'][1]:+.3f}] -> invariant: {a['tau0_invariant']}")
    print(f"    vascular mass exponent spans {a['vascular_mass_exponent_range'][0]:.2f}-{a['vascular_mass_exponent_range'][1]:.2f}; "
          f"{100*a['fraction_compatible']:.0f}% of the range is compatible with the observed {a['observed_range']}")
    print('P2  flow-pressure discordance')
    print(f"    relative effective resistance in the {b2['reconstructed']} reconstructed-outflow groups: "
          f"{[round(x,2) for x in b2['ratio_reconstructed']]}, all below 1: {b2['all_below_one']} "
          f"(sign test p = {b2['sign_test_p']:.3f}; smallest posterior probability of being below 1: {b2['min_prob_below_one']:.2f})")
    print(f"    their outcomes, reported separately: {b2['outcome_reconstructed_pct']}%")
    print('P3  within-cohort ordering (primary analysis: hierarchical, zP centred within study)')
    for k in ('primary', 'secondary', 'exploratory'):
        r = c3[k]
        print(f"    {k:12s} {r['groups']:2d} groups / {len(r['strata'])} strata: mu {r['mu_beta']:+.2f} "
              f"({r['mu_beta_ci'][0]:+.2f} to {r['mu_beta_ci'][1]:+.2f}), OR {r['OR_per_zP']:.2f}, "
              f"P(mu>0) = {r['prob_mu_positive']:.3f}, tau {r['tau_beta']:.2f}; R-hat {r['rhat']['max']:.3f}, "
              f"ESS {r['ess']['min']:.0f}, divergences {r['divergences']}, PPC {r['ppc_inside']}/{r['ppc_total']}")
    print(f"    every pair in the predicted direction: {c3['every_pair_in_predicted_direction']}")
    print('    prior on tau: ' + '; '.join(f"HN({p['tau_prior']}) mu {p['mu_beta']:+.2f} ({p['ci'][0]:+.2f},{p['ci'][1]:+.2f}) P>0 {p['prob_positive']:.3f}" for p in c3['prior_sensitivity']))
    print('    conditional validation (slope only, no intercept fitted on the held-out centre):')
    for v in c3['conditional_iecv']:
        print(f"      {v['centre']:10s} mu from the others {v['mu_from_others']:+.2f}, log-score gain over no association {v['log_score_gain']:+.2f}, directions correct {v['directions_correct']}")
    print(f"    specification curve: mu between {c3['mu_range_across_specifications'][0]:+.2f} and "
          f"{c3['mu_range_across_specifications'][1]:+.2f} across {len(c3['specification_curve'])} specifications; all positive: {c3['all_specifications_positive']}")
    if c3['measurement_monte_carlo']:
        mm = c3['measurement_monte_carlo']
        print(f"    with measurement uncertainty propagated: mu {mm['mu_beta']:+.2f} ({mm['mu_beta_ci'][0]:+.2f} to {mm['mu_beta_ci'][1]:+.2f}), P>0 {mm['prob_positive']:.2f}")
    if c3['design_recovery']:
        print('    design recovery: ' + '; '.join(f"true {r['mu_true']:.1f} -> {r['mean_estimate']:.2f} (coverage {r['coverage']:.2f}, power {r['power_to_exclude_zero']:.2f})" for r in c3['design_recovery']))
    print('P4  thresholds of three fields on one scale')
    for k, v in d4.items():
        if isinstance(v, dict) and 'source' in v: print(f"    {v['source']:45s} zP = {v['zP']:.2f}")
    print('    with a reference gradient of 4, 5 or 6 mmHg the same thresholds give:')
    for k, v in d4['reference_sensitivity'].items(): print(f"      {k:18s} {v}")


if __name__ == '__main__':
    main()
