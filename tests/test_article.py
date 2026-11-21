"""Every number the article quotes, checked against the results on disk (run_all.sh writes them).

One line per claim: where it is quoted, the value printed in the article, the value recomputed here and, for a number
that comes out of a sampler, how far it may move. If a result changes, this test names the sentence of the article
that has to change with it.

Two kinds of number:
  deterministic   counts, arithmetic on the table, the network model: compared exactly, at the digits the article prints;
  sampled         posterior means and interval limits (NUTS, Metropolis). At a fixed seed they are reproducible on one
                  version of the numerical libraries (requirements-lock.txt) and move in the second decimal on another:
                  jax 0.12 gives the primary slope +1.56 (+0.37 to +2.76) where the locked 0.11 gives +1.55 (+0.33 to
                  +2.79). The tolerances below are a few Monte Carlo standard errors: SLOPE for a posterior mean, LIMIT
                  for a limit of a 95% interval, WIDE for a limit of the interval of a new centre and for a log score.
LIPNET_EXACT=1 asks for the printed digits also of the sampled numbers; that is the test on the locked environment.
"""
import json, os, sys
import pandas as pd
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
import indices

SLOPE, LIMIT, WIDE = 0.05, 0.10, 0.30
EXACT = bool(os.environ.get('LIPNET_EXACT'))

# Results that are not kept in the repository and that this test needs. On a clone where run_all.sh has not been run
# (the continuous integration, a first look) each is written by the script that owns it, so that the test checks the
# article against the code and never fails for a missing file.
BUILDERS = {'model_human': ['src/model/human.py'], 'model_distension': ['src/model/distension.py'], 'ipd_ecological': ['src/ipd/run_ipd.py', 'ecological', '--quiet']}


def R(name):
    f = os.path.join(ROOT, 'results', name + '.json')
    if not os.path.exists(f) and name in BUILDERS:
        import subprocess
        os.makedirs(os.path.dirname(f), exist_ok=True)
        r = subprocess.run([sys.executable] + BUILDERS[name], cwd=ROOT, env=dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg'), capture_output=True, text=True)
        assert r.returncode == 0 and os.path.exists(f), f"{' '.join(BUILDERS[name])} did not write results/{name}.json:\n{r.stderr[-1500:]}"
    return json.load(open(f))


r2 = lambda x: round(float(x), 2)


def claims():
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t')); p = d.dropna(subset=['zF', 'zP'])
    P = R('predictions'); p1, p2, p3 = P['P1_allometry'], P['P2_discordance'], P['P3_within_cohort']; m, sw = p2['model_against_observed'], p2['reference_sweep']
    B = R('bayes_main'); pri, sec, exp = (B[k] for k in ('primary (SFSS or early dysfunction)', 'secondary (mortality or graft loss)', 'exploratory (all outcomes)'))
    prior = {x['tau_prior']: x for x in R('bayes_priors')}; H = R('model_human'); E = R('ipd_ecological'); S = R('spec_curve'); iecv = R('bayes_iecv'); w = p3['same_cohort_contrast'][0]
    wu = H['two_inflows_measured']['Wu 2011']; cast = H['cast_exponent']['portal vein']['all']; q = [x for x in H['portal_radii'] if abs(x['b'] - 0.75) < 1e-9][0]['model_over_cast']
    sweep = {s['sigma_mmHg']: s['mu_beta'] for s in E['sweep']}; X = R('model_distension'); bk = X['distension']['bracket']; pc = X['distension']['published_check']['states']; pf = p2['performance']; rb = X['robustness']['scenarios']; J = X['joint']
    return [
        ('Methods: groups, series', (32, 18), (len(d), d.study.nunique())),
        ('Results: resistance loads of the six groups', [0.31, 0.54, 0.55, 0.65, 0.72, 0.74], sorted(r2(x) for x in p.zR)),
        ('Results: flow loads 3.3 to 4.1', (3.3, 4.1), (round(p.zF.min(), 1), round(p.zF.max(), 1))),
        ('Results: cohorts below one, of', (5, 5), (p2['cohorts_below_one'], p2['cohorts'])),
        ('Results, Fig. 2B: exact sign test, one-sided', 0.031, round(p2['sign_test_cohorts_one_sided_p'], 3)),
        ('Results: network prediction for the groups, 0.81 to 0.90', (0.81, 0.90), (r2(min(g['zR_predicted'][0] for g in m['groups'])), r2(max(g['zR_predicted'][1] for g in m['groups'])))),
        ('Results: the cast network with the donor arterial flow gives 0.88 to 0.89 against an observed 0.58; root-mean-square error 0.34, 0.27 for the best single configuration', (0.88, 0.89, 0.58, 0.34, 0.27), (r2(min(pf['central_configuration']['zR'])), r2(max(pf['central_configuration']['zR'])), r2(m['mean_observed']), r2(pf['central']['rmse']), r2(pf['best_admissible']['rmse']))),
        ('Results: 16 to 73% of the departure', [0.16, 0.73], m['share_of_the_departure_predicted']),
        ('Results: from the flow loads alone every model predicted a ratio between the arms of 0.95 to 1.00, against 0.74 (0.61 to 0.88)', (0.95, 1.00, 0.74, 0.61, 0.88), (r2(min(v['ratio'][1] for v in X['wang_contrast']['models'].values())), r2(max(v['ratio'][1] for v in X['wang_contrast']['models'].values())), r2(X['wang_contrast']['observed_ratio'][1]), r2(X['wang_contrast']['observed_ratio'][0]), r2(X['wang_contrast']['observed_ratio'][2]))),
        ('Results: secondary outcome +1.79 with deaths counted directly', 1.79, r2(R('bayes_extra')['censored_survival']['secondary, Ogura 2010 counted directly']['mu_beta']), SLOPE),
        ('Abstract: published normal-liver relations give 0.39 to 0.79', (0.39, 0.79), (r2(min(g['zR_distensible_all'][0] for g in bk['groups'])), r2(max(g['zR_distensible_all'][1] for g in bk['groups'])))),
        ('Results: below one in 54 of 60 combinations', (54, 60), (sw['joint']['all_below_one'], sw['joint']['combinations'])),
        ('Results: below the network in 34 of 60, none at 4 mmHg, all at 6 mmHg', (34, 0, 20), (sw['joint']['all_below_model'], sw['joint']['by_normal_gradient']['4 mmHg']['all_below_model'], sw['joint']['by_normal_gradient']['6 mmHg']['all_below_model'])),
        ('Results: no manoeuvre 0.72 inside outflow reconstruction 0.31 to 0.74', ([0.72], 0.31, 0.74), ([r2(x) for x in p2['by_configuration']['none']], r2(min(p2['by_configuration']['outflow reconstruction'])), r2(max(p2['by_configuration']['outflow reconstruction'])))),
        ('Results: splenectomy 0.54 against 0.72 in its cohort', [0.54, 0.72], [r2(x) for x in p2['within_cohort_contrasts'][0]['zR']]),
        ('Results: sinusoids carry 3 to 15% of the drop', (0.03, 0.15), (r2(min(X['sinusoids']['share_of_the_drop'].values())), r2(max(X['sinusoids']['share_of_the_drop'].values())))),
        ('Results: relations that keep distending 0.39 to 0.58, the one that stops 0.61 to 0.79', (0.39, 0.58, 0.61, 0.79), (r2(min(g['zR_keeps_distending'][0] for g in bk['groups'])), r2(max(g['zR_keeps_distending'][1] for g in bk['groups'])), r2(min(g['zR_stops_distending'][0] for g in bk['groups'])), r2(max(g['zR_stops_distending'][1] for g in bk['groups'])))),
        ('Results: five groups within the relations, Chan 2011 below', (5, ['Chan 2011 / right lobe with MHV, no modulation']), (bk['within_the_relations'], bk['below_the_relations'])),
        ('Supplement: doubling the flow of a donor liver raises the gradient by 1.7 to 2.3 mmHg with the first two relations, 3.3 to 3.5 with the third', (1.7, 2.3, 3.3, 3.5), (round(min(X['distension']['flow_doubled_from_the_donor'][l]['rise_mmHg'][0] for l in ('cube law', 'compliance')), 1), round(max(X['distension']['flow_doubled_from_the_donor'][l]['rise_mmHg'][1] for l in ('cube law', 'compliance')), 1), round(X['distension']['flow_doubled_from_the_donor']['pressure-flow']['rise_mmHg'][0], 1), round(X['distension']['flow_doubled_from_the_donor']['pressure-flow']['rise_mmHg'][1], 1))),
        ('Results: pig livers at 2.7 times the flow, resistance 0.54 and 0.59 of the starting value; the compliance relation 0.53 to 0.63; the relation that stops distending 0.78 to 0.81', (0.54, 0.59, 0.53, 0.63, 0.78, 0.81), (r2(pc['in situ']['measured_resistance_ratio']), r2(pc['isolated, perfused']['measured_resistance_ratio']), r2(min(pc[k]['predicted']['compliance']['resistance_ratio'][0] for k in ('in situ', 'isolated, perfused'))), r2(max(pc[k]['predicted']['compliance']['resistance_ratio'][1] for k in ('in situ', 'isolated, perfused'))), r2(min(pc[k]['predicted']['pressure-flow']['resistance_ratio'][0] for k in ('in situ', 'isolated, perfused'))), r2(max(pc[k]['predicted']['pressure-flow']['resistance_ratio'][1] for k in ('in situ', 'isolated, perfused'))))),
        ('Results: the further 28 to 35% by which resistance would have to fall; supplement, 8 to 11% in radius', (0.28, 0.35, 1.08, 1.11), (r2(X['required']['reduction_missing'][0]), r2(X['required']['reduction_missing'][1]), r2(X['required']['whole_bed_radius_factor'][0]), r2(X['required']['whole_bed_radius_factor'][1]))),
        ('Results: vessel by vessel, below 0.65 only if the portal tree distends; the hepatic venous tree alone 0.67 to 0.92', (True, 0.67, 0.92), (all(v['keeps_distending_span'][0] < 0.65 for k, v in rb.items() if k != 'site: postsinusoidal site') and rb['site: postsinusoidal site']['keeps_distending_span'][0] > 0.65, r2(rb['site: postsinusoidal site']['keeps_distending_span'][0]), r2(rb['site: postsinusoidal site']['keeps_distending_span'][1]))),
        ('Results: pig livers, the second relation within 1.2 mmHg, the third 2.6 to 5.0 mmHg too high', (1.2, 2.6, 5.0), (round(max(abs(x) for k in ('in situ', 'isolated, perfused') for x in pc[k]['predicted']['compliance']['portal_pressure_error_mmHg']) + 0.049, 1), round(min(x for k in ('in situ', 'isolated, perfused') for x in pc[k]['predicted']['pressure-flow']['portal_pressure_error_mmHg']), 1), round(max(x for k in ('in situ', 'isolated, perfused') for x in pc[k]['predicted']['pressure-flow']['portal_pressure_error_mmHg']), 1))),
        ('Results: every uncertainty drawn together, Chan 2011 below all three in >99% of draws, the others within in 31 to 66%', (True, 0.31, 0.66), ([g for g in J['groups'] if g['study'] == 'Chan 2011'][0]['prob_observed_below_every_relation'] > 0.99, r2(min(g['prob_observed_within_the_relations'] for g in J['groups'] if g['study'] != 'Chan 2011')), r2(max(g['prob_observed_within_the_relations'] for g in J['groups'] if g['study'] != 'Chan 2011')))),
        ('Results: Wu 2011, arterial flow 0.66, portal 3.25-fold, model 0.82 to 0.86, observed 0.65', (0.66, 3.25, [0.82, 0.86], 0.65), (r2(wu['arterial_load']), r2(wu['flow_load']), [r2(x) for x in wu['zR_model']], r2(wu['zR_observed']))),
        ('Results: exponent of tau0 between -0.077 and +0.071', (-0.077, 0.071), tuple(round(x, 3) for x in p1['tau0_exponent_range'])),
        ('Results: blood volume scaling met in 5 of 28', 5, round(p1['fraction_compatible'] * 28)),
        ('Results: cast radii within 10% at b = 0.75', True, bool(0.895 <= min(q) and max(q) < 1.095)),
        ('Results: radius law of the cast, b = 0.73 (0.71 to 0.74)', (0.73, [0.71, 0.74]), (r2(cast['b']), [r2(x) for x in cast['leave_one_out']])),
        ('Results: primary slope +1.55 (+0.33 to +2.79)', (1.55, [0.33, 2.79]), (r2(pri['mu_beta']), [r2(x) for x in pri['mu_beta_ci']]), (SLOPE, LIMIT)),
        ('Abstract, Results: odds ratio 4.70 per unit of zP, 1.36 per mmHg', (4.70, 1.36), (r2(pri['OR_per_zP']), r2(pri['OR_per_mmHg'])), (4.70 * SLOPE, 0.02)),
        ('Results: secondary +1.74, all outcomes +1.79 (+0.97 to +2.64)', (1.74, 1.79, [0.97, 2.64]), (r2(sec['mu_beta']), r2(exp['mu_beta']), [r2(x) for x in exp['mu_beta_ci']]), (SLOPE, SLOPE, LIMIT)),
        ('Results: 15 specifications between +0.92 and +1.95', (15, 0.92, 1.95), (len(S), r2(min(s['mu_beta'] for s in S)), r2(max(s['mu_beta'] for s in S))), (0, LIMIT, LIMIT)),
        ('Results: widest prior, +1.46 (-0.07 to +2.84)', (1.46, [-0.07, 2.84]), (r2(prior[1.0]['mu_beta']), [r2(x) for x in prior[1.0]['mu_beta_ci']]), (SLOPE, LIMIT)),
        ('Results: slope in a new centre, -0.02 to +3.09', [-0.02, 3.09], [r2(x) for x in pri['prediction_interval']], WIDE),
        ('Results: same cohort, 19.7% against 9.7%, worse arm with higher pressure and lower flow', (19.7, 9.7, True, False), (round(w['outcome_pct'][0], 1), round(w['outcome_pct'][1], 1), w['worse_group_has_higher_pressure_load'], w['worse_group_has_higher_flow_load'])),
        ('Results: log-score gains +2.26, +2.74, +8.40', [2.26, 2.74, 8.40], sorted(r2(v['total_gain']) for v in iecv), WIDE),
        ('Results: slopes +1.58 to +1.74 for spreads of 1 to 4 mmHg, +1.56 at zero', (1.56, 1.58, 1.74), (r2(sweep[0.0]), r2(min(sweep[k] for k in (1.0, 2.0, 3.0, 4.0))), r2(max(sweep[k] for k in (1.0, 2.0, 3.0, 4.0)))), SLOPE),
        ('Discussion: odds about 0.21 times as high for a gradient 5 mmHg lower', 0.21, r2(1 / pri['OR_per_zP']), 0.02),
        ('Results: thresholds at zP = 2 and 2.4', (2.0, 2.0, 2.4), (P['P4_thresholds']['clinically_significant_portal_hypertension']['zP'], P['P4_thresholds']['graft_threshold']['zP'], P['P4_thresholds']['variceal_bleeding']['zP'])),
    ]


def flat(x):
    return [y for v in x for y in flat(v)] if isinstance(x, (list, tuple)) else [x]


def spread(tol, said):
    """One tolerance per printed number: a scalar for all of them, or the same nesting as the claim."""
    if not isinstance(tol, (list, tuple)): return [tol] * len(flat(said))
    return [t for part, tl in zip(said, tol) for t in spread(tl, part)]


def agrees(said, got, tol=0):
    a, b = flat(said), flat(got)
    if len(a) != len(b): return False
    return all(x == y if (EXACT or not t or isinstance(x, bool)) else abs(x - y) <= t + 1e-12 for x, y, t in zip(a, b, spread(tol, said)))


def test_every_number_of_the_article_is_the_one_on_disk():
    wrong = [f'{c[0]}: the article says {c[1]}, the results give {c[2]}' + (f' (tolerance {c[3]})' if len(c) > 3 and not EXACT else '') for c in claims() if not agrees(*c[1:])]
    assert not wrong, '\n'.join(wrong)


def test_the_tolerance_is_only_for_sampled_numbers():
    c = claims(); sampled = [x for x in c if len(x) > 3]
    assert 8 <= len(sampled) < len(c) / 2, 'most claims are deterministic and compared exactly'
    assert all(max(flat(x[3])) <= max(WIDE, 4.70 * SLOPE) for x in sampled)
    assert agrees((1.55, [0.33, 2.79]), (1.56, [0.37, 2.76]), (SLOPE, LIMIT)) or EXACT, 'the same model on another version of jax'
    assert not agrees((1.55, [0.33, 2.79]), (1.75, [0.33, 2.79]), (SLOPE, LIMIT)), 'and a different result is still caught'
    assert not agrees((5, 5), (4, 5)), 'a count is compared exactly'


if __name__ == '__main__':
    for c in claims(): print('ok  ' if agrees(*c[1:]) else 'FAIL', 'sampled      ' if len(c) > 3 else 'deterministic', c[0], '|', c[2])
