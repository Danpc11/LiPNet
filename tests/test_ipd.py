"""Checks for the patient-level layer: the Kaplan-Meier reconstruction has to recover a
cohort whose truth we know, and the integrated model has to reproduce the published fit
when the within-group spread is set to zero.

    python -m pytest tests/test_ipd.py -q      or      python tests/test_ipd.py
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))

import indices
from ipd import km_reconstruct as km
from ipd import ecological as ec

# the published NUTS fit of the primary outcome, results/bayes_main.json
NUTS_MU, NUTS_GROUPS, NUTS_PATIENTS, NUTS_EVENTS = 1.547, 7, 573, 76


def cohort(n=200, rate=0.10, seed=1, horizon=20.0):
    rng = np.random.default_rng(seed)
    T = rng.exponential(1.0 / rate, n)
    C = rng.uniform(0.0, horizon, n)
    return np.minimum(T, C), (T <= C).astype(int)


def primary():
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    if 'in_main_analysis' in d:
        d = d[d.in_main_analysis == True]
    return d[d.outcome_type == 'SFSS_or_dysfunction']


# ------------------------------------------------------------------ Kaplan-Meier
def test_kaplan_meier_is_the_textbook_estimator():
    t = np.array([1., 2., 3., 4., 5.])
    s = np.array([1, 0, 1, 1, 0])
    ts, ss, nr = km.kaplan_meier(t, s)
    assert list(ts) == [1., 3., 4.]
    assert np.allclose(ss, [4 / 5, (4 / 5) * (2 / 3), (4 / 5) * (2 / 3) * (1 / 2)])
    assert list(nr) == [5, 3, 2]


def test_reconstruction_recovers_a_cohort_it_was_not_given():
    """The whole point: from the step function and the numbers at risk, and nothing else,
    get back the patients. Cohort size is fixed exactly by the risk table; the event count
    and the curve come back to within the digitisation error."""
    for n, rate, seed in ((200, 0.10, 1), (120, 0.25, 2), (400, 0.05, 3)):
        t, s = cohort(n, rate, seed)
        pts = km.digitise(t, s, n_points=40, t_max=18.0)
        rk = km.risk_table(t, n_times=6, t_max=18.0)
        r = km.reconstruct(pts, rk, total_events=int(s.sum()))
        assert r['n_patients'] == n, 'cohort size must come back exactly'
        assert abs(r['n_events'] - int(s.sum())) <= max(3, 0.05 * s.sum()), \
            'events %d vs %d' % (r['n_events'], s.sum())
        assert r['max_abs_error'] < 0.03, 'curve off by %.3f' % r['max_abs_error']
        assert r['rmse'] < 0.01


def test_reconstruction_needs_no_event_total():
    """Papers do not always print the number of events. Without it the curve must still be
    reproduced, since the steps themselves carry the events."""
    t, s = cohort(200, 0.10, 1)
    pts = km.digitise(t, s, n_points=40, t_max=18.0)
    rk = km.risk_table(t, n_times=6, t_max=18.0)
    r = km.reconstruct(pts, rk)
    assert r['max_abs_error'] < 0.05
    assert 0.5 * s.sum() < r['n_events'] < 1.5 * s.sum()


def test_more_digitised_points_give_a_closer_curve():
    t, s = cohort(300, 0.12, 7)
    rk = km.risk_table(t, n_times=6, t_max=18.0)
    errs = [km.reconstruct(km.digitise(t, s, k, t_max=18.0), rk,
                           total_events=int(s.sum()))['rmse'] for k in (12, 30, 60)]
    assert errs[2] <= errs[0], 'reading more points off the figure must not hurt: %s' % errs


# ------------------------------------------------------------------ patient-level data
def test_the_transcription_reproduces_the_published_summaries():
    """The control on the transcription: every statistic Botha et al. report in their own
    tables and text is recomputed from the transcribed rows. One disagreement is expected
    and is in the paper, not in the transcription: Table 4 gives the minimum pre-shunt
    gradient as 12 mmHg while Table 5 lists a patient at 10."""
    from ipd import patient_level as pl
    v = pl.verify()
    bad = v[~v.agrees]
    assert len(v) >= 12
    assert set(bad.quantity) <= {'minimum gradient before'},         'unexpected disagreement: %s' % bad.to_dict('records')


def test_the_shunt_lowers_every_recipients_pressure_load():
    from ipd import patient_level as pl
    e = pl.shunt_effect()
    assert e['n_paired'] == 10
    assert e['all_fell'], 'the shunt should lower the gradient in every paired recipient'
    assert e['zP_before_median'] > 3.0 > e['zP_after_median']


def test_the_one_failure_is_the_highest_pressure_load():
    from ipd import patient_level as pl
    for outcome in ('sfss', 'graft_loss_or_death'):
        r = pl.event_rank_test(outcome=outcome)
        assert r['events'] == 1
        assert r['events_are_the_top_k'], '%s did not fall on the highest gradient' % outcome
        assert abs(r['p_events_in_top_k'] - 0.1) < 1e-9
        assert r['event_gradients'] == [15.0]


def test_ogura_dispersions_are_standard_errors():
    """Ogura's Results section gives the MELD of the same cohort as 23.2 +/- 1.1 (n = 86)
    and 22.9 +/- 1.7 (n = 43), while its Table 1 gives 22.9 +/- 10.6 for that cohort.
    Multiplying by sqrt(n) recovers the Table 1 value, which identifies the smaller
    figures as standard errors and licenses the same conversion for the gradient."""
    assert abs(ec.sem_to_sd(1.1, 86) - 10.6) < 0.6
    assert abs(ec.sem_to_sd(1.7, 43) - 10.6) < 0.6
    sds = ec.corrected_sigma_mmHg()
    assert min(sds) > 2.0, 'after correction no series should look artificially tight: %s' % sds
    assert max(sds) < 6.0
    assert 0.4 < ec.pooled_sigma(None if False else primary()) < 1.0


def test_chan_gives_an_outcome_contrast_with_no_pressure_separation():
    """Chan et al. report flow and pressure separately for the recipients who had a
    Clavien grade 2+ complication and those who did not: 12 events in 46 recipients,
    against the 2 in-hospital deaths the main extraction scores. The two groups have the
    same median portal pressure, so the series cannot enter a slope model, and the
    extremes of pressure sat among the uncomplicated recipients. It is a null, and it
    belongs in the record as one."""
    from ipd import patient_level as pl
    c = pl.chan_outcome_contrast()
    assert c['n_total'] == 46 and c['events_total'] == 12
    assert c['pressure_separation_mmHg'] == 0.0
    assert not c['usable_in_the_slope_model']
    assert c['extremes_were_uncomplicated']
    assert abs(c['zF_uncomplicated'] - c['zF_complicated']) < 0.2
    assert 0.25 < c['zR_uncomplicated'] < 0.35


def test_a_range_recovers_a_plausible_dispersion():
    from ipd import patient_level as pl
    assert abs(pl.sd_from_range(0.0, 4.0 * 1.0, 46) - 4.0 / 4.42) < 0.05
    assert 2.0 < pl.sd_from_range(-1.0, 11.0, 46) < 3.5


def test_ogura_table_three_is_a_patient_level_outcome_record():
    """Every death of that series is listed individually with its final portal pressure.
    Sorting them into the two pressure bands must reproduce the counts the paper states in
    its text, which is what would license digitising its Fig. 4."""
    from ipd import patient_level as pl
    o = pl.ogura_deaths()
    assert all(c['agrees'] for c in o['checks'].values()), o['checks']
    assert o['pressures'] == [8, 12, 13, 14, 15, 15, 16, 17, 17, 17, 17, 17, 18, 18, 18, 26]
    assert o['rate_at_or_above_15'] > 5 * o['rate_below_15']


# ------------------------------------------------------------------ ecological bias
def test_the_sampler_reproduces_the_published_fit_at_zero_spread():
    """With sigma = 0 the integrated model collapses to the group-mean model of
    src/bayes.py. This is what licenses the rest: the NumPy sampler is checked against the
    published NUTS fit on the same data before it is used to say anything new."""
    r = ec.fit(primary(), sigma=0.0, chains=4, warmup=6000, samples=15000)
    assert r['groups'] == NUTS_GROUPS
    assert r['patients'] == NUTS_PATIENTS
    assert r['events'] == NUTS_EVENTS
    assert abs(r['mu_beta'] - NUTS_MU) < 0.2, 'mu %.3f vs %.3f' % (r['mu_beta'], NUTS_MU)
    assert r['mu_beta_ci'][0] > 0, 'the published interval excludes zero'
    assert r['rhat'] < 1.05, 'R-hat %.3f' % r['rhat']


def test_every_series_agrees_on_the_size_of_the_within_group_spread():
    """Once Ogura's standard errors are converted and Botha's table is used, all seven
    dispersions fall between 2.6 and 4.4 mmHg, which is what licenses one pooled value."""
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    assert d.gradient_sd.notna().sum() >= 5
    sds = ec.corrected_sigma_mmHg()
    assert len(sds) >= 7, 'every series with a dispersion should be represented: %s' % sds
    assert max(sds) / min(sds) < 2.0, \
        'after correction the dispersions should agree in magnitude: %s' % sds
    assert 0.0 < ec.pooled_sigma(primary()) < 2.0


def test_within_group_spread_does_not_overturn_the_association():
    """The slope stays positive across every assumed spread, and over the range the series
    actually report it is no smaller than the published estimate."""
    sw = ec.sweep(primary(), (0.0, 0.4, 0.8), chains=4, warmup=6000, samples=15000)
    mu = np.array([s['mu_beta'] for s in sw])
    assert np.all(mu > 0), 'slope changed sign: %s' % mu
    assert np.all(np.array([s['prob_positive'] for s in sw]) > 0.9)
    assert mu[1] >= mu[0] - 0.1, 'ecological bias must not inflate the published estimate'
    assert np.all(np.array([s['rhat'] for s in sw]) < 1.05)


if __name__ == '__main__':
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                fn()
                print('ok   ', name)
            except AssertionError as e:
                fails += 1
                print('FAIL ', name, '--', e)
    raise SystemExit(fails)
