"""Blood viscosity (src/model/viscosity.py) and the published hepatectomies (src/model/hepatectomy.py)."""
import json, os, sys
import numpy as np, pandas as pd
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from model import viscosity as vs, hepatectomy as hp, distension as ds, compartments as cm
import indices


def _d():
    return indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))


def test_the_viscosity_law_and_its_limits():
    assert abs(float(vs.relative_viscosity(1e4, 0.45)) - 3.2) < 2e-3 and abs(float(vs.relative_viscosity(1e4, 0.0)) - 1.0) < 2e-3     # large vessel: 3.2 at 45%, plasma at 0
    assert float(vs.relative_viscosity(8.0, 0.45)) > float(vs.relative_viscosity(40.0, 0.45)) < float(vs.relative_viscosity(1e3, 0.45))  # the minimum between capillaries and large vessels
    for b in vs.B_HUMAN:
        assert abs(vs.factor(0.42, 0.42, b) - 1) < 1e-12
        f = [vs.factor(H, 0.42, b) for H in (0.20, 0.25, 0.30, 0.35, 0.42, 0.50)]; assert all(x < y for x, y in zip(f, f[1:]))        # resistance rises with the haematocrit
    lo, hi = vs.factor_range(); assert 0.55 < lo < 0.65 and 0.85 < hi < 0.92
    s = vs.factor(0.30, 0.42, segments=True); assert abs(s['total'] - vs.factor(0.30, 0.42)) < 1e-12 and min(s.values()) > 0.6 and max(s.values()) < 0.85
    x = vs.draw(np.random.default_rng(0), 4000); assert lo - 1e-9 <= x.min() and x.max() <= hi + 1e-9


def test_the_law_holds_in_perfused_livers_it_was_not_derived_from():
    c = vs.published_check(); t = c['text_check']
    assert all(abs(t['given_in_the_text'][k] - t['digitised'][k]) < 0.003 for k in t['given_in_the_text'])          # the digitised figure gives back the values of the text
    k = pd.read_csv(os.path.join(ROOT, 'data', 'kamikado2004.tsv'), sep='\t')
    assert len(k) == 20 and np.allclose(k.R_portal + k.R_sinusoidal + k.R_hepatic_venous, k.R_total, atol=2e-3) and set(k.basis) == {'digitised'}
    for name, st in c['states'].items():
        o = st['observed_mean']['total']; p = st['predicted']['network at rat size']['total']
        assert 0.80 < o / p < 1.10, name                                                                           # within 20% at every haematocrit, and the law never overstates the fall
    z = c['states']['haematocrit 0 over 0.30']; assert abs(z['observed_mean']['total'] - z['predicted']['network at rat size']['total']) < 0.03
    assert z['observed_mean']['sinusoidal'] < z['observed_mean']['portal'] and z['predicted']['human network']['sinusoidal'] < z['predicted']['human network']['portal']   # sinusoids change most, in the livers and in the law


def test_viscosity_enters_the_comparison_without_changing_the_one_without_it():
    d = _d(); a = ds.bracket(d); b = ds.bracket(d, viscosity=(1.0, 1.0)); assert a == b
    v = vs.factor_range(); c = ds.bracket(d, viscosity=v)
    for x, y in zip(a['groups'], c['groups']):
        assert abs(y['zR_rigid'][0] - x['zR_rigid'][0] * v[0]) < 2e-3 and abs(y['zR_rigid'][1] - x['zR_rigid'][1] * v[1]) < 2e-3
        assert y['zR_keeps_distending'][0] < x['zR_keeps_distending'][0] and y['zR_observed'] == x['zR_observed']
    assert a['counts']['rigid network']['point_inside'] == 0 and c['counts']['rigid network']['point_inside'] == 5      # five of six inside the rigid network with the blood of the recipient
    assert c['counts']['all relations']['interval_reaches'] == 6
    j0 = ds.joint(d, draws=3000); j1 = ds.joint(d, draws=3000, haematocrit=((0.42, 0.42), (0.42, 0.42))); j2 = ds.joint(d, draws=3000, haematocrit=(vs.HCT_RECIPIENT, vs.HCT_REFERENCE))
    for x, y, z in zip(j0['groups'], j1['groups'], j2['groups']):
        assert x['zR_observed'] == y['zR_observed'] == z['zR_observed']                                             # the other draws are the same draws
        assert np.allclose(x['models']['rigid network']['zR'], y['models']['rigid network']['zR'], atol=1e-3)       # same haematocrit: factor one
        assert z['models']['rigid network']['zR'][1] < x['models']['rigid network']['zR'][1]
        assert z['models']['rigid network']['prob_observed_below'] <= x['models']['rigid network']['prob_observed_below']


def test_the_two_clinical_thresholds_are_one_point_of_the_plane():
    t = vs.thresholds(); assert abs(t['zF'] - 250 / 90) < 1e-3 and t['zP'] == 2.0 and abs(t['zR_where_they_agree'] - 0.72) < 1e-9 and t['flow_at_which_donor_resistance_reaches_the_gradient'] == 180.0
    g = t['gradient_at_the_flow_threshold_mmHg']; assert g['normal_blood'][0] > t['gradient_mmHg'] > g['haematocrit 0.25'][1]          # above the pressure threshold with normal blood, below it at 25%
    assert g['haematocrit 0.25'][0] < g['haematocrit 0.3'][0] < g['haematocrit 0.35'][0] < g['normal_blood'][0]
    r = t['rise_of_the_gradient_from_0p25_to_0p35']; assert 1.2 < r[0] <= r[1] < 1.35
    k = t['curve']; assert all(a <= b for a, b in zip(k['gradient_low'], k['gradient_high'])) and all(x < y for x, y in zip(k['gradient_low'], k['gradient_low'][1:]))
    d = _d(); s = vs.clinical(d, draws=2000)['span']; assert s['blood_of_the_recipient']['rigid'][0] < s['normal_blood']['rigid'][0] and s['blood_of_the_recipient']['keeps_distending'][0] < s['blood_of_the_recipient']['stops_distending'][0]


def test_the_hepatectomy_table_has_a_source_for_every_state():
    d = hp.load(); assert len(d) == 10 and d.notes.notna().all() and set(d.dispersion) <= {'SD', 'SEM'} and set(d.flow_basis) <= {'per 100 g', 'absolute'}
    assert set(d.out_basis) <= {'hepatic vein pressure', 'central venous pressure', 'not reported'} and set(d.shunt) <= {'none', 'side-to-side'}
    assert (d.flow_before > 0).all() and (d.PVP_after > 0).all() and ((d.resected_fraction >= 0) & (d.resected_fraction < 1)).all()
    r = d[d.gradient_before.notna()]; assert np.allclose(r.PVP_before - r.out_before, r.gradient_before, atol=1.0) and np.allclose(r.PVP_after - r.out_after, r.gradient_after, atol=0.2)   # reported gradient against portal minus caval pressure (0.9 mmHg apart before the 90% resection)
    assert d[d.out_basis == 'not reported'].out_before.isna().all() and d[d.out_basis != 'not reported'].out_before.notna().all()


def test_a_remnant_in_the_rigid_network_and_the_state_before_resection():
    G = cm.geometries(); a = cm.A_RATIO
    assert np.allclose(hp.rigid(1.0, 1.0, 1.0, a, G), 1.0)                                                          # nothing resected, same flows: the liver before
    lo, hi = hp.rigid(3.0, 0.25, 1.0, a, G); assert 0.75 < lo < hi < 1.0
    assert hp.rigid(3.0, 0.25, 1.0, a, G)[0] < hp.rigid(3.0, 0.75, 1.0, a, G)[0]                                    # the trunk carries less in a smaller remnant
    assert hp.rigid(3.0, 0.25, 0.0, a, G)[1] < hp.rigid(3.0, 0.25, 3.0, a, G)[0]                                    # arterial flow raises the load behind the sinusoid
    p = hp.predicted(1.0, 1.0, 1.0, a, 5.0, 4.0, 1.0, G); assert all(abs(x - 1) < 1e-6 for v in p.values() for x in v)
    q = hp.predicted(4.0, 0.25, 1.0, a, 5.0, 4.0, 1.0, G); assert q['cube law'][1] < q['pressure-flow'][0] < q['rigid network'][0]


def test_resistance_falls_with_the_flow_load_in_the_published_hepatectomies():
    s = hp.summary(draws=4000); res = [x for x in s['states'] if x['shunt'] == 'none']; assert len(res) == 7 == s['resections_without_shunt']
    assert s['flow_load_range'][0] < 1.5 and s['flow_load_range'][1] > 8                                            # a range of flow loads the clinical groups do not have
    for x in res:
        assert x['zR_observed_ci'][0] <= x['zR_observed'] <= x['zR_observed_ci'][1] and x['model']['rigid network'][0] > 0.75
        assert all(abs(x['observed_over_model'][m][0] - x['zR_observed'] / x['model'][m][1]) < 2e-3 for m in hp.MODELS)
        if x['zF'] > 3: assert x['placement']['rigid network']['interval'] == 'below' and x['placement']['keeps distending']['interval'] == 'reaches'
    g16 = [x for x in res if x['study'] == 'Golriz 2016']; assert [x['viscosity_factor'] == 1.0 for x in g16] == [True, True, False] and 0.88 < g16[2]['viscosity_factor'] < 0.95
    low = [x for x in s['states'] if x['shunt'] != 'none' and x['zF'] < 1]; assert low and all(x['zR_observed'] > 1 and x['model']['pressure-flow'][0] > x['model']['rigid network'][1] for x in low)
    f = os.path.join(ROOT, 'results', 'model_hepatectomy.json')
    if os.path.exists(f):
        st = json.load(open(f)); assert [x['zR_observed'] for x in st['states']] == [x['zR_observed'] for x in s['states']] and st['trend']['spearman_rho'] == s['trend']['spearman_rho'] and st['errors']['to_the_range'] == s['errors']['to_the_range']


def test_the_trend_is_tested_on_one_state_per_group_of_animals():
    """Seven states come from five groups of animals: no P value for the seven, an exact one for each choice of one state per group."""
    e = hp.exact_p([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]); assert e['rho'] == -1.0 and abs(e['p_two_sided'] - 2 / 120) < 1e-4 and abs(e['p_one_sided'] - 1 / 120) < 1e-4
    e = hp.exact_p([1, 2, 3, 4, 5], [4, 5, 3, 2, 1]); assert abs(e['rho'] + 0.9) < 1e-9 and abs(e['p_one_sided'] - 5 / 120) < 1e-4
    T = hp.summary(draws=3000)['trend']; assert T['p_value'] is None and T['n'] == 7 and T['groups_of_animals'] == 5 and T['spearman_rho'] < -0.8
    assert len(T['one_state_per_group']) == 3 and all(x['n'] == 5 and x['rho'] <= -0.9 and len(set(s.split(',')[0] for s in x['states'])) == 4 for x in T['one_state_per_group'])     # Bucur 2017 has two groups of animals
    assert T['p_two_sided_one_state_per_group'][0] < 0.02 < T['p_two_sided_one_state_per_group'][1] < 0.1
    for c in ('0', '0.5', '0.8'): assert T['draws'][c]['all_states']['share_of_draws_negative'] > 0.99 and T['draws'][c]['one_state_per_group']['rho'][2] < 0


def test_the_correlation_between_before_and_after_narrows_the_intervals_and_moves_no_state_across_the_rigid_network():
    s = hp.summary(draws=4000); c = s['covariance_sensitivity']; by = c['by_correlation']; assert c['correlations'] == [0.0, 0.5, 0.8]
    assert by['0']['widest_interval_ratio'] > by['0.5']['widest_interval_ratio'] > by['0.8']['widest_interval_ratio'] > 1.5
    for i in range(len(s['states'])):
        w = [by[k]['zR_observed_ci'][i][1] / by[k]['zR_observed_ci'][i][0] for k in ('0', '0.5', '0.8')]; assert w[0] > w[1] > w[2]
    assert c['unchanged_against_the_rigid_network_and_the_two_families'] and all(by[k]['counts']['rigid network']['interval_below'] == 3 for k in by)
    assert all(x['model'] in ('cube law', 'compliance') for x in c['placements_that_change'])
    d = hp.load(); rng = np.random.default_rng(1); g = hp.cm.geometries(); r = d.iloc[0]; b = hp.baseline(rng, 2000)
    a1 = hp.state(r, 2000, np.random.default_rng(5), g, corr=0.0, base=b, private=True)[1]; a2 = hp.state(d.iloc[1], 2000, np.random.default_rng(6), g, corr=0.0, base=b, private=True)[1]
    assert np.corrcoef(a1['zR_d'], a2['zR_d'])[0, 1] > 0.3                                                          # two states of the same animals share the state before


def test_no_error_is_measured_from_the_middle_of_a_range():
    s = hp.summary(draws=2000); e = s['errors']; by = e['by_configuration']; assert 'rms_log_error_of_the_midrange' not in s and 'observed_over_midrange' not in s
    assert all(v[0] <= v[1] for v in by.values()) and by['cube law'][1] < by['pressure-flow'][0] and by['compliance'][1] < by['pressure-flow'][0] and by['pressure-flow'][1] < by['rigid network'][0]
    assert all(e['to_the_range'][m] <= by[m][1] + 1e-9 for m in hp.MODELS) and e['to_the_range']['keeps distending'] <= min(e['to_the_range']['cube law'], e['to_the_range']['compliance']) + 1e-9
    G = hp.cm.geometries(); c = hp.configurations(4.0, 0.25, 1.0, 0.25, 5.0, 4.0, 1.0, G); p = hp.predicted(4.0, 0.25, 1.0, 0.25, 5.0, 4.0, 1.0, G)
    assert all(abs(min(c[m]) - p[m][0]) < 1e-9 and abs(max(c[m]) - p[m][1]) < 1e-9 for m in hp.MODELS) and len(c['rigid network']) == len(G) and len(c['compliance']) == 2 * len(G)
    v = vs.clinical(_d(), draws=2000); assert 'rmse_rigid_midrange' not in v and v['rms_distance_to_the_rigid_range']['blood_of_the_recipient'] < v['rms_distance_to_the_rigid_range']['normal_blood']


def test_a_lower_haematocrit_in_small_vessels_does_not_change_which_groups_the_network_contains():
    assert np.allclose(vs.vessel_haematocrit([1e4, 70.0, 10.0, 5.0], 0.4, 0.5), [0.4, 0.4, 0.2, 0.2]) and abs(float(vs.vessel_haematocrit(26.46, 0.4, 0.5)) - 0.3) < 1e-3
    assert abs(vs.factor(0.30, 0.42, micro=1.0) - vs.factor(0.30, 0.42)) < 1e-12 and vs.factor(0.30, 0.42, micro=0.25) > vs.factor(0.30, 0.42, micro=0.5) > vs.factor(0.30, 0.42)
    m = vs.microvascular(_d()); sc = m['scenarios']; assert list(sc) == ['1', '0.75', '0.5', '0.25'] and m['groups_inside_the_rigid_network'] == [5, 5]
    assert all(x['groups_reaching_the_relations'] == 6 for x in sc.values()) and sc['0.25']['factor'][0] > sc['1']['factor'][0] and sc['0.25']['factor'][1] < 0.95
    err = {k: abs(x['perfused_livers']['0']['observed'] - x['perfused_livers']['0']['law']) for k, x in sc.items()}; assert err['1'] < err['0.75'] < err['0.5'] < err['0.25']      # the perfused livers favour the systemic haematocrit in every vessel
