"""The network model applied to the liver (src/model/liver.py) and the audit of the table and fits (src/audit.py)."""
import os, sys
import numpy as np, pandas as pd
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from model import liver, scaling as sc
import indices, audit


def test_canopy_is_the_published_optimum():
    for b in (1.0, 0.75, 0.5):
        a, r = liver.canopy(12, b), sc.canopy_optimum(12, b, 3)
        assert abs(a['D'] / r['D'] - 1) < 1e-10 and abs(a['V'] / r['V'] - 1) < 1e-10
        assert abs(a['share_inlet'] + a['share_sinusoid'] + a['share_outlet'] - 1) < 1e-12


def test_wall_shear_is_not_tau0_unless_b_is_one():
    """D_e / c_e is uniform at the optimum; the Poiseuille shear f / r^3 is tau0 m_e^((b-1)/2)."""
    for b in (1.0, 0.75, 0.5):
        c = liver.canopy(liver.HUMAN_GENERATIONS, b); tau0 = np.sqrt(c['D'])
        assert np.allclose(c['tau'], tau0 * c['m'] ** ((b - 1) / 2))
    c1 = liver.canopy(liver.HUMAN_GENERATIONS, 1.0); assert np.allclose(c1['tau'], c1['tau'][0])       # Murray: uniform shear
    c2 = liver.canopy(liver.HUMAN_GENERATIONS, 0.75); assert c2['tau_sinusoid'] / c2['tau_trunk'] > 5      # not uniform otherwise


def test_diameter_ratio_and_the_exponent_of_each_cast():
    assert abs(liver.diameter_ratio(1.0) - 3 ** (-1 / 3)) < 1e-12                                       # Murray's law, n = 3
    c = liver.canopy(12, 0.75); assert np.allclose(c['r'][1:5] / c['r'][0:4], liver.diameter_ratio(0.75))
    b = {k: liver.exponent_from_cast(**v) for k, v in liver.CASTS.items()}
    assert 0.85 < b['portal vein'] < 1.0 and 0.40 < b['hepatic vein'] < 0.50 and b['hepatic artery'] > 1.0
    assert b['portal vein'] - b['hepatic vein'] > 0.3, 'one exponent does not describe both trees'


def test_a_partial_graft_is_a_subnetwork():
    z = [liver.lobar_graft(liver.HUMAN_GENERATIONS, b, 1)['zR'] for b in (1.0, 0.75, 0.5)]
    assert all(0.90 < x < 1.0 for x in z) and z[0] < z[1] < z[2]
    assert liver.lobar_graft(liver.HUMAN_GENERATIONS, 1.0, 2)['zR'] < z[0], 'a smaller graft leaves more trunk behind'
    G = liver.graft_on_graph(R=4.0, b=1.0, seed=0)
    assert all(0.85 < g['zR'] < 1.0 for g in G['grafts']) and abs(sum(g['fraction_of_lobules'] for g in G['grafts']) - 1) < 1e-9
    assert all(abs(g['sinusoid'] - 1) < 1e-9 for g in G['grafts']), 'a rigid sinusoid keeps its drop per unit flow'


def test_two_inflows():
    assert abs(liver.two_inflows(4.0, 0.45, arterial_load=4.0)['zR'] - 1) < 1e-12           # artery scaling with the portal vein: nothing changes
    assert abs(liver.two_inflows(4.0, 0.0, arterial_load=0.0)['zR'] - (1 - liver.ARTERIAL_FRACTION)) < 1e-12
    assert liver.two_inflows(4.0, 0.45, arterial_load=1.0)['zR'] < 1 and liver.two_inflows(4.0, 1.0, arterial_load=0.0)['zR'] == 1.0


def test_allometry_reports_the_physical_shear():
    r = liver.allometry(1.0, 0.77, 0.0, 0.89)
    assert abs(r['shear_sinusoid'] - r['tau0']) < 1e-9 and abs(r['tau0'] + 0.0766) < 1e-3   # b = 1: the scale is the shear
    r = liver.allometry(0.5, 0.88, 0.10, 0.59)
    assert abs(r['tau0'] - 0.0705) < 1e-3 and r['shear_sinusoid'] < 0 and r['shear_trunk'] < -0.15
    s = liver.summary(graphs=False)['allometry_ranges']
    assert -0.08 < s['shear_sinusoid'][0] and s['shear_sinusoid'][1] < 0 and s['shear_trunk'][0] < -0.25
    assert all(x['flow'] == 'portal flow' and x['lobule_exp'] == 0.10 and x['b'] >= 0.75 for x in s['inside_observed']) and len(s['inside_observed']) == 5


def test_audit_of_the_table_and_the_paired_groups():
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    t = audit.table_checks(d)
    assert t['pct_disagrees_with_events_over_n'] == 0 and t['gradient_disagrees_with_PVP_minus_CVP'] == 0 and t['groups_with_both'] == 6
    p = audit.paired(d); assert len(p['groups']) == 6
    w = [g for g in p['groups'] if g['study'] == 'Wu 2011'][0]; assert abs(w['zR'] - (10.6 / 5) / (275.53 / 84.66)) < 1e-9, 'Wu 2011 is read against its own donor flow'
    for g in p['groups']:
        assert g['zR_ci_sampling'][0] < g['zR'] < g['zR_ci_sampling'][1] < 1
        assert abs(g['donor_gradient_at_which_zR_is_one_mmHg'] - 5 * g['zR']) < 1e-12
    assert max(g['donor_gradient_at_which_zR_is_one_mmHg'] for g in p['groups']) > 3.5, 'the claim needs a donor gradient above this'
    l = audit.loads(d); s = l['share_of_var_log_zP']
    assert abs(s['flow'] + s['resistance'] + s['covariance'] - 1) < 1e-9 and s['flow'] < 0.2 and abs(l['spearman_zF_zP']) < 0.3
    sl = audit.slope(d); u = [r for r in sl['strata'] if r['stratum'].startswith('Uemura')][0]
    assert u['span_mmHg'] < 0.5 and u['own_slope_se'] > 5, 'Uemura 2016 carries no information about the slope'
    m = audit.model_against_data(d); assert m['observed_below_lowest_model_value'] == 6 and 0.75 < m['lowest_model_value'] < 0.9


if __name__ == '__main__':
    for k, v in list(globals().items()):
        if k.startswith('test_'): v(); print('ok', k)


def test_the_passive_relations_span_the_observations_without_a_fitted_number():
    """src/model/distension.py: what the rigid network leaves out, from published relations for the normal liver only."""
    from model import distension as ds
    assert ds.CUBE == 3.0 and ds.DOUBLING_MMHG == (8.0, 13.0) and ds.GREENWAY == dict(intercept=2.2, slope=(0.046, 0.030), control_flow=139.0), 'the published numbers, not tuned'
    assert ds.FLOW_FIXED_SHARE == (0.256, 0.345)
    for law in ds.LAWS:                                   # the donor is a fixed point of every relation, and pressure rises far less than flow
        for p in ds.PARAMS[law]:
            r = ds.distensible(1.0, 1.0, 5.0, law, p); assert abs(r['zR'] - 1) < 1e-9 and abs(r['gradient_mmHg'] - 5) < 1e-9
            assert ds.distensible(4.0, 1.0, 5.0, law, p)['zR'] < ds.distensible(2.0, 1.0, 5.0, law, p)['zR'] < 1
    s = ds.summary(); B = s['distension']['bracket']; fd = s['distension']['flow_doubled_from_the_donor']
    assert s['sinusoids']['zR_floor_relative_to_rigid'][0] > 0.8 and s['outflow']['zR_floor_relative_to_rigid'] == 0.95, 'sinusoids and outlet cannot close the distance alone'
    assert s['required']['reduction_missing'] == [0.28, 0.348], 'against the lowest and the highest network, not the midpoint'
    assert all(1.5 < x < 2.3 for law in ('cube law', 'compliance') for x in fd[law]['rise_mmHg']) and all(3 < x < 4 for x in fd['pressure-flow']['rise_mmHg']), 'about 2 mmHg for a doubling, as the sources state, except the relation that stops distending'
    assert len(B['groups']) == 6 and B['within_the_relations'] == 5 and B['below_the_relations'] == ['Chan 2011 / right lobe with MHV, no modulation'] and not B['above_the_relations']
    assert B['below_the_rigid_network'] == 5 and (B['value_in_a_relation_that_keeps_distending'], B['value_in_the_relation_that_stops_distending']) == (1, 3)
    for g in B['groups']:
        assert g['zR_keeps_distending'][1] < g['zR_stops_distending'][0] < g['zR_stops_distending'][1] < 0.83 and g['zR_rigid'][1] <= 1.0
    band = s['distension']['band_in_the_plane']; assert all(lo <= hi for lo, hi in zip(band['zP_low'], band['zP_high']))


def test_distension_site_by_site_and_the_choices_it_rests_on():
    """src/model/compartments.py: every vessel at its own pressure, only some sites distending, other distending
    pressures, and every uncertainty drawn together."""
    import pandas as pd, indices
    from model import distension as ds, compartments as cm
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    for g in ds.groups(d):                                      # the rigid limit of the site-by-site model is the network of P2
        lo, hi = ds.rigid_range(g['zF']); z = [cm.rigid(g['zF'], x, phi) for x in cm.geometries() for phi in (0.0, 1.0, g['zF'])]
        assert abs(min(z) - lo) < 1e-9 and abs(max(z) - hi) < 1e-9
    donor = cm.geometries()[0][:3] + (0.0,)                     # the whole donor network, no trunk removed
    for site in cm.SITES:                                       # the donor is a fixed point at every site
        for law, p in (('cube law', None), ('compliance', 10.0)):
            assert abs(cm.solve(1.0, donor, 1.0, 5.0, law, p, site)['zR'] - 1.0) < 1e-8
            assert cm.solve(3.5, donor, 1.0, 5.0, law, p, site)['zR'] < cm.rigid(3.5, donor, 1.0), 'a higher flow distends'
    R = ds.robustness(d); S = R['scenarios']
    assert all(v['below'][0].startswith('Chan 2011') for v in S.values()), 'Chan 2011 lies below the relations however they are applied'
    assert all(v['within'] == 5 for k, v in S.items() if k != 'site: postsinusoidal site')
    post = S['site: postsinusoidal site']; assert post['within'] == 4 and post['keeps_distending_span'][0] > 0.6, 'distension of the hepatic venous tree alone does not reach the groups'
    assert all(0.36 < R['share_that_distends']['postsinusoidal site'][i] < 0.49 for i in (0, 1)) and R['share_that_distends']['every vessel'] == [1.0, 1.0]
    J = ds.joint(d, draws=4000)
    assert 0.75 < J['prob_below_the_rigid_network'][0] < 0.85 and J['prob_below_the_rigid_network'][1] > 0.999
    chan = [r for r in J['groups'] if r['study'] == 'Chan 2011'][0]; assert chan['prob_observed_below_every_relation'] > 0.99
    others = [r['prob_observed_within_the_relations'] for r in J['groups'] if r['study'] != 'Chan 2011']; assert all(0.3 < x < 0.75 for x in others)


def test_the_published_relations_hold_in_a_species_they_were_not_derived_from():
    """Pig livers of Ricciardi et al. 2001, numbers as printed. The cube law takes the mean of the pressures on either side
    of the resistance as the distending pressure, as its authors define it."""
    from model import distension as ds
    assert abs(ds.LAWS['cube law'](10.0, 5.0, 5.0) - (7.5 / 10.0) ** 3) < 1e-12, 'Pd = CVP + G / 2, not the portal pressure'
    assert abs(ds.LAWS['compliance'](13.0, 5.0, 5.0, 8.0) - 0.25) < 1e-12, 'volume doubled for 8 mmHg, resistance as volume^-2'
    c = ds.published_check(); st = c['states']
    for state in ('in situ', 'isolated, perfused'):
        m = st[state]['measured_resistance_ratio']; comp = st[state]['predicted']['compliance']['resistance_ratio']; pf = st[state]['predicted']['pressure-flow']['resistance_ratio']
        assert comp[0] - 0.03 <= m <= comp[1], 'the compliance relation holds in normal pig livers'
        assert m < pf[0] - 0.15, 'and the relation that stops distending is too stiff for them'
        assert st[state]['predicted']['compliance']['within_two_se_of_the_measured_ratio'] and not st[state]['predicted']['pressure-flow']['within_two_se_of_the_measured_ratio']
    assert all(v['caval_pressure_consistent'] for v in st.values()), 'the caval pressure implied by the printed values is the same at both flows'
    e = c['largest_portal_pressure_error_mmHg']
    assert e['compliance']['normal_livers'] < 1.2 and e['compliance']['after_cold_ischaemia'] > 2 and e['pressure-flow']['normal_livers'] > 4
    assert all(x < 0 for x in st['after 2 h of cold ischaemia']['predicted']['compliance']['portal_pressure_error_mmHg']), 'after cold ischaemia the liver distends less than the relation gives'
    assert abs(ds.RICCIARDI['after 2 h of cold ischaemia']['hvr'][0] / ds.RICCIARDI['isolated, perfused']['hvr'][0] - 2.0) < 0.01, 'after cold ischaemia the resistance is twice as high at the same flow'
