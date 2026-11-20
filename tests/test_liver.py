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


def test_the_distensible_bed_brackets_the_observations_without_a_fitted_number():
    """src/model/distension.py: what the rigid network leaves out, as bounds. Published relations only."""
    from model import distension as ds
    assert ds.CUBE == 3.0 and ds.DOUBLING_MMHG == 8.0, 'the two published relations, not tuned'
    # a rigid bed is recovered when the relation is switched off, and the donor state is a fixed point of both relations
    for law in ds.LAWS:
        r = ds.distensible(1.0, 1.0, 5.0, law); assert abs(r['zR'] - 1) < 1e-9 and abs(r['gradient_mmHg'] - 5) < 1e-9
        hi = ds.distensible(4.0, 1.0, 5.0, law); assert 0.3 < hi['zR'] < 0.6 and hi['zP'] < 2.5, 'pressure rises far less than flow'
        assert ds.distensible(2.0, 1.0, 5.0, law)['zR'] > hi['zR'], 'and the resistance load falls as the flow load rises'
        assert 1.0 < ds.distensible(2.0, 1.0, 5.0, law)['gradient_mmHg'] - 5 < 3.0, 'doubling the flow of a donor liver moves the gradient by a few mmHg, as the source states'
    s = ds.summary(); B = s['distension']['bracket']
    assert s['sinusoids']['zR_floor_relative_to_rigid'][0] > 0.8 and s['outflow']['zR_floor_relative_to_rigid'] == 0.95, 'sinusoids and outlet cannot close the distance alone'
    assert 0.28 < s['required']['reduction_missing'] < 0.36 and s['required']['by_share'][0]['radius_factor'] is None and 1.05 < s['required']['by_share'][-1]['radius_factor'] < 1.15
    assert len(B['groups']) == 6 and (B['between_the_two'], B['at_the_distensible_bed']) == (3, 2) and B['below_the_distensible_bed'] == ['Chan 2011 / right lobe with MHV, no modulation'] and B['above_the_bed_as_point_values'] == 5
    for g in B['groups']:
        assert g['zR_distensible_all'][1] < g['zR_rigid'][0] < g['zR_rigid'][1] <= 1.0, 'the distensible bed lies below the rigid network in every group'
        a, b = g['zR_distensible']['cube law'], g['zR_distensible']['compliance']; assert abs(a[0] - b[0]) < 0.10 and abs(a[1] - b[1]) < 0.10, 'two independent published relations give nearly the same bed'
        assert g['position'] != 'between the two' or g['zR_observed_ci'][0] > g['zR_distensible_all'][1], 'a group counted as between the two has its whole sampling interval above the bed under both relations'
    assert B['mean_distensible_all'][1] < B['mean_observed'] < B['mean_rigid'][0], 'the mean observation lies between the two'
    band = s['distension']['band_in_the_plane']; assert all(lo <= hi for lo, hi in zip(band['zP_low'], band['zP_high'])) and band['zP_high'][-1] < 2.5, 'at five times the flow the pressure load of a distensible bed stays under 2.5'


def test_the_published_relations_hold_in_a_species_they_were_not_derived_from():
    """Pig livers of Ricciardi et al. 2001, numbers as printed. The cube law takes the mean of the pressures above and
    below the resistance site as the distending pressure, as its authors define it."""
    from model import distension as ds
    assert abs(ds.LAWS['cube law'](10.0, 5.0, 5.0) - (7.5 / 10.0) ** 3) < 1e-12, 'Pd = CVP + G / 2, not the portal pressure'
    assert abs(ds.LAWS['compliance'](13.0, 5.0, 5.0) - 0.25) < 1e-12, 'volume doubled for 8 mmHg, resistance as volume^-2'
    c = ds.published_check()
    for state in ('in situ', 'isolated, perfused'):
        m = c[state]['measured_resistance_ratio']; p = [v['resistance_ratio'] for v in c[state]['predicted'].values()]
        assert 0.5 < m < 0.6 and all(0 <= m - x < 0.12 for x in p), 'a normal liver lies at the distensible bed or just above it'
    cold = c['after 2 h of cold ischaemia']; assert cold['resistance_relative_to_the_isolated_liver_at_low_flow'] > 1.1 and abs(ds.RICCIARDI['after 2 h of cold ischaemia']['hvr'][0] / ds.RICCIARDI['isolated, perfused']['hvr'][0] - 2.0) < 0.01, 'after cold ischaemia the resistance is twice as high at the same flow'
    fd = ds.summary()['distension']['flow_doubled_from_the_donor']; assert all(1.5 < v['rise_mmHg'] < 2.2 for v in fd.values()), 'about 2 mmHg for a doubling of portal flow, as published'
