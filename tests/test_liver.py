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
