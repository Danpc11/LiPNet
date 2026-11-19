"""The network model in physical units for a human liver, its growth and the regrowth of a graft (src/model/human.py)."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from model import human, liver
from fem import params as P


def test_the_lobule_and_the_network_hold_the_same_drop():
    f = P.Flow(); assert abs(P.to_mmHg(f.p_portal - f.p_central) - P.LOBULE_DROP_MMHG) < 1e-9
    for b in (0.6, 0.75, 0.93):
        c = human.liver(b)
        assert abs(c['lobule_drop_mmHg'] - P.LOBULE_DROP_MMHG) < 1e-9 and abs(c['dP_mmHg'].sum() - human.GRADIENT_MMHG) < 1e-9
        assert abs(c['share_inlet'] + c['share_sinusoid'] + c['share_outlet'] - 1) < 1e-12
    # perfusion of the lobule layer at that drop against the flow per volume the network sends through a lobule
    from fem import mesh as mesh_mod, darcy
    p = P.default(); m = mesh_mod.Mesh(p.geom); fem = darcy.solve(m, p.flow).Q / m.volume()
    c = human.liver(); network = c['f'][c['g']] / (c['territory_side_mm'] * 1e-3) ** 3
    assert 0.7 < fem / network < 1.0, 'both layers carry the hepatic flow per volume at that drop, within the perfusion target of the lobule layer'


def test_physical_canopy_is_the_canopy_of_liver_py():
    """With the lobule left to the optimum and one flow in both trees it is liver.canopy in other units."""
    for b in (1.0, 0.75, 0.5):
        a = liver.canopy(liver.HUMAN_GENERATIONS, b); c = human.liver(b, lobule_drop=None, arterial_fraction=0.0)
        assert abs(a['share_inlet'] - c['share_inlet']) < 1e-9 and abs(a['share_sinusoid'] - c['share_sinusoid']) < 1e-9
        assert np.allclose(a['r'][1:] / a['r'][:-1], c['r_mm'][1:] / c['r_mm'][:-1])
        # Poiseuille recomputed from radii, lengths and flows returns the gradient
        dP = 8 * human.VISCOSITY * c['l_mm'] * 1e-3 * c['f'] / (np.pi * (c['r_mm'] * 1e-3) ** 4) / P.MMHG
        assert abs(dP.sum() - human.GRADIENT_MMHG) < 1e-9


def test_portal_radii_of_the_cast_fix_the_exponent():
    e = human.portal_exponent()
    assert e['interval'][0] <= human.B_HUMAN <= e['interval'][1] and 0.65 < e['interval'][0] and e['interval'][1] < 0.85
    assert e['rms_at_half'] > 0.25 and e['rms_at_one'] > 0.20, 'neither end of [0.5, 1] reproduces the cast'
    r = human.portal_radii(human.B_HUMAN); q = r['model_over_cast']; assert len(q) == 11 and 0.88 < min(q) and max(q) < 1.12
    assert human.portal_exponent(kappa=1.0)['interval'] is None, 'branches as long as their territory do not reproduce it at any b'
    assert 0.60 < human.portal_exponent(g=13)['b'] < 0.80, 'a third of the lobules moves it little'
    c = human.cast_exponent('portal vein'); assert 0.70 < c['b'] < 0.76 and c['leave_one_out'][1] - c['leave_one_out'][0] < 0.06     # eleven generations
    c = human.cast_exponent('portal vein', 'meso'); assert 0.85 < c['b'] < 0.93 and 0.30 < min(c['length_over_territory_side']) and max(c['length_over_territory_side']) < 0.55
    assert 0.95 < human.cast_exponent('hepatic artery')['b'] < 1.10 and human.cast_exponent('hepatic vein')['b'] < 0.5
    a = human.liver(); assert 45 < a['portal_trunk_length_mm'] < 60 and 6.5 < a['portal_trunk_radius_mm'] < 7.5 and 10 < a['portal_trunk_velocity_cm_s'] < 20


def test_the_hepatic_veins_of_the_cast_are_not_those_of_a_perfused_liver():
    pv, hv = human.cast_poiseuille('portal vein'), human.cast_poiseuille('hepatic vein')
    meso = human.cast_poiseuille('portal vein', part='meso')
    assert 0.9 < meso['drop_ratio_per_generation'] < 1.1 and 0.6 < meso['drop_measured_generations'] < 0.9 and 1.0 < pv['drop_measured_generations'] < 1.3
    assert hv['generations'] == [2, 3, 4, 5, 6, 7, 8] and hv['drop_ratio_per_generation'] > 1.4 and hv['drop_measured_generations'] > 3.0
    assert all(1.5 < human.cast_tree_drop('portal vein', G) < 2.5 for G in human.CAST_GENERATIONS['portal vein'])
    assert all(human.cast_tree_drop('hepatic vein', G, 1.2) > human.GRADIENT_MMHG for G in human.CAST_GENERATIONS['hepatic vein'])


def test_growth_with_the_body():
    for b in (0.6, 0.75, 0.93):
        for s in (0.0, 0.5, 1.0):
            g = human.growth(b, s); e = g['exponents']
            assert 0.33 < e['portal_diameter_mm'] < 0.41 and 0.20 < e['portal_velocity_cm_s'] < 0.33, 'the trunk does not tell how the liver grows'
            assert abs(e['lobule_drop_mmHg'] - 2 * s / 3) < 1e-6 and abs(e['sinusoidal_velocity'] - s / 3) < 1e-6
            assert abs(g['rows'][-1]['lobule_drop_mmHg'] - P.LOBULE_DROP_MMHG) < 1e-9
    m = human.measured_growth()
    assert m['portal_diameter']['interval'][0] < 0.37 < m['portal_diameter']['interval'][1]
    assert abs(human.urata_volume(70) - 706.2 * 287 / 160 - 2.4) < 1e-9


def test_regrowth_of_a_graft():
    for cal in ('unchanged', 'proportional'):
        r = human.regeneration(1.0, zF0=3.3, zR0=0.55, calibre=cal); assert abs(r['zF'] - 3.3) < 1e-12 and abs(r['zR'] - 0.55) < 1e-12
    u, p = human.regeneration(1.58, calibre='unchanged'), human.regeneration(1.58, calibre='proportional')
    assert abs(u['zF'] - 1 / 1.58) < 1e-12 and u['zP'] > 1 > p['zP'] and u['zR'] > p['zR'] > 1
    assert abs(float(u['zP']) - (0.8 * 1.58 ** (1 / 3) + 0.2 * 1.58 ** (-1 / 3))) < 1e-9          # lobule share of one fifth
    r = human.regeneration(3.0, zF0=3.0, calibre='proportional')
    assert abs(float(r['zF']) - 1) < 1e-12 and 1.15 < float(r['zP']) < 1.30 and abs(float(r['sinusoidal_velocity']) - 3 ** (1 / 3)) < 1e-12


def test_two_inflows_measured():
    m = human.two_inflows_measured(); w = m['Wu 2011']
    assert abs(w['flow_load'] - 275.53 / 84.66) < 1e-9 and abs(w['arterial_load'] - (157.69 / 1713.98 * 275.53) / (237.89 / 527.81 * 84.66)) < 1e-9
    assert 0.60 < w['arterial_load'] < 0.70 and 0.28 < w['arterial_share_donor'] < 0.34 and abs(w['zR_observed'] - 2.12 / w['flow_load']) < 1e-9
    assert w['zR_observed'] < w['zR_model_all_resistance_after_the_junction'] < w['zR_model'][0] < w['zR_model'][1] < 0.9
    v = m['Garcia-Valdecasas 2003']; assert 0.75 < v['arterial_load'] < 0.85 and 'zR_observed' not in v and 0.85 < v['zR_model'][0]


def test_the_first_week():
    w = human.first_week()
    assert abs(w['unchanged']['pressure_load'] - w['inflow'] * (0.8 * w['volume'] ** (1 / 3) + 0.2 * w['volume'] ** (-1 / 3))) < 1e-9
    assert w['proportional']['pressure_load'] < 0.5 < w['unchanged']['pressure_load'] < 0.8
    lo = min(w['observed_gradient_if_caval_pressure_8'].values()); assert lo > w['unchanged']['pressure_load'] > w['proportional']['pressure_load']
    assert w['resistance_implied'][0] > w['unchanged']['resistance'] > 1 > w['proportional']['resistance']
    r = human.regeneration(2.0, inflow=0.5, calibre='unchanged'); assert abs(float(r['zF']) - 0.25) < 1e-12
    d = [x[3] for x in human.GRAFT_VOLUME if x[2] <= 15]; assert 1.5 < min(d) and max(d) < 2.0


if __name__ == '__main__':
    for k, v in sorted(globals().items()):
        if k.startswith('test_'): v(); print('ok   ', k)
