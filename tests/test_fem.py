"""Checks that matter for the lobule FEM: the discretisation is right (patch test), the
solver conserves what it should, the fields land where physiology says they land, and the
coupling back to the network reduces to the rigid Poiseuille edge when the feedback is off.

    python -m pytest tests/test_fem.py -q      or      python tests/test_fem.py
"""
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))

from fem import params as P
from fem import mesh as mesh_mod, fem_core as fc, darcy, transport as tr, no_wall, coupling


def small(nr=24, nt=36, **kw):
    p = P.default()
    p.geom.nr, p.geom.nt = nr, nt
    for k, v in kw.items():
        setattr(p.geom, k, v)
    return p, mesh_mod.build(p.geom)


# ------------------------------------------------------------------ mesh
def test_mesh_is_valid():
    p, m = small()
    assert np.all(m.area > 0), 'every element must be positively oriented'
    assert m.n == (m.nr + 1) * m.nt
    assert m.nt % 6 == 0, 'the six portal corners must fall on nodes'
    hexagon = 3 * np.sqrt(3) / 2 * p.geom.R ** 2 - np.pi * p.geom.r_cv ** 2
    assert abs(m.cross_area / hexagon - 1) < 0.02


def test_boundary_rings_are_consistent():
    p, m = small()
    for ring, sign in ((m.ring_outer, +1.0), (m.ring_inner, -1.0)):
        for k in range(0, len(ring['a']), 7):
            t = m.tri[ring['elem'][k]]
            assert ring['a'][k] in t and ring['b'][k] in t, 'edge not in its element'
            mid = 0.5 * (m.xy[ring['a'][k]] + m.xy[ring['b'][k]])
            assert sign * np.dot(mid, ring['normal'][k]) > 0, 'normal points the wrong way'
    assert abs(m.ring_outer['length'].sum() / (6 * p.geom.R) - 1) < 1e-9
    assert abs(m.ring_inner['length'].sum() / (2 * np.pi * p.geom.r_cv) - 1) < 1e-2


# ------------------------------------------------------------------ discretisation
def test_patch_test_reproduces_a_linear_field():
    """P1 elements must reproduce a linear field exactly, which is what fails loudly if
    the element gradients or the assembly are wrong. It is run with an isotropic tensor on
    purpose: the anisotropic K of the lobule rotates with the polar angle, so
    div(K(x) grad p) does not vanish for a linear p and there would be nothing exact to
    check against."""
    p, m = small()
    exact = 3.0 + 1e4 * m.xy[:, 0] - 7e3 * m.xy[:, 1]
    A = fc.assemble_aniso(m, 1.0, 1.0)
    bnd = np.concatenate([m.outer, m.inner])
    Ad, bd = fc.apply_dirichlet(A, np.zeros(m.n), bnd, exact[bnd])
    num = fc.solve_linear(Ad, bd)
    assert np.max(np.abs(num - exact)) / np.max(np.abs(exact)) < 1e-10


def test_darcy_conserves_mass_and_is_linear_in_pressure():
    import copy
    p, m = small()
    s = darcy.solve(m, p.flow)
    assert s.mass_error < 1e-9, 'inflow and outflow must balance to machine precision'
    f2 = copy.deepcopy(p.flow)
    f2.p_portal = f2.p_central + 2.0 * s.dP
    s2 = darcy.solve(m, f2)
    assert abs(s2.Q / s.Q - 2.0) < 1e-9, 'Darcy without feedback is linear in dP'


def test_perfusion_is_physiological():
    p, m = small(48, 72)
    s = darcy.solve(m, p.flow)
    per_g = s.Q * 6e7 / (m.volume() * 1e6)          # mL/min per g, tissue density ~1
    assert 0.5 < per_g < 2.0, 'perfusion %.2f mL/min/g is off scale' % per_g
    assert 3.0 < s.transit_time() < 15.0
    mean_shear = float(np.average(s.wall_shear, weights=m.area))
    assert 0.05 < mean_shear < 0.6, 'sinusoidal shear %.3f Pa is off scale' % mean_shear


# ------------------------------------------------------------------ transport
def test_a_constant_is_preserved_without_reaction():
    """With no reaction and a uniform inlet the only solution is the constant. Any
    spurious source in the advection operator shows up here."""
    p, m = small()
    s = darcy.solve(m, p.flow)
    react = lambda c: (np.zeros_like(c), np.zeros_like(c))
    c = fc.solve_adr(m, s.u, 1e-9, react, [(m.portal, 1.0)], c0=np.zeros(m.n),
                     qn=s.nodal_flux)
    assert np.max(np.abs(c - 1.0)) < 1e-10, 'the nodal-flux boundary term must be exact'


def test_solutes_stay_between_zero_and_the_inlet():
    p, m = small(32, 48)
    s = darcy.solve(m, p.flow)
    for res, cin in ((tr.ethanol(m, s, p.etoh), p.etoh.c_in),
                     (tr.drug(m, s, p.drug), p.drug.c_in),
                     (tr.oxygen(m, s, p.o2), p.o2.p_in)):
        assert res.c_min > -1e-9, '%s undershoots' % res.name
        assert res.c_max <= cin * (1 + 1e-9), '%s overshoots its inlet' % res.name


def test_mass_balance_improves_with_the_mesh():
    """influx = outflux + consumption holds exactly by construction, so what is checked
    here is the independent geometric evaluation of the same influx, which is only O(h).
    It must converge: this is what caught the non-conservative advection form, where it
    diverged under refinement instead."""
    err = []
    for nr, nt in ((24, 36), (48, 72), (96, 144)):
        p, m = small(nr, nt)
        s = darcy.solve(m, p.flow)
        err.append(tr.oxygen(m, s, p.o2).balance_error)
    assert err[0] > err[1] > err[2], 'balance error must fall under refinement: %s' % err
    assert err[-1] < 0.02


def test_oxygen_zonation_matches_physiology():
    p, m = small(48, 72)
    s = darcy.solve(m, p.flow)
    o = tr.oxygen(m, s, p.o2)
    assert 50 < o.c_max <= p.o2.p_in
    assert 25 < o.c_cv < 40, 'pericentral pO2 %.1f mmHg outside 25-40' % o.c_cv
    assert 0.15 < o.extraction < 0.45, 'O2 extraction %.2f outside 15-45%%' % o.extraction
    assert o.sat_in > o.sat_out > 0.4
    z, prof = m.radial_profile(o.c, 12)
    assert np.all(np.diff(prof) < 1e-9), 'pO2 must fall monotonically towards the centre'


def test_pericentral_enzymes_put_the_damage_pericentrally():
    p, m = small(32, 48)
    s = darcy.solve(m, p.flow)
    et = tr.ethanol(m, s, p.etoh)
    ac = tr.acetaldehyde(m, s, et, p.acet)
    dr = tr.drug(m, s, p.drug)
    assert ac.peak_zeta > 0.7, 'acetaldehyde should peak pericentrally'
    assert dr.reactive_peak_zeta > 0.7, 'CYP-derived metabolite should peak pericentrally'
    assert 0.0 < et.cyp_fraction < 0.5, 'at 20 mM most ethanol still goes through ADH'
    assert ac.extraction > 0.8, 'most acetaldehyde is cleared before leaving the lobule'


# ------------------------------------------------------------------ two phases
def test_two_phase_reduces_to_one_phase_when_the_membrane_is_open():
    """With Da_x >> 1 the hepatocyte equilibrates with the blood and the two-continuum
    model must return the single-continuum answer. The residual difference is the SUPG
    parameter, which is built on the vascular diffusivity phi D rather than on D, so it is
    an O(h) stabilisation difference and not a modelling one."""
    from fem import two_phase as tp
    p, m = small(32, 48)
    s = darcy.solve(m, p.flow)
    one = tr.oxygen(m, s, p.o2)
    grads = []
    for perm in (2e-3, 2e-2, 2e-1):
        two = tp.oxygen(m, s, p.o2, perm)
        assert two.Da_exchange > 1e3
        assert abs(two.extraction / one.extraction - 1) < 0.01
        grads.append(abs(two.gradient))
    assert grads[0] > grads[1] > grads[2], 'the blood-tissue gap must close as P rises'
    assert grads[2] / grads[0] < 0.02, 'and it must close like 1/P'


def test_a_tight_membrane_limits_clearance_not_vmax():
    """The point of the module: below Da_x ~ 1 the lobule clears what the membrane
    delivers, whatever the enzyme could do."""
    from fem import two_phase as tp
    p, m = small(32, 48)
    s = darcy.solve(m, p.flow)
    sweep = tp.permeability_sweep(m, s, p.drug, [1e-3, 1e-5, 1e-6, 1e-7])
    ext = np.array(sweep['extraction'])
    ratio = np.array(sweep['hepatocyte_over_blood'])
    assert np.all(np.diff(ext) <= 1e-9), 'extraction must fall as the membrane tightens'
    assert np.all(np.diff(ratio) < 0), 'so must the hepatocyte-to-blood ratio'
    assert ratio[0] > 0.99 and ratio[-1] < 0.8
    assert ext[-1] < 0.85 * ext[0], 'the tight membrane must cost real clearance'
    fast = np.array(sweep['Da_exchange'])
    assert fast[0] > 100 > 1 > fast[-1], 'the sweep must cross Da_x = 1'


def test_two_phase_conserves():
    from fem import two_phase as tp
    p, m = small(32, 48)
    s = darcy.solve(m, p.flow)
    r = tp.drug(m, s, p.drug, 1e-6)
    assert r.balance_error < 0.05
    assert r.c_v.min() > -1e-9 and r.c_t.min() > -1e-9
    assert r.c_v.max() <= p.drug.c_in * (1 + 1e-9)


def test_a_matched_membrane_clears_more_than_a_mismatched_one():
    """Zonal fenestration matters only through its alignment with the enzyme, and only at
    the margin: the diagonal wins, but by about 10%, against the five-fold effect of the
    permeability level itself."""
    from fem import two_phase as tp
    p, m = small(32, 48)
    s = darcy.solve(m, p.flow)
    mat = np.array(tp.zonation_matrix(m, s, p.drug, 1e-7)['extraction'])
    for i in range(3):
        assert np.argmax(mat[i]) == i, 'row %d: the matched membrane should win' % i
    spread = (mat.max() - mat.min()) / mat.mean()
    assert spread < 0.25, 'zonation is a second-order effect, not a first-order one'


# ------------------------------------------------------------------ the network
def test_sinusoidal_share_is_a_minority_and_shrinks():
    from fem import network_link as nl
    sc = nl.share_scaling((3.0, 4.5))
    sh = [r['sinusoidal_share'] for r in sc['rows']]
    assert all(0.0 < x < 0.5 for x in sh), 'share out of range: %s' % sh
    assert sh[1] < sh[0], 'the share must fall as the tree deepens'
    assert sc['log_slope'] < 0


def test_a_rigid_lobule_leaves_the_graft_indices_untouched():
    from fem import network_link as nl
    for h in (1.0, 2.0, 4.0):
        g = nl.graft(h, 0.25, None)
        assert abs(g['z_R_graft'] - 1.0) < 1e-9
        assert abs(g['zF'] - h) < 1e-9


def test_only_the_sinusoidal_share_of_the_lobule_effect_is_measured():
    """1 - z_R_graft must be the sinusoidal share of 1 - z_R_lobule, weighted by the part
    of the series the measurement spans. This is the whole point of the module."""
    from fem import network_link as nl
    load = dict(z_F=[1.0, 2.0, 4.0], z_R=[1.0, 0.75, 0.5])
    curve = nl.zr_interpolator(load)
    assert abs(curve(2.0) - 0.75) < 1e-9
    r_in = r_out = 0.05
    r_lob = 1 - r_in - r_out
    for share in (0.1, 0.25, 0.5):
        g = nl.graft(4.0, share, curve, r_in=r_in, r_out=r_out)
        expected = share * r_lob / (r_lob + r_out)
        got = (1 - g['z_R_graft']) / (1 - g['z_R_lobule'])
        assert abs(got - expected) < 1e-6, 'share %.2f: %.4f vs %.4f' % (share, got, expected)


def test_the_injury_branch_flips_the_sign():
    """The porcine small-for-size data report a higher portal resistance, not a lower one.
    With the injury branch on, the lobule must cross to z_R > 1, and the reference state
    must still be its own reference."""
    p, _ = small()
    p.fb.injury_shear = 0.4
    out = coupling.resistance_load(p, z_F=(1.0, 2.0, 4.0))
    zr = np.array(out['z_R'])
    assert abs(zr[0] - 1.0) < 1e-6
    assert np.all(zr[1:] > 1.0), 'injury must raise the resistance: %s' % zr
    assert np.all(np.diff(zr) > 0)


# ------------------------------------------------------------------ nitric oxide
def test_no_is_a_boundary_layer_not_a_field():
    p, _ = small()
    assert np.sqrt(p.no.D / p.no.k_rbc) < 1e-6, 'lumen diffusion length must be sub-micron'
    rc, cq = no_wall.radial_profile(p.no, p.flow)
    a = 0.5 * p.flow.d_sin
    lumen = np.interp(0.0, rc, cq)
    wall = np.interp(a, rc, cq)
    assert lumen < 1e-4 * wall, 'haemoglobin must annihilate NO in the lumen'
    tissue = cq[rc > a + p.no.wall]
    assert np.all(np.diff(tissue) < 1e-12), 'NO must decay monotonically into the tissue'
    assert no_wall.geometric_factor(p.no, p.flow) > 0


def test_no_follows_shear_and_shear_follows_flow():
    import copy
    p, m = small()
    s1 = darcy.solve(m, p.flow)
    f2 = copy.deepcopy(p.flow)
    f2.p_portal = f2.p_central + 2.0 * s1.dP
    s2 = darcy.solve(m, f2)
    c1, _ = no_wall.field(m, s1, p.no)
    c2, _ = no_wall.field(m, s2, p.no)
    ratio = float(np.average(c2, weights=m.area) / np.average(c1, weights=m.area))
    assert abs(ratio - 2.0 ** p.no.shear_gain) < 1e-6


# ------------------------------------------------------------------ coupling
def test_without_feedback_the_lobule_is_a_rigid_resistor():
    p, _ = small()
    p.fb.gain = 0.0
    out = coupling.resistance_load(p, z_F=(1.0, 2.0, 4.0))
    assert np.allclose(out['z_R'], 1.0, atol=1e-6), 'gain 0 must reproduce the Poiseuille edge'


def test_feedback_makes_resistance_fall_with_flow():
    p, _ = small()
    out = coupling.resistance_load(p, z_F=(1.0, 2.0, 3.0))
    zr = np.array(out['z_R'])
    assert abs(zr[0] - 1.0) < 1e-6, 'the reference state must be its own baseline'
    assert np.all(np.diff(zr) < 0), 'z_R must fall as z_F rises: %s' % zr
    assert zr[-1] > 0.3, 'the clamp should keep the effect physiological'
    shear = np.array(out['shear_mean_Pa'])
    assert np.allclose(shear / shear[0], out['z_F'], rtol=1e-3), \
        'at prescribed flow the shear scales with z_F whatever the calibre does'


def test_shape_factor_is_order_one():
    p, _ = small()
    k = coupling.shape_factor(p)
    assert 0.3 < k['kappa'] < 3.0, 'kappa %.3f is not a shape factor' % k['kappa']


def test_graph_patch_touches_only_sinusoids():
    kind = np.array(['in', 'sin', 'out', 'sin'])
    G = dict(kind=kind, L=np.array([2.0, 3.0, 4.0, 5.0]))
    H = coupling.apply_to_graph(G, 1.5)
    assert np.allclose(H['L'], [2.0, 4.5, 4.0, 7.5])
    assert np.allclose(G['L'], [2.0, 3.0, 4.0, 5.0]), 'the original graph must not move'
    H2 = coupling.apply_to_graph(G, 1.5, mode='conductance')
    assert np.allclose(H2['w_mult'], [1.0, 1 / 1.5, 1.0, 1 / 1.5])


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
