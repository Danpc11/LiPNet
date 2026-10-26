"""What the lobule FEM gives back to the LiPNet network.

LiPNet represents a whole lobule as ONE edge of conductance w = r^4 / l, i.e. a straight
Poiseuille tube from a portal triad to the central vein. That is a defensible lumped model
for an optimal-network calculation, but it fixes two things it should not:

1. the shape factor. A hexagonal bed fed from its whole perimeter, drained at a point and
   anisotropic in between does not have the resistance of the idealised radial path that
   the single edge assumes. The honest comparison is against an equal-area radial annulus
   with the same permeability,

        R_annulus = (mu / k_r) ln(b/a) / (2 pi H),   pi b^2 = area of the hexagon

   because the absolute conductance scale in LiPNet comes from the cost budget, not from
   physiology: only the RELATIVE correction

        kappa = R_FEM / R_annulus

   is transferable. It is a pure number that depends on the hexagonal outline, the
   anisotropy k_t/k_r and the size of the central vein, and nothing else. It is geometry,
   not fitting.

2. the linearity. A Poiseuille edge has a resistance independent of its flow. A real
   sinusoidal bed does not: raising the flow raises the wall shear stress, which raises
   eNOS output, which dilates the bed and lowers its resistance. That is a resistance that
   falls when flow rises, which is the sign that the z_R < 1 of a hyperperfused graft
   needs. The lobule cannot by itself explain the published z_R values, which come mostly
   from outflow reconstruction, but it says how much of the effect the microcirculation
   can contribute, and it predicts the shape of R(Q) rather than assuming it.

    z_R_lobule(z_F) = R(z_F Q0) / R(Q0)

is computed here by driving the lobule at a prescribed flow instead of a prescribed
pressure, with the shear-NO loop closed at each flow. Feeding it into the network is
`apply_to_graph`.
"""
import copy
import numpy as np
from . import darcy, mesh as mesh_mod, no_wall


class Reference(object):
    """Baseline state: the lobule at its reference pressure drop with no adaptation.
    Everything else is measured against this, so the feedback is inactive at baseline."""

    def __init__(self, params):
        self.params = params
        self.mesh = mesh_mod.build(params.geom)
        self.sol = darcy.solve(self.mesh, params.flow)
        self.no, self.chi = no_wall.field(self.mesh, self.sol, params.no)
        # the reference is elementwise, not a lobule average: the adaptation is a
        # deviation from each element's own baseline operating point, so that the
        # permeability multiplier is exactly 1 everywhere at baseline and z_R(1) = 1.
        self.no_ref = np.maximum(self.no, 1e-300)
        self.shear_ref = self.sol.wall_shear
        self.Q0 = self.sol.Q
        self.R0 = self.sol.R_lobule
        self.dP0 = self.sol.dP


def adapt(ref, dP, params=None):
    """Darcy <-> NO <-> permeability fixed point at a prescribed pressure drop.

    Returns (solution, iterations, converged).
    """
    p = params or ref.params
    fb, no = p.fb, p.no
    flow = copy.deepcopy(p.flow)
    flow.p_portal = flow.p_central + dP
    m = np.ones(len(ref.mesh.area))
    sol = None
    for it in range(fb.max_iter):
        sol = darcy.solve(ref.mesh, flow, perm_scale=m)
        c, _ = no_wall.field(ref.mesh, sol, no)
        target = 1.0 + fb.gain * (c / ref.no_ref - 1.0)
        if fb.injury_shear > 0:
            # measured against the baseline, like the NO term, so that the reference
            # state stays its own reference and z_R(1) = 1 by construction
            over = (np.maximum(sol.wall_shear / fb.injury_shear - 1.0, 0.0)
                    - np.maximum(ref.shear_ref / fb.injury_shear - 1.0, 0.0))
            target = target - fb.injury_gain * over
        target = np.clip(target, fb.k_min, fb.k_max)
        new = (1.0 - fb.relax) * m + fb.relax * target
        err = float(np.max(np.abs(new - m)))
        m = new
        if err < fb.tol:
            return sol, it + 1, True
    return sol, fb.max_iter, False


def solve_at_flow(ref, Q_target, params=None, tol=1e-6, max_iter=40):
    """Drive the lobule at a prescribed flow.

    With the shear-NO loop closed, Q(dP) is strongly superlinear: more pressure gives more
    shear gives more NO gives more permeability gives more flow. A plain secant
    dP <- dP Q*/Q overshoots and oscillates. Q(dP) is monotone and close to a power law,
    so the secant is run on log Q against log dP, with the step limited to a factor e^0.5.
    """
    x = np.log(ref.dP0 * Q_target / ref.Q0)
    sol, _, _ = adapt(ref, float(np.exp(x)), params)
    y = np.log(sol.Q / Q_target)
    if abs(y) < tol:
        return sol
    x_prev, y_prev = x, y
    x = x - 0.3
    for _ in range(max_iter):
        sol, _, _ = adapt(ref, float(np.exp(x)), params)
        y = np.log(sol.Q / Q_target)
        if abs(y) < tol:
            return sol
        dy = y - y_prev
        step = -y * (x - x_prev) / dy if abs(dy) > 1e-14 else -0.1 * np.sign(y)
        x_prev, y_prev = x, y
        x = x + float(np.clip(step, -0.5, 0.5))
    sol.flow_converged = False
    return sol


def resistance_load(params, z_F=(0.5, 1.0, 1.5, 2.0, 3.0, 4.0)):
    """The lobule's own contribution to z_R as a function of the flow load z_F.

    Returns dict with z_F, z_R, shear, and the oxygen delivery per lobule, which is the
    other thing that changes when a graft is hyperperfused.
    """
    ref = Reference(params)
    out = dict(z_F=[], z_R=[], shear_mean_Pa=[], dP_mmHg=[], perm_mean=[])
    from .params import to_mmHg
    for z in z_F:
        sol = solve_at_flow(ref, z * ref.Q0, params)
        out['z_F'].append(float(z))
        out['z_R'].append(float(sol.R_lobule / ref.R0))
        out['shear_mean_Pa'].append(float(np.average(sol.wall_shear, weights=ref.mesh.area)))
        out['dP_mmHg'].append(float(to_mmHg(sol.dP)))
        out['perm_mean'].append(float(np.mean(sol.perm_scale)) if sol.perm_scale is not None else 1.0)
    out['no_feedback_gain'] = params.fb.gain
    out['injury_shear_Pa'] = params.fb.injury_shear
    out['R0_Pa_s_m3'] = ref.R0
    out['Q0_m3_s'] = ref.Q0
    return out


def shape_factor(params):
    """kappa = R_FEM / R_annulus, the correction LiPNet's single sinusoidal edge is missing."""
    ref = Reference(params)
    R_ann = darcy.annulus_equivalent_resistance(ref.mesh, params.geom, params.flow)
    return dict(kappa=float(ref.R0 / R_ann),
                R_fem_Pa_s_m3=float(ref.R0),
                R_annulus_Pa_s_m3=float(R_ann),
                Q0_m3_s=float(ref.Q0),
                Q0_mL_min=float(ref.Q0 * 6e7),
                note='multiply the resistance of every sinusoidal edge of the graph by kappa')


def apply_to_graph(G, kappa, mode='length'):
    """Rescale the sinusoidal edges of a LiPNet graph (model.network.build_liver_graph).

    mode='length'      L_sin <- kappa * L_sin. Resistance and the cost weight a_e = L^(3b/2)
                       both follow, which is the consistent choice if the lobule really is
                       kappa times longer in the hydraulic sense.
    mode='conductance' leaves L alone and returns a multiplier to apply to w on the
                       sinusoidal edges, for when the cost should not move.

    Returns a shallow copy of G; the original is not touched.
    """
    H = dict(G)
    kind = np.asarray(G['kind'])
    sin = kind == 'sin'
    if mode == 'length':
        L = np.asarray(G['L'], dtype=float).copy()
        L[sin] = L[sin] * kappa
        H['L'] = L
    elif mode == 'conductance':
        mult = np.ones(len(kind))
        mult[sin] = 1.0 / kappa
        H['w_mult'] = mult
    else:
        raise ValueError("mode must be 'length' or 'conductance'")
    H['kappa_fem'] = float(kappa)
    return H
