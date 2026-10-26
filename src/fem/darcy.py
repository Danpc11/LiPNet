"""Perfusion of one lobule as anisotropic Darcy flow.

    div(q) = 0,   q = -(K/mu) grad p,   K = diag(k_r, k_t) in the local polar frame

with p = p_portal on the six triad windows, p = p_central on the central vein, and no flux
on the rest of the hexagon boundary (the symmetry plane shared with the neighbouring
lobule). q is the superficial (Darcy) flux; the interstitial velocity seen by a red cell is
q/phi.

What comes out, and why each piece is here:

    Q           flow through the lobule [m^3/s]        -> validates against 1 mL/min/g
    R_lobule    dP/Q [Pa s m^-3]                       -> this is what LiPNet needs
    tau_transit phi V / Q [s]                          -> validates against 5-10 s
    tau_wall    sinusoidal wall shear stress [Pa]      -> drives eNOS in the NO sub-model

The flow is obtained from the Dirichlet reactions (A p restricted to the constrained nodes),
which is exact at the discrete level: inflow and outflow balance to machine precision.
"""
import numpy as np
from . import fem_core as fc


class DarcySolution(object):
    def __init__(self, mesh, p, u, flow, perm_scale=None):
        self.mesh, self.p, self.u, self.flow = mesh, p, u, flow
        self.perm_scale = perm_scale

    @property
    def speed(self):
        """Darcy speed per element [m/s]."""
        return np.linalg.norm(self.u, axis=1)

    @property
    def interstitial(self):
        """Mean sinusoidal velocity per element [m/s]."""
        return self.speed / self.flow.phi

    @property
    def wall_shear(self):
        """Sinusoidal wall shear stress, Poiseuille estimate tau = 8 mu v / d [Pa]."""
        return 8.0 * self.flow.mu * self.interstitial / self.flow.d_sin

    def transit_time(self):
        return self.flow.phi * self.mesh.volume() / self.Q

    def summary(self):
        from .params import to_mmHg
        tw = self.wall_shear
        return dict(Q_m3_s=self.Q,
                    Q_mL_min=self.Q * 6e7,
                    dP_mmHg=to_mmHg(self.dP),
                    R_lobule=self.R_lobule,
                    transit_s=self.transit_time(),
                    shear_mean_Pa=float(np.average(tw, weights=self.mesh.area)),
                    shear_p90_Pa=float(np.percentile(tw, 90)),
                    v_sinusoid_um_s=float(np.average(self.interstitial, weights=self.mesh.area)) * 1e6)


def solve(mesh, flow, perm_scale=None):
    """Darcy solve. perm_scale is an optional per-element permeability multiplier."""
    A = fc.assemble_aniso(mesh, flow.k_r / flow.mu, flow.k_t / flow.mu, scale=perm_scale)
    b = np.zeros(mesh.n)
    nodes = np.concatenate([mesh.portal, mesh.inner])
    vals = np.concatenate([np.full(mesh.portal.shape, flow.p_portal),
                           np.full(mesh.inner.shape, flow.p_central)])
    Ad, bd = fc.apply_dirichlet(A, b, nodes, vals)
    p = fc.solve_linear(Ad, bd)

    # Darcy flux per element (the sign convention: q points from portal to centre)
    th = np.arctan2(mesh.centroid[:, 1], mesh.centroid[:, 0])
    ct, st = np.cos(th), np.sin(th)
    kxx = flow.k_r * ct * ct + flow.k_t * st * st
    kyy = flow.k_r * st * st + flow.k_t * ct * ct
    kxy = (flow.k_r - flow.k_t) * ct * st
    if perm_scale is not None:
        kxx, kyy, kxy = kxx * perm_scale, kyy * perm_scale, kxy * perm_scale
    gp = np.einsum('ei,eia->ea', p[mesh.tri], mesh.grad)       # grad p per element
    u = -np.column_stack([(kxx * gp[:, 0] + kxy * gp[:, 1]),
                          (kxy * gp[:, 0] + kyy * gp[:, 1])]) / flow.mu

    react = A.dot(p)                                           # nodal fluxes
    Q_in = float(np.sum(react[mesh.portal])) * mesh.geom.H     # per lobule, 2D -> prism
    Q_out = -float(np.sum(react[mesh.inner])) * mesh.geom.H

    sol = DarcySolution(mesh, p, u, flow, perm_scale)
    # Nodal outflow, positive where fluid leaves the lobule. By construction
    # (A p)_i = - int grad w_i . u, so qn_i = -(A p)_i is exactly the boundary term the
    # transport operator needs for a constant to be preserved node by node. Using the
    # geometric edge integral instead leaves an O(h) error at the central vein, where the
    # element-wise Darcy flux jumps most.
    sol.nodal_flux = -react
    sol.Q = 0.5 * (Q_in + Q_out)
    sol.mass_error = abs(Q_in - Q_out) / max(abs(Q_in), 1e-300)
    sol.dP = flow.p_portal - flow.p_central
    sol.R_lobule = sol.dP / sol.Q
    return sol


def annulus_equivalent_resistance(mesh, geom, flow):
    """Resistance of the idealised radial path: an annulus of the same cross-sectional
    area, the same axial permeability and the same central vein, with purely radial flow,

        R = (mu / k_r) ln(b / r_cv) / (2 pi H),    pi b^2 = hexagon area + central vein.

    Dividing the FEM resistance by this isolates what the hexagonal outline, the
    anisotropy and the point drain do, independently of mu and of the permeability scale,
    which is the only part of the lobule that LiPNet can take over: its absolute
    conductance scale is fixed by the cost budget, not by physiology."""
    b = np.sqrt((mesh.cross_area + np.pi * geom.r_cv ** 2) / np.pi)
    return (flow.mu / flow.k_r) * np.log(b / geom.r_cv) / (2.0 * np.pi * geom.H)


def calibrate(mesh, flow, Q_target, transit_target=None):
    """Two observables, two knobs: the permeability sets the flow (Q is linear in k at
    fixed dP) and the porosity sets the transit time. Returns a copy of `flow`."""
    import copy
    f = copy.deepcopy(flow)
    sol = solve(mesh, f)
    scale = Q_target / sol.Q
    f.k_r *= scale
    f.k_t *= scale
    if transit_target is not None:
        f.phi = Q_target * transit_target / mesh.volume()
    return f
