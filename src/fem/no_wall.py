"""Nitric oxide as a wall sub-model, not as a lobule-scale field.

Why it is not solved on the lobule mesh. With haemoglobin in the sinusoidal lumen the
scavenging rate is k[Hb] ~ 5e7 M^-1 s^-1 x 1e-2 M ~ 5e5 s^-1, so the diffusion length
sqrt(D/k) is about 0.1 um. Even the perivascular tissue, with k ~ 1e2 s^-1, gives
sqrt(D/k) ~ 6 um, which is smaller than one element of a lobule mesh. An advection-
diffusion-reaction field for NO on the lobule mesh would be resolving nothing: the answer
is set entirely by a boundary layer around each sinusoid, so that is what gets solved.

The sub-model is one cylinder per element:

    0 < r < a - delta      red-cell core, scavenging k_rbc
    a - delta < r < a      plasma free layer, no scavenging (this layer is most of the
                           reason NO reaches the tissue at all)
    a < r < a + w          endothelium: eNOS releases q = q0 (tau_w / tau_ref)^gain
    a + w < r < b          tissue, scavenging k_tissue, b = a / sqrt(phi) from the
                           sinusoidal volume fraction, zero flux at r = b by symmetry

Because every consumption term is first order, the whole problem is LINEAR in the
production q. So the radial solve is done once, giving a geometric factor

    chi = [NO](r_target) / q        [mol m^-3 / (mol m^-2 s^-1)]

and the NO field over the lobule is then chi * q0 * (tau_w/tau_ref)^gain evaluated element
by element. One tridiagonal solve for the whole lobule instead of one per element.
"""
import numpy as np
from scipy.linalg import solve_banded


def radial_profile(no, flow):
    """Solve the 1D cylinder for unit production. Returns (r, C/q)."""
    a = 0.5 * flow.d_sin
    b = a / np.sqrt(flow.phi)
    n = int(no.n_radial)
    rf = np.linspace(0.0, b, n + 1)                 # cell faces
    rc = 0.5 * (rf[1:] + rf[:-1])                   # cell centres
    vol = 0.5 * (rf[1:] ** 2 - rf[:-1] ** 2)        # per unit length / (2 pi)

    k = np.where(rc < a - no.free_layer, no.k_rbc,
                 np.where(rc < a, 0.0, no.k_tissue))

    shell = (rc >= a) & (rc < a + no.wall)
    if not np.any(shell):
        shell = np.argmin(np.abs(rc - a)) == np.arange(n)
    src = np.zeros(n)
    src[shell] = a / np.sum(vol[shell])             # unit flux at r = a, spread in the shell

    lo, di, up = np.zeros(n), np.zeros(n), np.zeros(n)
    dr = rc[1] - rc[0]
    for i in range(n):
        di[i] = k[i] * vol[i]
        if i > 0:
            t = no.D * rf[i] / dr
            di[i] += t
            lo[i] = -t
        if i < n - 1:
            t = no.D * rf[i + 1] / dr
            di[i] += t
            up[i] = -t
    ab = np.zeros((3, n))
    ab[0, 1:] = up[:-1]
    ab[1, :] = di
    ab[2, :-1] = lo[1:]
    c = solve_banded((1, 1), ab, src * vol)
    return rc, c


def geometric_factor(no, flow):
    """chi = [NO] at the target radius per unit eNOS flux [s / m]."""
    rc, cq = radial_profile(no, flow)
    a = 0.5 * flow.d_sin
    rt = a + no.wall + no.r_target
    return float(np.interp(rt, rc, cq))


def field(mesh, sol, no):
    """Perivascular [NO] per element [mol/m^3], driven by the local wall shear stress."""
    chi = geometric_factor(no, sol.flow)
    tau = sol.wall_shear
    q = no.q0 * (np.maximum(tau, 0.0) / no.tau_ref) ** no.shear_gain
    return chi * q, chi


def summary(mesh, sol, no):
    c, chi = field(mesh, sol, no)
    w = mesh.area
    return dict(no_chi_s_per_m=chi,
                no_mean_nM=float(np.average(c, weights=w)) * 1e6,
                no_periportal_nM=float(np.average(c[mesh.zeta[mesh.tri].mean(axis=1) < 0.25],
                                                  weights=w[mesh.zeta[mesh.tri].mean(axis=1) < 0.25])) * 1e6,
                no_pericentral_nM=float(np.average(c[mesh.zeta[mesh.tri].mean(axis=1) > 0.75],
                                                   weights=w[mesh.zeta[mesh.tri].mean(axis=1) > 0.75])) * 1e6)
