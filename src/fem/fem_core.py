"""P1 finite elements on the lobule mesh: anisotropic Darcy and a SUPG-stabilised,
Newton-solved advection-diffusion-reaction operator.

Two routines do all the work.

assemble_aniso(mesh, k_r, k_t)
    stiffness of -div(K grad p) with K = k_r e_r e_r^T + k_t e_th e_th^T, rotated at the
    centroid of each element.

solve_adr(mesh, u, diff, react, dirichlet, g, dg)
    steady   div(u g(c)) - div(diff grad c) + react(c) = 0,   strong Dirichlet data.

    g is the advective capacity: identity for a free solute, total O2 (free plus the pool
    bound to haemoglobin) for oxygen, which is what makes the oxygen problem nonlinear
    even though its reaction is a plain Michaelis-Menten.

    The advection is assembled in CONSERVATIVE form,

        int_K div(u g) w  =  - int_K grad w . u g  +  (boundary term),

    rather than as u . grad g. The two are equivalent only if the discrete velocity is
    pointwise divergence free, which a P1 Darcy flux is not: it is divergence free in the
    weak sense, with jumps across faces. With u . grad g the mass lost scales with the
    magnitude of g rather than with its variation, and for oxygen g is about 8 mol/m^3
    while its variation across the lobule is about 2, so the balance error grows instead
    of converging under refinement.

    The boundary term is the nodal Darcy outflow qn = -(A p), not a geometric edge
    integral. Per element -int grad phi_i . u + int_dK phi_i u.n = 0 exactly, but summing
    that over elements leaves the jumps of u across interior faces, so an edge integral
    preserves a constant only to O(h), with the error piled up at the central vein where
    the flux jumps most. The nodal flux is by definition the quantity that cancels
    -int grad w_i . u node by node, so with it a uniform inlet and no reaction give the
    constant solution to machine precision, and the lobule balance
    influx = outflux + consumption is an identity rather than an approximation.

Reactions use a lumped mass matrix, which keeps the operator monotone at the Damkohler
numbers of the drug and ethanol problems. SUPG weight
tau = (h / 2|u|)(coth Pe - 1/Pe), Pe = |u| h / (2 diff); with P1 the second derivative in
the strong residual vanishes elementwise, so the stabilisation only sees advection and
reaction.
"""
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import spsolve


# ---------------------------------------------------------------- Darcy / diffusion
def assemble_aniso(mesh, k_r, k_t, scale=None):
    """Stiffness of -div(K grad .). scale is an optional per-element multiplier."""
    th = np.arctan2(mesh.centroid[:, 1], mesh.centroid[:, 0])
    ct, st = np.cos(th), np.sin(th)
    kxx = k_r * ct * ct + k_t * st * st
    kyy = k_r * st * st + k_t * ct * ct
    kxy = (k_r - k_t) * ct * st
    if scale is not None:
        kxx, kyy, kxy = kxx * scale, kyy * scale, kxy * scale
    G = mesh.grad
    KG = np.empty_like(G)
    KG[:, :, 0] = kxx[:, None] * G[:, :, 0] + kxy[:, None] * G[:, :, 1]
    KG[:, :, 1] = kxy[:, None] * G[:, :, 0] + kyy[:, None] * G[:, :, 1]
    Ke = np.einsum('eia,eja->eij', G, KG) * mesh.area[:, None, None]
    return scatter(mesh, Ke)


def scatter(mesh, Ke):
    """Element 3x3 blocks -> global sparse matrix."""
    rows = np.repeat(mesh.tri, 3, axis=1).ravel()
    cols = np.tile(mesh.tri, (1, 3)).ravel()
    return coo_matrix((Ke.ravel(), (rows, cols)), shape=(mesh.n, mesh.n)).tocsr()


def apply_dirichlet(A, b, nodes, values):
    """Strong Dirichlet rows. Returns new (A, b)."""
    A = A.tolil(copy=True)
    nodes = np.asarray(nodes)
    for i in nodes:
        A.rows[i] = [i]
        A.data[i] = [1.0]
    b = np.asarray(b, dtype=float).copy()
    b[nodes] = values
    return A.tocsr(), b


def solve_linear(A, b):
    return spsolve(csr_matrix(A), b)


# ---------------------------------------------------------------- transport
def supg_tau(mesh, u, diff):
    """(M,) stabilisation parameter."""
    un = np.linalg.norm(u, axis=1)
    out = np.zeros_like(un)
    m = un > 1e-30
    if not np.any(m):
        return out
    if diff <= 0:
        out[m] = mesh.h[m] / (2.0 * un[m])
        return out
    pe = un[m] * mesh.h[m] / (2.0 * diff)
    xi = np.where(pe < 1e-6, pe / 3.0,
                  1.0 / np.tanh(np.maximum(pe, 1e-12)) - 1.0 / np.maximum(pe, 1e-12))
    out[m] = mesh.h[m] / (2.0 * un[m]) * xi
    return out



def solve_adr(mesh, u, diff, react, dirichlet, g=None, dg=None, c0=None,
              supg=True, qn=None, tol=1e-9, max_iter=50, verbose=False):
    """Steady advection-diffusion-reaction by damped Newton.

    Parameters
    ----------
    u         : (M, 2) Darcy (superficial) velocity per element [m/s]
    diff      : scalar diffusivity, in the units of the unknown
    react     : callable c -> (R, dR/dc), nodal arrays; +R is a sink
    dirichlet : list of (node array, value) pairs
    g, dg     : advective capacity and its derivative (default: identity)
    qn        : (N,) nodal outflow from the Darcy solve (sol.nodal_flux), positive where
                fluid leaves the lobule. Required for the advective boundary term.
    tol       : convergence threshold on ||residual||, RELATIVE to the first residual.
                An absolute threshold is meaningless here: the residual is a molar flux
                over elements of ~1e-10 m^2, so its natural scale is ~1e-9 whatever the
                error is.
    """
    n = mesh.n
    g = (lambda c: c) if g is None else g
    dg = (lambda c: np.ones_like(c)) if dg is None else dg
    qn = np.zeros(n) if qn is None else np.asarray(qn, dtype=float)

    dirichlet = list(dirichlet) if dirichlet else []
    nodes = (np.concatenate([np.asarray(k, dtype=np.int64).ravel() for k, _ in dirichlet])
             if dirichlet else np.array([], dtype=np.int64))
    vals = (np.concatenate([np.broadcast_to(np.asarray(v, dtype=float),
                                            np.asarray(k).ravel().shape).ravel()
                            for k, v in dirichlet]) if dirichlet else np.array([]))

    c = np.zeros(n) if c0 is None else np.asarray(c0, dtype=float).copy()
    if nodes.size:
        c[nodes] = vals

    Kd = assemble_aniso(mesh, diff, diff) if diff > 0 else csr_matrix((n, n))
    b_e = np.einsum('ea,eia->ei', u, mesh.grad)          # (M, 3) = u . grad phi_i
    tau = supg_tau(mesh, u, diff) if supg else np.zeros(len(mesh.area))
    wsup = (tau * mesh.area)[:, None] * b_e              # SUPG test-function increment
    free = np.ones(n, dtype=bool)
    free[nodes] = False
    tri = mesh.tri

    def residual(cv):
        gc = g(cv)
        gm = gc[tri].mean(axis=1)
        R = react(cv)[0]
        res = Kd.dot(cv) + mesh.lumped * R + qn * gc                  # + boundary outflow
        np.add.at(res, tri, -(mesh.area * gm)[:, None] * b_e)        # -int grad w . u g
        if supg:                                                      # strong residual
            adv = np.einsum('ei,ei->e', b_e, gc[tri])
            np.add.at(res, tri, wsup * (adv + R[tri].mean(axis=1))[:, None])
        return res

    def jacobian(cv):
        dgc = dg(cv)
        R, dR = react(cv)
        Je = -(mesh.area[:, None, None] / 3.0) * b_e[:, :, None] * dgc[tri][:, None, :]
        if supg:
            Je = Je + wsup[:, :, None] * (b_e * dgc[tri])[:, None, :]
            Je = Je + wsup[:, :, None] * (dR[tri] / 3.0)[:, None, :]
        J = (scatter(mesh, Je) + Kd + diags(mesh.lumped * dR + qn * dgc)).tolil()
        for i in nodes:
            J.rows[i] = [i]
            J.data[i] = [1.0]
        return J.tocsr()

    res = residual(c)
    res[nodes] = c[nodes] - vals
    nrm0 = max(float(np.linalg.norm(res[free])), 1e-300)
    for it in range(max_iter):
        nrm = float(np.linalg.norm(res[free]))
        if verbose:
            print('  newton %2d  |res| %.3e  (rel %.2e)' % (it, nrm, nrm / nrm0))
        if nrm <= tol * nrm0:
            break
        dc = spsolve(jacobian(c), -res)
        if not np.all(np.isfinite(dc)):
            raise RuntimeError('Newton step is not finite (check the parameters)')
        step, ok = 1.0, False
        for _ in range(14):                                           # backtracking
            trial = c + step * dc
            rt = residual(trial)
            rt[nodes] = trial[nodes] - vals
            if float(np.linalg.norm(rt[free])) < nrm:
                c, res, ok = trial, rt, True
                break
            step *= 0.5
        if not ok:
            c = c + dc
            res = residual(c)
            res[nodes] = c[nodes] - vals
            break
    return c


