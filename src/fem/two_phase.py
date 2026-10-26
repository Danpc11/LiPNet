"""Sinusoid and hepatocyte as two coupled continua.

The single-continuum model in transport.py puts the solute and the enzyme in the same
place: whatever is in the blood is immediately available to be metabolised. For oxygen and
for ethanol that is close enough, because the sinusoidal endothelium is fenestrated and has
no basement membrane, so the barrier is a micron of water. For anything whose uptake is
transporter-limited it is wrong, and wrong in a way that matters: it makes the lobule look
like it clears the drug as fast as the enzyme can work, when in fact it clears it as fast
as the membrane can deliver it.

Two fields, vascular and tissue, with the metabolism moved to where the enzymes are:

    div(q g(c_v)) - div(phi D_v grad c_v)   + beta (c_v - c_t / K_p)            = 0
                  - div((1-phi) D_t grad c_t) - beta (c_v - c_t / K_p) + R(c_t) = 0

beta [1/s] is the exchange conductance per unit tissue volume. It is not a free parameter:

    beta = P a_v,     a_v = 4 phi / d_sin

with P the permeability of the endothelium plus the space of Disse [m/s] and a_v the
sinusoidal surface per unit volume, about 7.5e4 m^-1 for this lobule. For a solute crossing
by free diffusion through a micron of water, P ~ D/delta ~ 1e-3 m/s and beta ~ 1e2 1/s; for
a transporter-limited drug P can be 1e-7 m/s and beta ~ 1e-2 1/s.

The number that decides the regime is the exchange Damkohler number

    Da_x = beta * tau_transit

built on the sinusoidal transit time. Da_x >> 1 is the equilibrated limit, where
c_t -> K_p c_v and this model reduces exactly to transport.py: that reduction is checked in
tests/test_fem.py and is the reason to trust the block solver. Da_x << 1 is the delivery-
limited regime, where the apparent clearance is set by beta and not by V_max at all.

The vascular equation keeps the conservative advection and the nodal-flux boundary term of
fem_core, so influx = outflux + consumption remains an identity: the exchange term is
lumped at the nodes and appears with opposite signs in the two phases, so it cancels
exactly in the sum.
"""
import numpy as np
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import spsolve

from . import fem_core as fc
from .transport import surface_flux, zonation, _invert_capacity


def beta_from_permeability(permeability, flow):
    """beta = P a_v with a_v = 4 phi / d_sin, the sinusoidal surface per unit volume.
    Accepts a scalar or a nodal array."""
    a_v = 4.0 * flow.phi / flow.d_sin
    out = np.asarray(permeability, dtype=float) * a_v
    return float(out) if out.ndim == 0 else out


def zonal_permeability(mesh, p0, sharpness):
    """A membrane that is not the same everywhere.

    Sinusoidal fenestration is not uniform along the lobule: periportal sinusoids have more
    and smaller fenestrae, pericentral ones fewer and larger, and the reported net porosity
    is higher in zone 3. sharpness > 0 puts the open membrane pericentrally, < 0
    periportally, 0 is uniform. The lobule mean is held at p0, so a sweep in sharpness
    redistributes the exchange capacity without changing how much of it there is: any
    effect is then the zonation, not a net gain or loss of surface.
    """
    if abs(sharpness) < 1e-12:
        return np.full(mesh.n, float(p0))
    w = np.exp(sharpness * (mesh.zeta - 0.5))
    w = w / (np.sum(w * mesh.lumped) / np.sum(mesh.lumped))
    return float(p0) * w


class TwoPhaseResult(object):
    def __init__(self, name, mesh, c_v, c_t, rate, c_in, sol, D_v, unit, beta, kp,
                 g=None, capacity=1.0):
        H = mesh.geom.H
        self.name, self.mesh, self.unit = name, mesh, unit
        self.c_v, self.c_t, self.rate = c_v, c_t, rate
        self.c_in, self.kp = c_in, kp
        gg = (lambda x: x) if g is None else g
        qn = sol.nodal_flux
        self.consumed = float(np.sum(rate * mesh.lumped)) * H
        self.outflux = float(np.sum(np.maximum(qn, 0.0) * gg(c_v))) * H
        self.influx = self.outflux + self.consumed
        self.extraction = self.consumed / self.influx if self.influx > 0 else float('nan')
        # outflux is a flux of CONTENT, so for oxygen it has to be turned back into a
        # partial pressure through the same capacity curve
        content_out = self.outflux / sol.Q
        self.c_out = (content_out if g is None
                      else _invert_capacity(content_out, gg, 4.0 * max(c_in, 1.0)))
        # relaxation of the tissue phase towards the blood takes capacity / beta, where
        # capacity = dC/dc on the tissue side (1 for a free solute, the solubility alpha
        # when the unknown is a partial pressure). Comparing that with the transit time is
        # what decides whether the membrane matters.
        self.capacity = capacity
        self.beta = float(np.sum(np.asarray(beta) * mesh.lumped) / np.sum(mesh.lumped)) \
            if np.ndim(beta) else float(beta)
        self.beta_field = np.asarray(beta)
        self.Da_exchange = self.beta / capacity * sol.transit_time()
        w = mesh.lumped
        self.gradient = float(np.sum((c_v - c_t / kp) * w) / np.sum(w))
        self.c_v_mean = float(np.sum(c_v * w) / np.sum(w))
        self.c_t_mean = float(np.sum(c_t * w) / np.sum(w))
        geo = -surface_flux(mesh, mesh.ring_outer, sol.u, c_v, D_v, g, mesh.portal_edges())
        self.balance_error = abs(geo - self.influx) / max(abs(self.influx), 1e-300)

    def summary(self):
        return {self.name + '_in': self.c_in,
                self.name + '_out_mixed': self.c_out,
                self.name + '_sinusoid_mean': self.c_v_mean,
                self.name + '_hepatocyte_mean': self.c_t_mean,
                self.name + '_mean_gradient': self.gradient,
                self.name + '_extraction': self.extraction,
                self.name + '_beta_per_s': self.beta,
                self.name + '_Da_exchange': self.Da_exchange,
                self.name + '_balance_error': self.balance_error,
                self.name + '_unit': self.unit}


def solve(mesh, sol, react, c_in, D_v, D_t, beta, kp=1.0, g=None, dg=None,
          name='solute', unit='mol/m^3', capacity=1.0, supg=True, tol=1e-9, max_iter=50):
    """Block Newton on the two coupled fields. react acts on the TISSUE concentration."""
    n = mesh.n
    g = (lambda c: c) if g is None else g
    dg = (lambda c: np.ones_like(c)) if dg is None else dg
    phi = sol.flow.phi
    qn = sol.nodal_flux
    u = sol.u
    beta = np.broadcast_to(np.asarray(beta, dtype=float), (n,)).copy()

    Kv = fc.assemble_aniso(mesh, phi * D_v, phi * D_v) if D_v > 0 else csr_matrix((n, n))
    Kt = fc.assemble_aniso(mesh, (1 - phi) * D_t, (1 - phi) * D_t) if D_t > 0 else csr_matrix((n, n))
    b_e = np.einsum('ea,eia->ei', u, mesh.grad)
    tau = fc.supg_tau(mesh, u, phi * D_v) if supg else np.zeros(len(mesh.area))
    wsup = (tau * mesh.area)[:, None] * b_e
    tri, lump = mesh.tri, mesh.lumped

    nodes = np.asarray(mesh.portal, dtype=np.int64)
    free = np.ones(n, dtype=bool)
    free[nodes] = False

    c_v = np.full(n, float(c_in))
    c_t = np.full(n, float(c_in) * kp)

    def pieces(cv, ct):
        gc = g(cv)
        exch = beta * (cv - ct / kp)
        R = react(ct)[0]
        rv = Kv.dot(cv) + qn * gc + lump * exch
        np.add.at(rv, tri, -(mesh.area * gc[tri].mean(axis=1))[:, None] * b_e)
        if supg:
            adv = np.einsum('ei,ei->e', b_e, gc[tri])
            np.add.at(rv, tri, wsup * (adv + exch[tri].mean(axis=1))[:, None])
        rt = Kt.dot(ct) - lump * exch + lump * R
        return rv, rt

    def residual(cv, ct):
        rv, rt = pieces(cv, ct)
        rv[nodes] = cv[nodes] - c_in
        return np.concatenate([rv, rt])

    def jacobian(cv, ct):
        dgc = dg(cv)
        dR = react(ct)[1]
        Je = -(mesh.area[:, None, None] / 3.0) * b_e[:, :, None] * dgc[tri][:, None, :]
        if supg:
            Je = Je + wsup[:, :, None] * (b_e * dgc[tri])[:, None, :]
            Je = Je + wsup[:, :, None] * (beta[tri] / 3.0)[:, None, :]
        Jvv = (fc.scatter(mesh, Je) + Kv + diags(qn * dgc + lump * beta)).tolil()
        Jvt = diags(-lump * beta / kp).tolil()
        if supg:
            Jvt = (Jvt + fc.scatter(mesh, -(wsup[:, :, None]
                                            * (beta[tri] / kp / 3.0)[:, None, :]))).tolil()
        for i in nodes:
            Jvv.rows[i] = [i]
            Jvv.data[i] = [1.0]
            Jvt.rows[i] = []
            Jvt.data[i] = []
        Jtv = diags(-lump * beta)
        Jtt = Kt + diags(lump * (beta / kp + dR))
        return bmat([[Jvv.tocsr(), Jvt.tocsr()], [Jtv, Jtt]], format='csr')

    res = residual(c_v, c_t)
    nrm0 = max(float(np.linalg.norm(res)), 1e-300)
    for _ in range(max_iter):
        nrm = float(np.linalg.norm(res))
        if nrm <= tol * nrm0:
            break
        d = spsolve(jacobian(c_v, c_t), -res)
        if not np.all(np.isfinite(d)):
            raise RuntimeError('two-phase Newton step is not finite')
        step, ok = 1.0, False
        for _ in range(14):
            tv, tt = c_v + step * d[:n], c_t + step * d[n:]
            rt_ = residual(tv, tt)
            if float(np.linalg.norm(rt_)) < nrm:
                c_v, c_t, res, ok = tv, tt, rt_, True
                break
            step *= 0.5
        if not ok:
            c_v, c_t = c_v + d[:n], c_t + d[n:]
            res = residual(c_v, c_t)
            break
    return TwoPhaseResult(name, mesh, c_v, c_t, react(c_t)[0], c_in, sol,
                          phi * D_v, unit, beta, kp, g, capacity)


# ------------------------------------------------------------------ species wrappers
def drug(mesh, sol, d, permeability, kp=1.0, name='drug'):
    """The same drug as transport.drug, but taken up across a membrane of finite
    permeability. Metabolism now sees the hepatocyte concentration, not the blood one."""
    w = zonation(mesh, d.zonation)

    def react(c):
        cc = np.maximum(c, 0.0)
        R = d.vmax * w * cc / (d.km + cc)
        dR = np.where(c > 0, d.vmax * w * d.km / (d.km + cc) ** 2, 0.0)
        return R, dR

    beta = beta_from_permeability(permeability, sol.flow)
    res = solve(mesh, sol, react, d.c_in, d.D, d.D, beta, kp, name=name)
    res.reactive = d.toxic_yield * res.rate
    res.reactive_peak_zeta = float(mesh.zeta[int(np.argmax(res.reactive))])
    res.reactive_total = d.toxic_yield * res.consumed
    res.permeability = permeability
    return res


def oxygen(mesh, sol, o2, permeability=2.0e-3):
    """Oxygen with an explicit, if nearly irrelevant, membrane. Da_x is of order 1e3, so
    this must reproduce the single-continuum answer; it is the reduction test."""
    n, p50 = o2.hill, o2.p50
    w = zonation(mesh, -o2.zonation)

    def sat(p):
        pp = np.maximum(p, 0.0)
        return pp ** n / (p50 ** n + pp ** n)

    def gg(p):
        return o2.alpha * p + o2.c_hb * sat(p)

    def dgg(p):
        pp = np.maximum(p, 1e-6)
        den = p50 ** n + pp ** n
        return o2.alpha + o2.c_hb * n * pp ** (n - 1) * p50 ** n / den ** 2

    def react(p):
        pp = np.maximum(p, 0.0)
        return (o2.vmax * w * pp / (o2.km + pp),
                np.where(p > 0, o2.vmax * w * o2.km / (o2.km + pp) ** 2, 0.0))

    # in pO2 units the exchange flux is P a_v alpha (p_v - p_t), and the tissue capacity
    # is alpha, so the relaxation time is 1 / (P a_v) as it should be
    beta = beta_from_permeability(permeability, sol.flow) * o2.alpha
    return solve(mesh, sol, react, o2.p_in, o2.D * o2.alpha, o2.D * o2.alpha,
                 beta, 1.0, g=gg, dg=dgg, name='pO2', unit='mmHg', capacity=o2.alpha)


def permeability_sweep(mesh, sol, d, permeabilities, kp=1.0):
    """Extraction and zonal exposure against membrane permeability. The point of the whole
    module: where the curve bends is where V_max stops being the answer."""
    out = dict(permeability_m_s=[], beta_per_s=[], Da_exchange=[], extraction=[],
               hepatocyte_over_blood=[], reactive_peak_zeta=[])
    for P in permeabilities:
        r = drug(mesh, sol, d, P, kp)
        out['permeability_m_s'].append(float(P))
        out['beta_per_s'].append(r.beta)
        out['Da_exchange'].append(r.Da_exchange)
        out['extraction'].append(r.extraction)
        out['hepatocyte_over_blood'].append(r.c_t_mean / max(r.c_v_mean, 1e-300))
        out['reactive_peak_zeta'].append(r.reactive_peak_zeta)
    return out


def zonation_matrix(mesh, sol, d, permeability, sharpness=(-3.0, 0.0, 3.0)):
    """Extraction for every combination of enzyme zonation and membrane zonation, at a
    fixed lobule-mean permeability. The diagonal is the point: a membrane matched to the
    enzyme clears more than a mismatched one, and the penalty for mismatch is of order
    10%, an order of magnitude smaller than the effect of the permeability itself."""
    import copy
    out = dict(enzyme_sharpness=[float(x) for x in sharpness],
               membrane_sharpness=[float(x) for x in sharpness], extraction=[])
    for ez in sharpness:
        dd = copy.deepcopy(d)
        dd.zonation = float(ez)
        row = [drug(mesh, sol, dd, zonal_permeability(mesh, permeability, mz)).extraction
               for mz in sharpness]
        out['extraction'].append([float(v) for v in row])
    return out
