"""Solute transport across the lobule on the Darcy field.

Three problems, three reasons for being here:

oxygen
    The one that cannot be written as a diluted solute. Blood carries free O2 *and* the
    much larger haemoglobin-bound pool, so the advected quantity is
    C_tot(p) = alpha p + c_Hb S(p) with S a Hill curve, while diffusion acts on the free
    pool only. The capacity dC_tot/dp is tens of times alpha near p50 and falls at high
    pO2, which is what flattens the periportal part of the profile and steepens the
    pericentral part. A linear "diluted species" model misses this.

ethanol and acetaldehyde
    Two enzymes with Michaelis constants three orders of magnitude apart (ADH ~0.5 mM,
    CYP2E1 ~10 mM) and different zonal distributions. At a drinking concentration both are
    saturated, the kinetics are zero order, and ethanol falls almost linearly along the
    sinusoid. Acetaldehyde is a second species whose source is the ethanol sink.

drug
    A generic small molecule cleared by a pericentral CYP, i.e. the paracetamol geometry:
    the reactive-metabolite production rate maps the region at risk.

A zonation weight redistributes enzyme capacity along zeta (0 periportal, 1 pericentral)
without changing the lobule total, so sharpness can be swept at constant Vmax.

Every balance uses the discrete boundary flux returned by solve_adr, not Q * c_in: with a
Dirichlet inlet there is a diffusive inflow as well as an advective one, and ignoring it
inflates the extraction.
"""
import numpy as np
from . import fem_core as fc


def zonation(mesh, sharpness):
    """Nodal weight with area-weighted mean 1. sharpness > 0 is pericentral."""
    if abs(sharpness) < 1e-12:
        return np.ones(mesh.n)
    w = np.exp(sharpness * (mesh.zeta - 0.5))
    return w / (np.sum(w * mesh.lumped) / np.sum(mesh.lumped))


def surface_flux(mesh, ring, u, c, D, g=None, mask=None):
    """Outward solute flux through a boundary ring [mol/s per lobule].

    advective  (q . n) g(c)       with q the Darcy flux and g the advected content
    diffusive  -D (grad c . n)
    Both are exact for P1: grad c is constant on the element and q is constant too, so a
    single edge-midpoint value integrates the edge exactly.
    """
    g = (lambda x: x) if g is None else g
    e = ring['elem']
    gc = 0.5 * (g(c)[ring['a']] + g(c)[ring['b']])
    gradc = np.einsum('ei,eia->ea', c[mesh.tri], mesh.grad)[e]
    un = np.einsum('ea,ea->e', u[e], ring['normal'])
    dn = np.einsum('ea,ea->e', gradc, ring['normal'])
    f = (un * gc - D * dn) * ring['length'] * mesh.geom.H
    return float(np.sum(f[mask])) if mask is not None else float(np.sum(f))


class Result(object):
    """A solved field plus the integral quantities worth reporting.

    influx/outflux are per lobule [mol/s]; c_out is the mixed-cup outlet concentration.
    """

    def __init__(self, name, mesh, c, rate, c_in, Q, unit, sol, D, g=None, produced=0.0):
        H = mesh.geom.H
        self.name, self.mesh, self.c, self.rate, self.unit = name, mesh, c, rate, unit
        self.c_in, self.Q = c_in, Q
        self.consumed = float(np.sum(rate * mesh.lumped)) * H
        gg = (lambda x: x) if g is None else g
        qn = sol.nodal_flux
        self.outflux = float(np.sum(np.maximum(qn, 0.0) * gg(c))) * H
        self.produced = produced
        supply = self.outflux + self.consumed            # exact: the operator is conservative
        self.influx = supply - produced
        self.extraction = self.consumed / supply if supply > 0 else float('nan')
        # independent check: the same influx from the geometric boundary integral. It is
        # an O(h) evaluation of a quantity the operator gets exactly, so it is a mesh
        # diagnostic, not a conservation error.
        geo = -surface_flux(mesh, mesh.ring_outer, sol.u, c, D, g, mesh.portal_edges())
        scale = max(abs(self.influx), abs(produced), 1e-300)
        self.balance_error = abs(geo - self.influx) / scale
        self.c_out = self.outflux / Q if Q > 0 else float('nan')
        self.c_cv = float(np.mean(c[mesh.inner]))
        self.c_min, self.c_max = float(np.min(c)), float(np.max(c))

    def summary(self):
        return {self.name + '_in': self.c_in,
                self.name + '_out_mixed': self.c_out,
                self.name + '_at_central_vein': self.c_cv,
                self.name + '_min': self.c_min,
                self.name + '_extraction': self.extraction,
                self.name + '_turnover_mol_s': self.consumed,
                self.name + '_balance_error': self.balance_error,
                self.name + '_unit': self.unit}


# ------------------------------------------------------------------ oxygen
def oxygen(mesh, sol, o2):
    """c = pO2 [mmHg]; rate = O2 consumption [mol m^-3 s^-1]."""
    n, p50 = o2.hill, o2.p50
    w = zonation(mesh, -o2.zonation)        # o2.zonation > 0 puts consumption periportal

    def sat(p):
        pp = np.maximum(p, 0.0)
        return pp ** n / (p50 ** n + pp ** n)

    def g(p):
        return o2.alpha * p + o2.c_hb * sat(p)

    def dg(p):
        pp = np.maximum(p, 1e-6)
        den = p50 ** n + pp ** n
        return o2.alpha + o2.c_hb * n * pp ** (n - 1) * p50 ** n / den ** 2

    def react(p):
        pp = np.maximum(p, 0.0)
        R = o2.vmax * w * pp / (o2.km + pp)
        dR = np.where(p > 0, o2.vmax * w * o2.km / (o2.km + pp) ** 2, 0.0)
        return R, dR

    p = fc.solve_adr(mesh, sol.u, o2.D * o2.alpha, react,
                     [(mesh.portal, o2.p_in)], g=g, dg=dg,
                     c0=np.full(mesh.n, o2.p_in), qn=sol.nodal_flux)
    res = Result('pO2', mesh, p, react(p)[0], o2.p_in, sol.Q, 'mmHg', sol, o2.D * o2.alpha, g)
    res.c_out = _invert_capacity(res.outflux / sol.Q, g, 4.0 * o2.p_in)
    res.saturation = sat(p)
    res.sat_in = float(sat(np.array([o2.p_in]))[0])
    res.sat_out = float(sat(np.array([res.c_out]))[0])
    return res


def _invert_capacity(ctot, g, p_hi):
    """The pO2 that carries a given total O2 content, by bisection."""
    lo, hi = 0.0, float(p_hi)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if g(np.array([mid]))[0] < ctot:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------------ ethanol
def ethanol(mesh, sol, et):
    """ADH (uniform) plus CYP2E1 (pericentral). c in mol/m^3 = mM."""
    wc = zonation(mesh, et.cyp_zonation)

    def react(c):
        cc = np.maximum(c, 0.0)
        R = et.vmax_adh * cc / (et.km_adh + cc) + et.vmax_cyp * wc * cc / (et.km_cyp + cc)
        dR = np.where(c > 0,
                      et.vmax_adh * et.km_adh / (et.km_adh + cc) ** 2
                      + et.vmax_cyp * wc * et.km_cyp / (et.km_cyp + cc) ** 2, 0.0)
        return R, dR

    c = fc.solve_adr(mesh, sol.u, et.D, react, [(mesh.portal, et.c_in)],
                     c0=np.full(mesh.n, et.c_in), qn=sol.nodal_flux)
    rate = react(c)[0]
    res = Result('ethanol', mesh, c, rate, et.c_in, sol.Q, 'mol/m^3', sol, et.D)
    cc = np.maximum(c, 0.0)
    cyp = et.vmax_cyp * wc * cc / (et.km_cyp + cc)
    res.cyp_fraction = float(np.sum(cyp * mesh.lumped) / max(np.sum(rate * mesh.lumped), 1e-300))
    return res


def acetaldehyde(mesh, sol, etoh_result, ac):
    """Source = the ethanol oxidation rate; sink = ALDH2, nearly first order at these levels."""
    src = etoh_result.rate

    def react(c):
        cc = np.maximum(c, 0.0)
        R = ac.vmax_aldh * cc / (ac.km_aldh + cc) - src
        dR = np.where(c > 0, ac.vmax_aldh * ac.km_aldh / (ac.km_aldh + cc) ** 2, 0.0)
        return R, dR

    c = fc.solve_adr(mesh, sol.u, ac.D, react, [(mesh.portal, ac.c_in)],
                     c0=np.full(mesh.n, 1e-3), qn=sol.nodal_flux)
    cc = np.maximum(c, 0.0)
    clear = ac.vmax_aldh * cc / (ac.km_aldh + cc)
    produced = float(np.sum(src * mesh.lumped)) * mesh.geom.H
    res = Result('acetaldehyde', mesh, c, clear, ac.c_in, sol.Q, 'mol/m^3', sol, ac.D,
                 produced=produced)
    res.peak = float(np.max(c))
    res.peak_zeta = float(mesh.zeta[int(np.argmax(c))])
    return res


# ------------------------------------------------------------------ drug
def drug(mesh, sol, d):
    w = zonation(mesh, d.zonation)

    def react(c):
        cc = np.maximum(c, 0.0)
        R = d.vmax * w * cc / (d.km + cc)
        dR = np.where(c > 0, d.vmax * w * d.km / (d.km + cc) ** 2, 0.0)
        return R, dR

    c = fc.solve_adr(mesh, sol.u, d.D, react, [(mesh.portal, d.c_in)],
                     c0=np.full(mesh.n, d.c_in), qn=sol.nodal_flux)
    rate = react(c)[0]
    res = Result('drug', mesh, c, rate, d.c_in, sol.Q, 'mol/m^3', sol, d.D)
    res.reactive = d.toxic_yield * rate
    res.reactive_peak_zeta = float(mesh.zeta[int(np.argmax(res.reactive))])
    res.reactive_total = d.toxic_yield * res.consumed
    return res
