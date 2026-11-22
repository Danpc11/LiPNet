"""The network model in physical units for a human liver, its growth with the body and the regrowth of a graft.

src/model/liver.py works in the units of the optimisation. This module puts the same optimum in millimetres,
mL/min, mmHg and Pa, so that it can be set against what has been measured in people, and it makes the three layers of
the repository describe one liver.

1. Calibration. Four numbers fix the canopy of a human liver: its volume, the portal flow per 100 g (90 mL/min, the
   manuscript's donor reference), the portocaval gradient (5 mmHg) and the viscosity (3.5 mPa s, as in src/fem).
   Two more come from the human cast of Debbaut et al. (J Anat 2014;224:509, one liver of 1.9 kg):

   BRANCH_LENGTH_FACTOR  a branch is 0.35 to 0.50 times the side of the territory it supplies (portal generations
                         7 to 11), not one side. Every length of the canopy is multiplied by 0.44. With it the portal
                         trunk is 54 mm long (49 mm in adults, pooled; Chau et al. 2023), without it 122 mm.
   LOBULE_DROP_MMHG      the drop from portal triad to central vein at donor flow, imposed on the sinusoidal edge
                         instead of optimised, and shared with src/fem/params.py. Left to the optimum the lobule
                         would hold 0.25 to 0.4 mmHg. See `liver`.

   With these and nothing fitted, `portal_radii` compares the absolute radii of the portal tree with the eleven
   measured generations of the cast, and `cast_exponent` reads b from the cast's own radii, lengths and counts.

2. The hepatic veins of the cast are not used to fix b. `cast_poiseuille` shows why: their seven measured generations
   hold 3.5 mmHg at donor flow and the drop grows 1.5-fold per generation, which no liver with a 5 mmHg gradient can
   have. The authors of the cast met the same problem: extrapolated as measured, their model put the sinusoidal
   pressure above the portal pressure (Debbaut et al. 2011).

3. Growth with the body. `growth` scales the calibrated liver between a newborn and an adult at invariant gradient
   and invariant flow per gram, with one free number the human data do not yet fix: the share s of the growth in
   volume that is enlargement of the lobules (s = 1) rather than addition of lobules (s = 0).

4. Regrowth of a graft. After transplantation a graft grows by enlarging its lobules, not by adding them (rat:
   Papp 2013; Tsomaia et al. 2020). `regeneration` gives the loads of a graft whose volume has grown rho-fold, for a
   given change of inflow, between two bounds for the calibre of its vessels: unchanged, or grown in proportion.
   `first_week` sets both against the portal pressure measured in recipients during the first week.

5. Two inflows, measured. `two_inflows_measured` evaluates the two-inflow relation of liver.py in the two series that
   measured portal and arterial flow in the same tissue before and after transplantation.

    python src/model/human.py      -> results/model_human.json
"""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import indices
from model import liver as net
from fem.params import LOBULE_DROP_MMHG, MMHG, Flow

VISCOSITY = Flow().mu                        # Pa s, the same as the lobule layer
PORTAL_PER_100G = indices.DEFAULT_DONOR_REF  # mL/min per 100 g, the donor reference of the loads
GRADIENT_MMHG = indices.NORMAL_GRADIENT      # portal vein to vena cava, the reference of the pressure load
DENSITY = 1.05                               # g/mL
ADULT_LIVER_G, CAST_LIVER_G, NEWBORN_LIVER_G = 1500.0, 1900.0, 150.0
BRANCH_LENGTH_FACTOR = 0.44
B_HUMAN = 0.75                               # the portal tree of the cast: 0.73 from its radius law, 0.74 from its absolute radii (0.68 to 0.80)
B_CALIBRATED = (0.6, 2 / 3, 0.75, 0.85, 0.93)  # the calibrated human livers of the rigid range (compartments.geometries)

# One human liver cast, generation -> (mean radius mm, mean length mm, number of vessels).
# Generations 1 to 5 or 6 were measured on the whole cast (Debbaut et al., IEEE Trans Biomed Eng 2011;58:25, Tables II-IV);
# the following ones on dissected samples, with estimated numbers of vessels (Debbaut et al., J Anat 2014;224:509, Table 1).
CAST_MACRO = {
    'portal vein': {1: (7.3, 35, 1), 2: (6.3, 35, 2), 3: (4.0, 33, 7), 4: (2.1, 32, 27), 5: (1.5, 27, 80), 6: (1.1, 16, 220)},
    'hepatic vein': {1: (13, 68, 1), 2: (5.0, 74, 3), 3: (2.4, 33, 18), 4: (1.7, 28, 33), 5: (1.1, 19, 76)},        # generation 1 is the vena cava
    'hepatic artery': {1: (3.5, 10, 1), 2: (2.9, 74, 2), 3: (1.8, 40, 7), 4: (1.1, 31, 36), 5: (0.83, 25, 120), 6: (0.59, 16, 270)}}
CAST_MESO = {
    'portal vein': {7: (0.718, 8.14, 432), 8: (0.482, 5.24, 1339), 9: (0.309, 3.16, 3715), 10: (0.231, 2.68, 10541), 11: (0.175, 1.48, 23976)},
    'hepatic vein': {6: (0.741, 8.56, 171), 7: (0.456, 6.59, 494), 8: (0.246, 4.09, 2109)},
    'hepatic artery': {7: (0.345, 7.22, 542), 8: (0.277, 3.51, 1536), 9: (0.241, 3.09, 3433), 10: (0.208, 2.36, 8943), 11: (0.160, 1.48, 15447)}}
CAST = {t: {**{g: v for g, v in CAST_MACRO[t].items() if not (t == 'hepatic vein' and g == 1)}, **CAST_MESO[t]} for t in CAST_MESO}
CAST_GENERATIONS = {'portal vein': (17, 18), 'hepatic vein': (13, 14)}       # estimated totals (J Anat 2014); the 2011 model used 16 and 14

# What has been measured in people, for `growth` and `regeneration` (sources in the comments below)
PORTAL_DIAMETER_MM = dict(children=(6.60, 5.38, 7.82), adults=(10.72, 10.25, 11.19))             # pooled mean and 95% CI, Chau et al. 2023
PORTAL_BY_WEIGHT = [((1.56, 5.0), (3.0, 4.0), (4.0, 12.0)), ((5.0, 15.0), (4.0, 7.0), (6.0, 12.0)),   # weight kg, calibre mm, mean velocity cm/s
                    ((15.0, 24.0), (4.0, 8.0), (9.0, 17.0)), ((24.0, 40.0), (6.0, 9.0), (10.0, 19.0))]   # ranges, Lopez Barrio et al. 1996
GRAFT_GROWTH = {0.5: 1.58, 1.0: 1.50, 3.0: 1.44, 6.0: 1.31}                                       # months -> volume / implanted volume, Zhang et al. 2022
# Volume of the graft over its implanted volume in adult recipients: (series, n, days after transplantation, ratio)
GRAFT_VOLUME = [('Yagi 2005, PVP < 20', 15, 14, 1.02 / 0.633), ('Yagi 2005, PVP < 20', 15, 28, 0.975 / 0.633), ('Yagi 2005, PVP < 20', 15, 90, 0.943 / 0.633),
                ('Yagi 2005, PVP < 20', 15, 180, 0.963 / 0.633), ('Yagi 2005, PVP >= 20', 17, 14, 1.24 / 0.638), ('Yagi 2005, PVP >= 20', 17, 28, 1.27 / 0.638),
                ('Yagi 2005, PVP >= 20', 17, 90, 1.29 / 0.638), ('Yagi 2005, PVP >= 20', 17, 180, 0.940 / 0.638), ('Eguchi 2003, cirrhosis', 9, 14, 1.0 / 0.601),
                ('Byun 2016', 32, 14, 1.7), ('Zhang 2022', 62, 15, 1.58), ('Zhang 2022', 62, 30, 1.50), ('Zhang 2022', 62, 90, 1.44), ('Zhang 2022', 62, 180, 1.31),
                ('Jiang 2009', 18, 30, 2.194), ('Garcia-Valdecasas 2003', 22, 60, 1 / (1 - 0.449)), ('Olthoff 2015', 308, 90, 1 + 549 / 989), ('Wu 2011', 64, 90, 2.093)]
# Both inflows measured in the same tissue, in the donor before retrieval and in the graft after reperfusion (transit-time flowmeter)
TWO_INFLOWS_MEASURED = {
    'Wu 2011': dict(n=64, donor_portal=84.66, donor_arterial=84.66 * 237.89 / 527.81, graft_portal=275.53, graft_arterial=275.53 * 157.69 / 1713.98,   # mL/min per 100 g, means
                    gradient_mmHg=10.6, graft_g=629.7),
    'Garcia-Valdecasas 2003': dict(n=22, donor_portal=69.6, donor_arterial=20.2, graft_portal=243.7, graft_arterial=16.0, gradient_mmHg=None, graft_g=759.0)}   # medians
# The first week in the recipient
PORTAL_INFLOW_DAY10 = 1198.5 / 1700.0        # portal flow at day 10 over the day of transplantation, 64 recipients (Wu 2011)
VOLUME_DAY7 = 1.0 / 0.636                    # graft volume at one week over the implanted volume, 32 recipients (Yagi 2005: GV/SLV 0.64 to about 1.0)
PORTAL_PRESSURE_DAY7 = dict(high=(24.2, 22.0), low=(16.2, 14.8))   # mean portal pressure at 12 h and day 7, mmHg, read from Fig. 1 of Yagi 2005 (17 and 15 recipients)
DONOR_PORTAL_PER_100G = {'Troisi 2003': 91.0, 'Troisi 2005': 118.0, 'Chan 2011': 81.0, 'Wu 2011': 84.66, 'Garcia-Valdecasas 2003': 69.6}

ml_min = lambda q: q * 1e-6 / 60             # mL/min -> m3/s


def cast_arrays(tree, part='all'):
    """Generations, radii (mm), lengths (mm) and numbers of vessels of the cast: 'all' measured generations, the
    'macro' ones (whole cast) or the 'meso' ones (dissected samples)."""
    t = dict(all=CAST, macro={k: {g: v for g, v in d.items() if g in CAST[k]} for k, d in CAST_MACRO.items()}, meso=CAST_MESO)[part][tree]; g = sorted(t)
    return np.array(g), np.array([t[k][0] for k in g], float), np.array([t[k][1] for k in g], float), np.array([t[k][2] for k in g], float)


def cast_exponent(tree, part='all'):
    """b from the radius law of the optimum, log r = c - log N / (2 + b) + (1 - b) log l / (2 (2 + b)), with the flow
    of a generation shared equally among its N vessels. Returns the fit and its range leaving one generation out."""
    g, r, l, N = cast_arrays(tree, part)

    def fit(r, l, N):
        res = lambda p: np.log(r) - (p[0] - np.log(N) / (2 + p[1]) + (1 - p[1]) / (2 * (2 + p[1])) * np.log(l))
        return float(least_squares(res, [0.0, 0.9]).x[1])
    loo = [fit(np.delete(r, k), np.delete(l, k), np.delete(N, k)) for k in range(len(g))] if len(g) > 3 else []
    side = (CAST_LIVER_G / DENSITY * 1e3 / N) ** (1 / 3)
    return dict(tree=tree, part=part, generations=[int(g[0]), int(g[-1])], b=fit(r, l, N), leave_one_out=[float(min(loo)), float(max(loo))] if loo else None,
                murray_exponent=float(-1 / np.polyfit(np.log(N), np.log(r), 1)[0]), length_over_territory_side=[float(x) for x in l / side])


def cast_poiseuille(tree, radius_scale=1.0, liver_g=CAST_LIVER_G, part='all'):
    """Pressure drop (mmHg) and wall shear (Pa) of each measured generation at donor flow, Poiseuille."""
    portal = PORTAL_PER_100G * liver_g / 100; total = portal / (1 - net.ARTERIAL_FRACTION)
    Q = dict({'portal vein': portal, 'hepatic vein': total, 'hepatic artery': total - portal})[tree]
    g, r, l, N = cast_arrays(tree, part); r = r * radius_scale * 1e-3; l = l * 1e-3; f = ml_min(Q) / N
    dP = 8 * VISCOSITY * l * f / (np.pi * r ** 4) / MMHG; tau = 4 * VISCOSITY * f / (np.pi * r ** 3); k = len(g) - 1
    return dict(tree=tree, flow_mL_min=float(Q), generations=[int(x) for x in g], drop_mmHg=[float(x) for x in dP], shear_Pa=[float(x) for x in tau],
                drop_measured_generations=float(dP.sum()), drop_ratio_per_generation=float((dP[-1] / dP[0]) ** (1 / k)),
                shear_ratio_per_generation=float((tau[-1] / tau[0]) ** (1 / k)))


def cast_tree_drop(tree, generations, radius_scale=1.0):
    """Poiseuille drop (mmHg) of a whole tree of the cast at donor flow: the measured generations as they are and the
    unmeasured ones below them along the trend of the dissected samples (radius, length and number per generation)."""
    p = cast_poiseuille(tree, radius_scale); g, r, l, N = cast_arrays(tree, 'meso'); k = len(g) - 1
    q = (l[-1] / l[0]) ** (1 / k) / (N[-1] / N[0]) ** (1 / k) / (r[-1] / r[0]) ** (4 / k)        # drop of a generation over the one above it
    return float(p['drop_measured_generations'] + sum(p['drop_mmHg'][-1] * q ** j for j in range(1, generations - p['generations'][-1] + 1)))


def liver(b=B_HUMAN, g=net.HUMAN_GENERATIONS, liver_g=ADULT_LIVER_G, portal_per_100g=PORTAL_PER_100G, gradient=GRADIENT_MMHG,
          lobule_drop=LOBULE_DROP_MMHG, kappa=BRANCH_LENGTH_FACTOR, n=3.0, d=3.0, lobule_conductance=None, arterial_fraction=net.ARTERIAL_FRACTION):
    """The canopy of src/model/liver.py in physical units.

    Portal flow in the portal tree, portal plus arterial flow in the sinusoids and the hepatic veins; Poiseuille
    conductance pi r^4 / (8 mu l); the two trees at the optimum for the gradient left by the lobule.
    lobule_drop = None leaves the sinusoidal edge to the optimum as well (the model of liver.py);
    lobule_conductance (m3 s-1 Pa-1 per lobule) overrides lobule_drop, for a lobule of another size.
    """
    N = n ** (g - 1); V = liver_g / DENSITY * 1e-6; side = (V / N) ** (1 / 3); Lh = kappa * side
    portal = ml_min(portal_per_100g * liver_g / 100); total = portal / (1 - arterial_fraction)
    alpha = b / 2; i = np.arange(g); c1 = n ** i; l1 = Lh * n ** ((g - 1 - i) / d)
    cnt = np.concatenate([c1, [N], c1[::-1]]); l = np.concatenate([l1, [Lh], l1[::-1]]); f = np.concatenate([portal / c1, [total / N], total / c1[::-1]])
    u = (f ** 2 / l ** (1.5 * b)) ** (1 / (1 + alpha)); k = 8 * VISCOSITY / np.pi; target = gradient * MMHG
    if lobule_drop is None and lobule_conductance is None:
        w = u * k * np.sum(f / u) / target
    else:
        G = lobule_conductance if lobule_conductance is not None else f[g] / (lobule_drop * MMHG)
        rest = target - f[g] / G; tree = np.arange(len(f)) != g
        if rest <= 0: raise ValueError('the lobule alone exceeds the gradient')
        w = u * k * np.sum(f[tree] / u[tree]) / rest; w[g] = G * k
    r = (w * l) ** 0.25; dP = k * f / w; tau = 4 * VISCOSITY * f / (np.pi * r ** 3); vol = cnt * np.pi * r ** 2 * l
    return dict(b=float(b), g=int(g), n=float(n), N=float(N), liver_g=float(liver_g), territory_side_mm=float(side * 1e3), cnt=cnt, f=f, r_mm=r * 1e3, l_mm=l * 1e3,
                dP_mmHg=dP / MMHG, tau_Pa=tau, share_inlet=float(dP[:g].sum() / target), share_sinusoid=float(dP[g] / target), share_outlet=float(dP[g + 1:].sum() / target),
                trunks_removed=float((dP[0] + dP[-1]) / target), lobule_drop_mmHg=float(dP[g] / MMHG), lobule_conductance=float(f[g] / (dP[g])),
                portal_trunk_radius_mm=float(r[0] * 1e3), portal_trunk_length_mm=float(l[0] * 1e3), portal_trunk_velocity_cm_s=float(f[0] / (np.pi * r[0] ** 2) * 100),
                portal_trunk_shear_Pa=float(tau[0]), portal_terminal_shear_Pa=float(tau[g - 1]), hepatic_trunk_radius_mm=float(r[-1] * 1e3),
                tree_volume_mL=float((vol[:g].sum() + vol[g + 1:].sum()) * 1e6), tree_volume_fraction=float((vol[:g].sum() + vol[g + 1:].sum()) / V))


def radius_at(c, counts, side='portal'):
    """Model radius (mm) of the generation holding `counts` vessels, interpolated on logarithms."""
    g = c['g']; cnt, r = (c['cnt'][:g], c['r_mm'][:g]) if side == 'portal' else (c['cnt'][g + 1:][::-1], c['r_mm'][g + 1:][::-1])
    return np.exp(np.interp(np.log(counts), np.log(cnt), np.log(r)))


def portal_radii(b=B_HUMAN, **kw):
    """Absolute portal radii of the model against the cast, the eleven measured generations, with no fitted scale."""
    g, cast, l, counts = cast_arrays('portal vein')
    kw.setdefault('liver_g', CAST_LIVER_G); model = radius_at(liver(b, **kw), counts); q = model / cast
    return dict(b=float(b), generations=[int(x) for x in g], vessels=[float(x) for x in counts], cast_mm=[float(x) for x in cast], model_mm=[float(x) for x in model], model_over_cast=[float(x) for x in q],
                rms_log=float(np.sqrt(np.mean(np.log(q) ** 2))))


def portal_exponent(bs=np.round(np.arange(0.50, 1.001, 0.01), 2), tolerance=0.10, **kw):
    """The b that reproduces the absolute portal radii of the cast, and the interval within `tolerance` rms."""
    e = np.array([portal_radii(float(b), **kw)['rms_log'] for b in bs]); ok = bs[e < tolerance]
    return dict(b=float(bs[e.argmin()]), rms_log=float(e.min()), interval=[float(ok.min()), float(ok.max())] if len(ok) else None,
                rms_at_half=float(e[0]), rms_at_one=float(e[-1]))


def growth(b=B_HUMAN, s=1.0, sizes=None, **kw):
    """The calibrated liver from a newborn to an adult at invariant gradient and invariant flow per gram.

    s is the share of the growth in volume that is lobule enlargement: lobule volume ~ V^s, number of lobules
    ~ V^(1-s). The number of generations is kept and the branching number carries the change in N. The sinusoidal
    bed keeps its permeability, so the conductance of a lobule grows with its length (Darcy flow across an annulus
    of fixed proportions): the lobular drop at a fixed flow per gram grows as V^(2s/3).
    Returns the liver at each size and the exponents with liver volume between the two ends.
    """
    sizes = np.geomspace(NEWBORN_LIVER_G, ADULT_LIVER_G, 9) if sizes is None else np.asarray(sizes, float)
    adult = liver(b, liver_g=ADULT_LIVER_G, **kw); g = adult['g']; rows = []
    for m in sizes:
        x = m / ADULT_LIVER_G; N = adult['N'] * x ** (1 - s); n = N ** (1 / (g - 1))
        c = liver(b, liver_g=m, n=n, lobule_conductance=adult['lobule_conductance'] * x ** (s / 3), **kw)
        rows.append(dict(liver_g=float(m), lobules=float(N), lobule_side_mm=c['territory_side_mm'], portal_diameter_mm=2 * c['portal_trunk_radius_mm'],
                         portal_velocity_cm_s=c['portal_trunk_velocity_cm_s'], portal_trunk_shear_Pa=c['portal_trunk_shear_Pa'], lobule_drop_mmHg=c['lobule_drop_mmHg'],
                         sinusoidal_velocity=float((c['f'][g] / c['territory_side_mm'] ** 2) / (adult['f'][g] / adult['territory_side_mm'] ** 2)),
                         tree_volume_fraction=c['tree_volume_fraction']))
    ex = lambda k: float(np.log(rows[-1][k] / rows[0][k]) / np.log(sizes[-1] / sizes[0]))
    return dict(b=float(b), s=float(s), rows=rows, exponents={k: ex(k) for k in ('portal_diameter_mm', 'portal_velocity_cm_s', 'portal_trunk_shear_Pa', 'lobule_drop_mmHg',
                                                                                'sinusoidal_velocity', 'tree_volume_fraction')})


def urata_volume(weight_kg):
    """Standard liver volume (mL) = 706.2 BSA + 2.4 (Urata et al. 1995), with BSA from weight alone (Costeff)."""
    return 706.2 * (4 * weight_kg + 7) / (weight_kg + 90) + 2.4


def measured_growth():
    """Exponents of portal calibre and mean velocity with liver volume from the weight groups of Lopez Barrio et al.
    1996 (ranges only): central value from the mid-ranges, interval from the opposite ends of the first and last group."""
    V = [urata_volume(np.mean(w)) for w, _, _ in PORTAL_BY_WEIGHT]; dV = np.log(V[-1] / V[0]); out = dict(liver_volume_mL=[float(x) for x in V])
    for name, j in (('portal_diameter', 1), ('portal_velocity', 2)):
        lo0, hi0 = PORTAL_BY_WEIGHT[0][j]; lo1, hi1 = PORTAL_BY_WEIGHT[-1][j]
        out[name] = dict(exponent=float(np.log((lo1 + hi1) / (lo0 + hi0)) / dV), interval=[float(np.log(lo1 / hi0) / dV), float(np.log(hi1 / lo0) / dV)])
    return out


def regeneration(rho, zF0=1.0, zR0=1.0, b=B_HUMAN, calibre='unchanged', inflow=1.0, **kw):
    """Loads of a graft whose volume has grown rho-fold by enlarging its lobules, with `inflow` times the inflow of
    reperfusion.

    Lengths grow as rho^(1/3). The lobule keeps its permeability, so its resistance falls as rho^(-1/3). The
    resistance of the trees grows as rho^(1/3) if their calibres are unchanged and falls as 1/rho if the calibres
    grow in proportion to the lengths. zF0 and zR0 are the loads at reperfusion; loads are per unit of current mass.
    """
    c = liver(b, **kw); rho = np.asarray(rho, float); s_lob = c['share_sinusoid']
    tree = rho ** (1 / 3) if calibre == 'unchanged' else 1 / rho
    resistance = (1 - s_lob) * tree + s_lob * rho ** (-1 / 3)                 # of the whole graft, relative to reperfusion
    zF = zF0 * inflow / rho; zP = zF0 * zR0 * inflow * resistance
    return dict(rho=rho, zF=zF, zP=zP, zR=zP / zF, resistance=resistance, sinusoidal_velocity=zF0 * inflow * rho ** (-2 / 3), calibre=calibre, lobule_share=float(s_lob))


def two_inflows_measured(bs=B_CALIBRATED):
    """The two-inflow relation of liver.py with the arterial share and the arterial load measured, not assumed.

    a is the arterial-to-portal flow ratio of the donor tissue and phi the arterial flow per gram of the graft over
    that of the donor. Returns the measured loads and the resistance load the model gives for a lobar graft."""
    out = {}
    for name, d in TWO_INFLOWS_MEASURED.items():
        zF = d['graft_portal'] / d['donor_portal']; a = d['donor_arterial'] / d['donor_portal']; phi = d['graft_arterial'] / d['donor_arterial']
        total = (d['graft_portal'] + d['graft_arterial']) / (d['donor_portal'] + d['donor_arterial'])
        model = [net.two_inflows(zF, c['share_inlet'], 1 - c['trunks_removed'], arterial_fraction=a / (1 + a), arterial_load=phi)['zR'] for c in (liver(b) for b in bs)]
        r = dict(n=d['n'], flow_load=zF, arterial_to_portal_donor=a, arterial_share_donor=a / (1 + a), arterial_load=phi, total_flow_load=total,
                 portal_share_graft=d['graft_portal'] / (d['graft_portal'] + d['graft_arterial']), zR_model=[float(min(model)), float(max(model))],
                 zR_model_all_resistance_after_the_junction=float(net.two_inflows(zF, 0.0, arterial_fraction=a / (1 + a), arterial_load=phi)['zR']))
        if d['gradient_mmHg'] is not None: r.update(pressure_load=d['gradient_mmHg'] / GRADIENT_MMHG, zR_observed=d['gradient_mmHg'] / GRADIENT_MMHG / zF)
        out[name] = r
    return out


def first_week():
    """The graft one week after transplantation: what the model gives for its portal pressure, with the measured growth
    in volume and the measured fall in portal inflow, against the measured course of the portal pressure."""
    out = dict(volume=VOLUME_DAY7, inflow=PORTAL_INFLOW_DAY10)
    for cal in ('unchanged', 'proportional'):
        r = regeneration(VOLUME_DAY7, calibre=cal, inflow=PORTAL_INFLOW_DAY10)
        out[cal] = dict(flow_load=float(r['zF']), pressure_load=float(r['zP']), resistance=float(r['resistance']))
    out['observed_portal_pressure'] = {k: v[1] / v[0] for k, v in PORTAL_PRESSURE_DAY7.items()}
    out['observed_gradient_if_caval_pressure_8'] = {k: (v[1] - 8) / (v[0] - 8) for k, v in PORTAL_PRESSURE_DAY7.items()}
    lo, hi = min(out['observed_gradient_if_caval_pressure_8'].values()), max(out['observed_portal_pressure'].values())
    out['resistance_implied'] = [lo / PORTAL_INFLOW_DAY10, hi / PORTAL_INFLOW_DAY10]
    return out


def summary():
    out = dict(calibration=dict(branch_length_factor=BRANCH_LENGTH_FACTOR, lobule_drop_mmHg=LOBULE_DROP_MMHG, gradient_mmHg=GRADIENT_MMHG,
                                portal_mL_min_100g=PORTAL_PER_100G, viscosity_Pa_s=VISCOSITY, b=B_HUMAN))
    out['cast_exponent'] = {t: {part: cast_exponent(t, part) for part in ('all', 'macro', 'meso')} for t in CAST}
    out['cast_poiseuille'] = {t: cast_poiseuille(t) for t in CAST}
    out['cast_tree_drop_mmHg'] = {t: {str(s): [cast_tree_drop(t, G, s) for G in CAST_GENERATIONS[t]] for s in (1.0, 1.1, 1.2)} for t in CAST_GENERATIONS}
    out['portal_exponent'] = portal_exponent(); out['portal_exponent_lobule_at_optimum'] = portal_exponent(lobule_drop=None)
    out['portal_exponent_space_filling_lengths'] = portal_exponent(kappa=1.0)
    out['portal_radii'] = [portal_radii(b) for b in (0.5, 2 / 3, 0.75, 0.9, 1.0)]
    out['cast_lengths_mm'] = dict(model=[float(x) for x in liver(B_HUMAN, liver_g=CAST_LIVER_G)['l_mm'][:11]], cast=[float(x) for x in cast_arrays('portal vein')[2]])
    keep = ('share_inlet', 'share_sinusoid', 'share_outlet', 'trunks_removed', 'lobule_drop_mmHg', 'portal_trunk_radius_mm', 'portal_trunk_length_mm',
            'portal_trunk_velocity_cm_s', 'portal_trunk_shear_Pa', 'portal_terminal_shear_Pa', 'hepatic_trunk_radius_mm', 'tree_volume_mL', 'tree_volume_fraction')
    out['adult'] = []
    for b in B_CALIBRATED:
        for label, drop in (('imposed', LOBULE_DROP_MMHG), ('optimum', None)):
            c = liver(b, lobule_drop=drop); row = dict(b=float(b), lobule=label, zR_lobar_graft=1 - c['trunks_removed']); row.update({k: c[k] for k in keep})
            row['zR_two_inflows'] = [net.two_inflows(zF, c['share_inlet'], 1 - c['trunks_removed'], arterial_load=phi)['zR'] for zF in (3.3, 4.1) for phi in (0.0, 1.0)]
            out['adult'].append(row)
    out['lobule_drop_sensitivity'] = []
    for drop in (0.35, 1.0, 2.2):
        c = liver(B_HUMAN, lobule_drop=drop); p = portal_radii(B_HUMAN, lobule_drop=drop)
        two = [net.two_inflows(zF, c['share_inlet'], 1 - c['trunks_removed'], arterial_load=phi)['zR'] for zF in (3.3, 4.1) for phi in (0.0, 1.0)]
        out['lobule_drop_sensitivity'].append(dict(lobule_drop_mmHg=drop, share_inlet=c['share_inlet'], share_sinusoid=c['share_sinusoid'], share_outlet=c['share_outlet'],
                                                   portal_model_over_cast=[min(p['model_over_cast']), max(p['model_over_cast'])], zR_two_inflows=[min(two), max(two)]))
    out['growth'] = [dict(b=float(b), s=float(s), exponents=growth(b, s)['exponents'], newborn=growth(b, s)['rows'][0], adult=growth(b, s)['rows'][-1])
                     for b in (0.6, 0.75, 0.93) for s in (0.0, 0.5, 1.0)]
    out['measured_growth'] = measured_growth()
    out['measured_growth']['pooled_diameter_ratio_adults_over_children'] = PORTAL_DIAMETER_MM['adults'][0] / PORTAL_DIAMETER_MM['children'][0]
    out['regeneration'] = []
    for cal in ('unchanged', 'proportional'):
        for months, rho in GRAFT_GROWTH.items():
            r = regeneration(rho, calibre=cal)
            out['regeneration'].append(dict(calibre=cal, months=months, rho=rho, flow_load=float(r['zF']), pressure_load=float(r['zP']), resistance_load=float(r['zR']),
                                            sinusoidal_velocity=float(r['sinusoidal_velocity'])))
    out['two_inflows_measured'] = two_inflows_measured(); out['first_week'] = first_week()
    d = np.array([x[2] for x in GRAFT_VOLUME]); v = np.array([x[3] for x in GRAFT_VOLUME])
    out['graft_volume'] = dict(series=[dict(series=a, n=n, days=t, ratio=float(x)) for a, n, t, x in GRAFT_VOLUME],
                               range_by_day={str(t): [float(v[d == t].min()), float(v[d == t].max())] for t in sorted(set(d))})
    out['donor_portal_per_100g'] = DONOR_PORTAL_PER_100G
    r = regeneration(3.0, zF0=3.0, calibre='proportional')
    out['regeneration_restored_mass'] = dict(rho=3.0, zF0=3.0, pressure_load=float(r['zP']), sinusoidal_velocity=float(r['sinusoidal_velocity']))
    return out


if __name__ == '__main__':
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'); os.makedirs(f'{root}/results', exist_ok=True)
    s = summary(); json.dump(s, open(f'{root}/results/model_human.json', 'w'), indent=1, default=float)
    print('exponent b read from the human cast (radius law of the optimum, measured radii, lengths and counts)')
    for t, parts in s['cast_exponent'].items():
        for part, v in parts.items():
            print(f"   {t:15s} {part:5s} generations {v['generations'][0]:>2}-{v['generations'][1]:<2}: b = {v['b']:.2f}" + (f" ({v['leave_one_out'][0]:.2f} to {v['leave_one_out'][1]:.2f} leaving one out)" if v['leave_one_out'] else '')
                  + f"; branch length / territory side {min(v['length_over_territory_side']):.2f} to {max(v['length_over_territory_side']):.2f}")
    print('Poiseuille on the measured generations at donor flow (mmHg per generation)')
    for t, v in s['cast_poiseuille'].items():
        print(f"   {t:15s} {' '.join(f'{x:.2f}' for x in v['drop_mmHg'])}  (x{v['drop_ratio_per_generation']:.2f} per generation; shear {v['shear_Pa'][0]:.2f} to {v['shear_Pa'][-1]:.2f} Pa)")
    for t, v in s['cast_tree_drop_mmHg'].items(): print(f"   whole {t}: {v['1.0'][0]:.1f} to {v['1.0'][1]:.1f} mmHg; {v['1.2'][0]:.1f} to {v['1.2'][1]:.1f} with radii 20% larger")
    p = s['portal_exponent']; print(f"absolute portal radii, nothing fitted: b = {p['b']:.2f} (rms {p['rms_log']:.3f}; below 0.10 for b in {p['interval'][0]:.2f} to {p['interval'][1]:.2f}; "
                                    f"{p['rms_at_half']:.2f} at b = 0.5, {p['rms_at_one']:.2f} at b = 1)")
    for r in s['portal_radii']: print(f"   b = {r['b']:.2f}: model / cast " + ' '.join(f'{x:.2f}' for x in r['model_over_cast']))
    print('adult liver of 1500 g, lobule drop imposed')
    for r in s['adult']:
        if r['lobule'] == 'imposed':
            print(f"   b = {r['b']:.2f}: drop {r['share_inlet']:.2f} / {r['share_sinusoid']:.2f} / {r['share_outlet']:.2f}; portal trunk radius {r['portal_trunk_radius_mm']:.1f} mm, length {r['portal_trunk_length_mm']:.0f} mm, "
                  f"velocity {r['portal_trunk_velocity_cm_s']:.0f} cm/s; shear {r['portal_trunk_shear_Pa']:.2f} to {r['portal_terminal_shear_Pa']:.2f} Pa; trees {r['tree_volume_fraction'] * 100:.0f}% of the volume; "
                  f"zR lobar graft {r['zR_lobar_graft']:.3f}, with two inflows {min(r['zR_two_inflows']):.2f} to {max(r['zR_two_inflows']):.2f}")
    print('growth from 150 g to 1500 g: exponents with liver volume')
    for r in s['growth']:
        e = r['exponents']; print(f"   b = {r['b']:.2f}, s = {r['s']:.1f}: portal diameter {e['portal_diameter_mm']:+.2f} ({r['newborn']['portal_diameter_mm']:.1f} to {r['adult']['portal_diameter_mm']:.1f} mm), "
                                  f"velocity {e['portal_velocity_cm_s']:+.2f} ({r['newborn']['portal_velocity_cm_s']:.0f} to {r['adult']['portal_velocity_cm_s']:.0f} cm/s), lobule drop {e['lobule_drop_mmHg']:+.2f}, "
                                  f"sinusoidal velocity {e['sinusoidal_velocity']:+.2f}")
    m = s['measured_growth']; print(f"   measured (weight groups, ranges): portal diameter {m['portal_diameter']['exponent']:+.2f} ({m['portal_diameter']['interval'][0]:+.2f} to {m['portal_diameter']['interval'][1]:+.2f}), "
                                    f"velocity {m['portal_velocity']['exponent']:+.2f} ({m['portal_velocity']['interval'][0]:+.2f} to {m['portal_velocity']['interval'][1]:+.2f})")
    print('both inflows measured in the donor and in the graft')
    for k, v in s['two_inflows_measured'].items():
        print(f"   {k}: flow load {v['flow_load']:.2f}, arterial share of the donor {v['arterial_share_donor']:.2f}, arterial flow per gram x{v['arterial_load']:.2f}, total flow load {v['total_flow_load']:.2f}; "
              f"model zR {v['zR_model'][0]:.2f} to {v['zR_model'][1]:.2f}" + (f"; observed {v['zR_observed']:.2f}" if 'zR_observed' in v else ''))
    w = s['first_week']; print(f"first week: volume x{w['volume']:.2f}, portal inflow x{w['inflow']:.2f}; portal pressure x{w['unchanged']['pressure_load']:.2f} with calibres unchanged, x{w['proportional']['pressure_load']:.2f} in proportion; "
                               f"measured x{min(w['observed_gradient_if_caval_pressure_8'].values()):.2f} to x{max(w['observed_portal_pressure'].values()):.2f}")
    print('graft regrowing by lobule enlargement at constant inflow (loads relative to reperfusion)')
    for r in s['regeneration']: print(f"   calibres {r['calibre']:12s} {r['months']:>4} months, volume x{r['rho']:.2f}: flow load x{r['flow_load']:.2f}, pressure load x{r['pressure_load']:.2f}, resistance load x{r['resistance_load']:.2f}")
