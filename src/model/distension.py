"""The resistance the rigid network does not explain: what a passively distensible bed can account for.

    python src/model/distension.py      -> results/model_distension.json

The network of src/model/ is rigid: its calibres are those of the optimum and do not change with pressure. It gives a
partial graft a resistance load of 0.8 to 1 (liver.py, human.py), and the six published groups that report both loads
lie lower (0.31 to 0.74). This module asks how much of that distance each mechanism the network leaves out can carry.
Nothing here is fitted to the six groups: every number that enters is the model's own or a published one.

Three mechanisms, each as a bound.

sinusoids     Dilation or recruitment of sinusoids acts on the share of the pressure drop the sinusoids carry. In the
              optimum that share is 3 to 15% (b from 1 to 0.5, 14 generations), so with unlimited dilation the
              resistance load cannot fall below 0.85 to 0.97 of the rigid value. Not enough on its own.
outflow       Widening the outlet anastomosis acts on the 5% of the resistance that anastomosis holds
              (whole_graft.py): 0.95 at most. Not enough on its own.
distension    The venous resistance of the normal liver is passively distensible. Three published relations for the
              whole liver, all measured in cats and none on a graft:
                cube law     resistance proportional to 1 / Pd^3 (index of contractility, IC = R Pd^3), Pd the mean of
                             the pressures on either side of the resistance (Lautt, Greenway and Legare, Am J Physiol
                             1991;260:G325, validated over 10 to 50 mL/min/kg, a third to 1.7 times the control flow, and
                             outflow pressures up to 15 mmHg; Lautt and Legare, Am J Physiol 1992;263:G702);
                compliance   hepatic blood volume linear in intrahepatic pressure. It more than doubled for 9.4 mmHg in
                             Lautt and Greenway (Am J Physiol 1976;231:292: +32 on 27 mL/100 g), "8 mmHg will double"
                             in Lautt and Legare 1992, and a compliance of 2.6 to 3.1 mL/mmHg/100 g on 34 mL/100 g,
                             a doubling at 11 to 13 mmHg, in Greenway, Seaman and Innes (Am J Physiol 1985;248:H468).
                             The doubling pressure is therefore taken over 8 to 13 mmHg. With the length of the vessels
                             fixed and Poiseuille flow, resistance goes as volume^-2: this takes the resistance vessels
                             to distend as the capacitance vessels do, an assumption and an upper bound;
                pressure-flow intrahepatic pressure linear in total hepatic flow up to 240 mL/min/100 g of portal flow,
                             about 2.2 times the control total flow, with an intercept of 2.2 mmHg at zero flow and a
                             slope of 0.030 mmHg per mL/min/100 g relative to the caval pressure at a control flow of
                             139 (Greenway, Seaman and Innes 1985; nonlinear fits were worse). Expressed on the
                             portocaval gradient (see GREENWAY), the share of the control gradient that does not scale
                             with flow is f = 0.30 to 0.35, and G / g0 = f + (1 - f) q at a flow load q. It is the
                             relation measured at the highest flows, and it stops distending: resistance tends to the
                             slope as flow rises.
              The first two keep distending at the flows of a graft, the third does not. Their authors note that in the
              normal liver the changes of resistance are notable below 10 mmHg and minor above (Ingles, Legare and
              Lautt, Can J Physiol Pharmacol 1994;72:687), and a graft is at 14 to 22 mmHg. The three together are a
              range, not a prediction.
check         `published_check` sets the relations against pig livers they were not derived from (Ricciardi et al.,
              J Gastrointest Surg 2001;5:490): in situ and isolated, and the same livers after 2 h of cold ischaemia,
              when resistance was twice as high at every flow. For each state it gives the error of the predicted
              resistance ratio and of the predicted portal pressure at the higher flow, and checks that the caval
              pressure implied by the printed pressures and resistances is the same at both flows.

For a group with flow load zF and central venous pressure CVP each relation is solved self-consistently: the gradient G
satisfies G = 5 zF zR_rigid d(G), with d the relation and the donor at a gradient of 5 mmHg; the pressure-flow relation
is applied to the effective load q = zF zR_rigid. The inputs are the flow load, the rigid network and the central venous
pressure; the portal pressure of the recipient is the output and is not used. The rigid value runs from the network with
the two inflows (0.81) to the network alone (1.00), because the relations were measured with the hepatic artery in place.

What comes out is the range of zR the passive relations give, set beside the rigid network. Each group is placed with its
95% sampling interval: below the relations, within their range, or above it. That placement is descriptive. Three further
analyses ask how much of it depends on choices:
robustness    `robustness` (with src/model/compartments.py): the relations applied vessel by vessel at the pressure the
              network gives each, with only the pre- and postsinusoidal sites, only one of them, or every vessel
              distending, the share that distends fixed by the network; and the global relations with the distending
              pressure CVP - Pext + lambda G for lambda from 0.25 to 1.
joint         `joint`: measurement, references, network and parameters drawn together; for each group the probability
              that the observation lies below the rigid network, below every relation, within them or above them.
The relations describe normal cat livers, mostly below about 10 mmHg. A reperfused human graft differs by species,
ischaemia-reperfusion, endothelial injury and vascular tone, and the pig livers after cold ischaemia distended less than
the relations give for their pressure. The relations are therefore read as the largest plausible passive contribution,
subject to those differences, not as a prediction of the graft.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import indices
from model import liver, human

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
DOUBLING_MMHG = (8.0, 13.0)              # intrahepatic pressure that doubles hepatic blood volume (Lautt and Legare 1992; Greenway et al. 1985)
CUBE = 3.0                               # exponent of the index of contractility, IC = R Pd^3 (Lautt, Greenway and Legare 1991)
# Greenway et al. 1985 regressed the intrahepatic pressure, the mean of the portal and the lobar venous pressures, on
# total hepatic flow: intercept 2.2 mmHg at zero flow, slope 0.046 mmHg per mL/min/100 g, of which 0.015 is the outflow
# cannula raising the caval pressure with flow, so 0.030 relative to the caval pressure; control flow 139. The loads of
# this repository use the portocaval gradient, portal minus caval pressure, which exceeds the intrahepatic-to-caval
# difference by half the portal-to-lobar difference, 1 to 2 mmHg in that paper. With that half taken proportional to
# flow, as in a rigid bed, the share of the control portocaval gradient that does not scale with flow is
#     f = 2.2 / (2.2 + 0.030 x 139 + d / 2),   d from 0 to 2 mmHg          -> 0.30 to 0.35
# The slope of 0.046 is not used: it includes the rise of the caval pressure, which is not part of the gradient.
GREENWAY = dict(intercept=2.2, slope=0.030, slope_with_the_cannula=0.046, control_flow=139.0, portal_to_lobar=(0.0, 2.0))
FLOW_FIXED_SHARE = tuple(sorted(round(GREENWAY['intercept'] / (GREENWAY['intercept'] + GREENWAY['slope'] * GREENWAY['control_flow'] + dd / 2), 3) for dd in GREENWAY['portal_to_lobar']))
CVPS = (3.0, 5.0, 7.0, 9.0)              # central venous pressure where a series reports none, as in predictions.REF_CVPS
OUTLET_SHARE = 0.05                      # share of the resistance held by the outlet anastomosis (whole_graft.py, r_out)

# The three published relations. The first two give resistance relative to the reference state as a function of the
# gradient G, the central venous pressure and the reference gradient g0; the third is linear in flow and is solved directly.
#   cube law      distending pressure the mean of the pressures on either side of the resistance, CVP + G / 2
#   compliance    distending pressure the intrahepatic pressure, which the portal venous pressure estimates, CVP + G
#   pressure-flow G / g0 = f + (1 - f) q
LAWS = {
    'cube law': lambda G, cvp, g0, p=None: ((cvp + g0 / 2.0) / (cvp + G / 2.0)) ** CUBE,
    'compliance': lambda G, cvp, g0, p: (1.0 + (G - g0) / p) ** -2.0,
    'pressure-flow': None,
}
PARAMS = {'cube law': (None,), 'compliance': DOUBLING_MMHG, 'pressure-flow': FLOW_FIXED_SHARE}
FAMILY = {'cube law': 'keeps distending', 'compliance': 'keeps distending', 'pressure-flow': 'stops distending'}


def distensible(zF, zR_rigid, cvp, law, p=None, normal_gradient=indices.NORMAL_GRADIENT):
    """Gradient and loads of a bed that follows `law` (parameter p), at flow load zF."""
    if law == 'pressure-flow':
        f = FLOW_FIXED_SHARE[0] if p is None else p; G = normal_gradient * (f + (1 - f) * zF * zR_rigid)
        return dict(gradient_mmHg=float(G), zP=float(G / normal_gradient), zR=float(G / normal_gradient / zF))
    p = PARAMS[law][0] if p is None else p; d = lambda G: LAWS[law](G, cvp, normal_gradient, p); lo, hi = 0.0, normal_gradient * zF * zR_rigid * 2 + 1
    f = lambda G: normal_gradient * zF * zR_rigid * d(G) - G        # decreasing in G above the reference gradient
    for _ in range(80):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    G = 0.5 * (lo + hi); return dict(gradient_mmHg=float(G), zP=float(G / normal_gradient), zR=float(G / normal_gradient / zF))


def over(law, fn):
    """fn applied over the parameter range of a relation: (lowest, highest)."""
    v = [fn(p) for p in PARAMS[law]]; return min(v), max(v)


# Pressure and flow of normal livers of a third species, not used to derive any relation: the pig livers of Ricciardi
# et al. (J Gastrointest Surg 2001;5:490), in situ and isolated on a bypass circuit, and the same livers after 2 h of
# cold ischaemia. Portal flow in mL/min, portal venous pressure in mmHg, resistance (portal minus caval pressure over
# portal flow) in mmHg min/mL, at the lower and the higher flow, each with its standard error, as printed in the
# Results of that paper.
RICCIARDI = {
    'in situ': dict(n=9, flow=(400.0, 1100.0), pvp=(6.8, 9.6), pvp_se=(0.7, 0.24), hvr=(0.016, 0.0086), hvr_se=(0.0018, 0.0004)),
    'isolated, perfused': dict(n=7, flow=(300.0, 800.0), pvp=(5.14, 8.14), pvp_se=(0.3, 0.46), hvr=(0.017, 0.010), hvr_se=(0.001, 0.0006)),
    'after 2 h of cold ischaemia': dict(n=6, flow=(300.0, 800.0), pvp=(10.2, 15.8), pvp_se=(1.2, 1.3), hvr=(0.034, 0.020), hvr_se=(0.0039, 0.0016)),
}


def published_check():
    """The relations against measurements they were not derived from. For each state of the pig liver the flow is
    raised by the published factor from the published starting gradient and caval pressure, and each relation gives
    the resistance and the portal pressure at the higher flow. Reported for each state and relation: the predicted
    resistance ratio against the measured one (difference, and the measured ratio's standard error, independent
    errors assumed), the predicted portal pressure against the measured one (error in mmHg, and in standard errors of
    the measurement), and whether the caval pressure implied by pressure and resistance is the same at both flows,
    which it must be if the printed values are consistent."""
    out = {}
    for name, r in RICCIARDI.items():
        g0 = r['hvr'][0] * r['flow'][0]; cvp = r['pvp'][0] - g0; zF = r['flow'][1] / r['flow'][0]
        cvp_high = r['pvp'][1] - r['hvr'][1] * r['flow'][1]; meas = r['hvr'][1] / r['hvr'][0]
        meas_se = meas * float(np.hypot(r['hvr_se'][0] / r['hvr'][0], r['hvr_se'][1] / r['hvr'][1])); pred = {}
        for law in LAWS:
            sol = [distensible(zF, 1.0, cvp, law, p, normal_gradient=g0) for p in PARAMS[law]]
            z = [x['zR'] for x in sol]; pv = [cvp + x['gradient_mmHg'] for x in sol]
            err_r = [zz - meas for zz in z]; err_p = [v - r['pvp'][1] for v in pv]
            pred[law] = dict(resistance_ratio=[round(min(z), 3), round(max(z), 3)], resistance_ratio_error=[round(min(err_r), 3), round(max(err_r), 3)],
                             portal_pressure_mmHg=[round(min(pv), 2), round(max(pv), 2)], portal_pressure_error_mmHg=[round(min(err_p), 2), round(max(err_p), 2)],
                             portal_pressure_error_in_se=[round(min(err_p) / r['pvp_se'][1], 1), round(max(err_p) / r['pvp_se'][1], 1)],
                             within_two_se_of_the_measured_ratio=bool(min(abs(e) for e in err_r) <= 2 * meas_se or (min(err_r) <= 0 <= max(err_r))))
        out[name] = dict(n=r['n'], flow_factor=round(zF, 3), starting_gradient_mmHg=round(g0, 2), caval_pressure_mmHg=round(cvp, 2),
                         caval_pressure_at_the_higher_flow_mmHg=round(cvp_high, 2), caval_pressure_consistent=bool(abs(cvp_high - cvp) < 0.5),
                         measured_resistance_ratio=round(meas, 3), measured_resistance_ratio_se=round(meas_se, 3),
                         measured_portal_pressure_mmHg=r['pvp'][1], measured_portal_pressure_se_mmHg=r['pvp_se'][1], predicted=pred,
                         resistance_relative_to_the_isolated_liver_at_low_flow=round(r['hvr'][1] / RICCIARDI['isolated, perfused']['hvr'][0], 3))
    normal = [k for k in RICCIARDI if 'ischaemia' not in k]
    mae = {law: dict(normal_livers=round(float(np.mean([max(abs(e) for e in out[k]['predicted'][law]['portal_pressure_error_mmHg']) for k in normal])), 2),
                     after_cold_ischaemia=round(float(max(abs(e) for e in out['after 2 h of cold ischaemia']['predicted'][law]['portal_pressure_error_mmHg'])), 2)) for law in LAWS}
    return dict(states=out, largest_portal_pressure_error_mmHg=mae,
                note='normal livers: in situ and isolated; the ischaemic state is the same livers after 2 h at 4 C. Each relation is anchored at the '
                     'lower flow of each state, so a baseline doubled by ischaemia does not enter; what is tested is how resistance falls with flow')


def sinusoidal_share(bs=(1.0, 0.75, 0.5), g=liver.HUMAN_GENERATIONS):
    """Share of the pressure drop the sinusoids carry in the optimum, over the range of b."""
    return {b: float(liver.canopy(g, b)['share_sinusoid']) for b in bs}


def rigid_range(zF, alone=True):
    """The rigid network at flow load zF, over b in [0.5, 1] and the calibrated human liver: with the two inflows
    (arterial flow per lobule between none and the donor's) and, if `alone`, also the network by itself (arterial
    flow scaling with the portal one), which is its upper end."""
    lob = [liver.lobar_graft(liver.HUMAN_GENERATIONS, b, 1) for b in (1.0, 0.75, 0.5)]; cal = [human.liver(b) for b in (0.6, 2 / 3, 0.75, 0.85, 0.93)]
    phis = (0.0, 1.0, zF) if alone else (0.0, 1.0)
    z = [liver.two_inflows(zF, l['share_inlet'], l['zR'], arterial_load=phi)['zR'] for l in lob for phi in phis]
    z += [liver.two_inflows(zF, c['share_inlet'], 1 - c['trunks_removed'], arterial_load=phi)['zR'] for c in cal for phi in phis]
    return float(min(z)), float(max(z))


def groups(d):
    """The six groups that report both loads, with the central venous pressure each reports or implies."""
    rows = []
    for _, r in d.dropna(subset=['zF', 'zP']).iterrows():
        if pd.notna(r.gradient_mmHg) and pd.notna(r.PVP): cvp, basis = float(r.PVP - r.gradient_mmHg), 'portal pressure minus measured gradient'
        elif r.CVP_basis == 'imputed': cvp, basis = None, 'not reported'
        else: cvp, basis = float(r.CVP), str(r.CVP_basis)
        rows.append(dict(study=r.study, group=r.group, zF=float(r.zF), zP=float(r.zP), zR=float(r.zR), PVP=float(r.PVP) if pd.notna(r.PVP) else None, cvp=cvp, cvp_basis=basis))
    return rows


def _place(value, ci, region):
    """Where an observation stands against a region [lo, hi] of the model: the point estimate and the 95% sampling
    interval are reported apart, with one rule for every region (the interval reaches the region if they intersect)."""
    lo, hi = region
    return dict(point='below' if value < lo else ('above' if value > hi else 'inside'),
                interval='below' if ci[1] < lo else ('above' if ci[0] > hi else 'reaches'))


def bracket(d):
    """Each group against four regions of the model: the rigid network, the relations that keep distending, the one
    that stops, and the span of all three. For each region the point estimate and the sampling interval are placed
    with the same rule (`_place`). This is a description: a sampling interval reaching a region is not a probability
    that the mechanism of that region acts. `joint` gives the draws of every uncertainty together."""
    import audit
    CI = {(x['study'], x['group']): x['zR_ci_sampling'] for x in audit.paired(d, correlations=())['groups']}; rows = []
    for g in groups(d):
        lo, hi = rigid_range(g['zF']); cvps = (g['cvp'],) if g['cvp'] is not None else CVPS; sol = {}
        for law in LAWS:
            z = [distensible(g['zF'], zr, c, law, p)['zR'] for zr in (lo, hi) for c in cvps for p in PARAMS[law]]; sol[law] = [round(min(z), 3), round(max(z), 3)]
        fam = lambda which: [min(sol[l][0] for l in LAWS if FAMILY[l] == which), max(sol[l][1] for l in LAWS if FAMILY[l] == which)]
        regions = {'rigid network': [round(lo, 3), round(hi, 3)], 'keeps distending': fam('keeps distending'), 'stops distending': fam('stops distending'),
                   'all relations': [min(v[0] for v in sol.values()), max(v[1] for v in sol.values())]}
        ci = [round(float(CI[(g['study'], g['group'])][0]), 3), round(float(CI[(g['study'], g['group'])][1]), 3)]
        place = {k: _place(g['zR'], ci, v) for k, v in regions.items()}
        rows.append(dict(study=g['study'], group=g['group'], zF=round(g['zF'], 3), zR_observed=round(g['zR'], 3), cvp_mmHg=g['cvp'], cvp_basis=g['cvp_basis'],
                         zR_rigid=regions['rigid network'], zR_distensible=sol, zR_distensible_all=regions['all relations'],
                         zR_keeps_distending=regions['keeps distending'], zR_stops_distending=regions['stops distending'], zR_observed_ci=ci, placement=place,
                         position={'below': 'below the relations', 'above': 'above the relations', 'reaches': 'within the relations'}[place['all relations']['interval']],
                         relations_containing_the_value=[law for law, v in sol.items() if v[0] <= g['zR'] <= v[1]],
                         below_the_rigid_network=bool(place['rigid network']['interval'] == 'below')))
    obs = np.array([r['zR_observed'] for r in rows]); cnt = lambda reg, kind, val: int(sum(r['placement'][reg][kind] == val for r in rows))
    counts = {reg: dict(point_inside=cnt(reg, 'point', 'inside'), point_below=cnt(reg, 'point', 'below'), point_above=cnt(reg, 'point', 'above'),
                        interval_reaches=cnt(reg, 'interval', 'reaches'), interval_below=cnt(reg, 'interval', 'below'), interval_above=cnt(reg, 'interval', 'above'))
              for reg in ('rigid network', 'keeps distending', 'stops distending', 'all relations')}
    return dict(groups=rows, counts=counts, mean_observed=round(float(obs.mean()), 3),
                mean_rigid=[round(float(np.mean([r['zR_rigid'][k] for r in rows])), 3) for k in (0, 1)],
                mean_distensible={law: [round(float(np.mean([r['zR_distensible'][law][k] for r in rows])), 3) for k in (0, 1)] for law in LAWS},
                mean_distensible_all=[round(float(np.mean([r['zR_distensible_all'][k] for r in rows])), 3) for k in (0, 1)],
                within_the_relations=counts['all relations']['interval_reaches'],
                below_the_relations=[f"{r['study']} / {r['group']}" for r in rows if r['position'] == 'below the relations'],
                above_the_relations=[f"{r['study']} / {r['group']}" for r in rows if r['position'] == 'above the relations'],
                below_the_rigid_network=counts['rigid network']['interval_below'],
                value_in_a_relation_that_keeps_distending=counts['keeps distending']['point_inside'],
                value_in_the_relation_that_stops_distending=counts['stops distending']['point_inside'],
                note='placement is descriptive: point estimate and sampling interval against ranges of the model, one rule for every region')


def consistency(zFs=np.round(np.geomspace(0.25, 20, 121), 4), cvps=CVPS):
    """The three relations share one definition of pressure (the portocaval gradient, in units of the donor's) and of
    flow (the portal flow load, through the effective load q = zF zR_rigid of the rigid network), and are composed with
    the network in the same way: the relation acts on the load the rigid network would carry. Checked here:
    reference     with the whole donor network (zR_rigid = 1) at zF = 1 every relation gives zP = zF = zR = 1 at every CVP
                  and parameter, so the network and the relation do not both count the donor state
    no double count  the relation sees the rigid network only through q: two networks with the same q give the same
                  gradient, so architecture and arterial share enter once
    monotone      zP rises and zR falls with flow, with no jump on a fine grid of zF from 0.25 to 20
    asymptote     the local exponent d log zP / d log zF at zF = 20: 1/4 for the cube law and 1/3 for compliance
                  (resistance keeps falling), 1 for the pressure-flow relation (resistance tends to (1 - f) zR_rigid)"""
    out = {}
    for law in LAWS:
        ref = [distensible(1.0, 1.0, c, law, p) for c in cvps for p in PARAMS[law]]
        same_q = max(abs(distensible(2.0, 0.9, c, law, p)['gradient_mmHg'] - distensible(1.8, 1.0, c, law, p)['gradient_mmHg']) for c in cvps for p in PARAMS[law])
        mono, jump, expo = True, 0.0, []
        for c in cvps:
            for p in PARAMS[law]:
                zp = np.array([distensible(float(z), 1.0, c, law, p)['zP'] for z in zFs]); zr = zp / zFs
                mono = mono and bool(np.all(np.diff(zp) > 0) and np.all(np.diff(zr) < 0))
                jump = max(jump, float(np.max(np.abs(np.diff(np.log(zp))) / np.diff(np.log(zFs)))))
                expo.append(float(np.log(zp[-1] / zp[-2]) / np.log(zFs[-1] / zFs[-2])))
        out[law] = dict(reference_error=float(max(max(abs(r['zP'] - 1), abs(r['zR'] - 1)) for r in ref)), same_load_same_gradient_mmHg=round(same_q, 12),
                        monotone=mono, largest_local_exponent=round(jump, 3), exponent_at_the_highest_flow=[round(min(expo), 3), round(max(expo), 3)])
    return dict(relations=out, pressure='portocaval gradient, PVP - CVP, divided by the donor gradient', flow='portal flow per 100 g divided by the donor value',
                composition='each relation acts on q = zF zR_rigid, the gradient ratio of the rigid network; at q = 1 every relation returns the donor')


def heterogeneity(d, B=None):
    """Do the relations order the groups? They depend on the flow load and the CVP of each group and on parameters common
    to all, so they can set the level of zR but need not reproduce its spread or its order. Reported: the spread of the
    observed zR and of each relation at its central parameters (log scale), and the rank correlation between them."""
    from scipy.stats import spearmanr
    B = B or bracket(d); gs = {(g['study'], g['group']): g for g in groups(d)}; obs = np.array([r['zR_observed'] for r in B['groups']]); out = {}
    mid = {'cube law': None, 'compliance': float(np.mean(DOUBLING_MMHG)), 'pressure-flow': float(np.mean(FLOW_FIXED_SHARE))}
    for law in LAWS:
        z = []
        for r in B['groups']:
            g = gs[(r['study'], r['group'])]; lo, hi = rigid_range(g['zF'], alone=False); c = g['cvp'] if g['cvp'] is not None else 5.0
            z.append(distensible(g['zF'], 0.5 * (lo + hi), c, law, mid[law])['zR'])
        z = np.array(z); rho = spearmanr(obs, z)
        out[law] = dict(zR=[round(float(x), 3) for x in z], range=[round(float(z.min()), 3), round(float(z.max()), 3)], sd_log=round(float(np.std(np.log(z), ddof=1)), 3),
                        spearman_rho=round(float(rho.statistic), 3), spearman_p=round(float(rho.pvalue), 3))
    return dict(groups=[f"{r['study']} / {r['group']}" for r in B['groups']], observed=dict(zR=[round(float(x), 3) for x in obs], range=[round(float(obs.min()), 3), round(float(obs.max()), 3)], sd_log=round(float(np.std(np.log(obs), ddof=1)), 3)),
                relations=out, note='central parameters, the rigid network at the middle of its two-inflow range, CVP 5 mmHg where none is reported')


def wang_contrast(d, draws=20000, seed=11):
    """The one contrast inside a cohort: the two arms of Wang 2014, splenectomy against none. Both arms are given the same
    CVP, the same network and the same parameters in each draw (one cohort, one centre), and each relation is solved from
    the flow load of each arm only; their pressures are not used. Compared: the ratio of zR between the arms, observed
    (sampling error of the four means and the CVP common to both) and from each relation."""
    import audit
    rng = np.random.default_rng(seed); arms = {r.group: r for _, r in d.dropna(subset=['zF', 'zP']).iterrows() if r.study == 'Wang 2014'}
    A = {(g['study'], g['group']): g for g in audit.paired(d, correlations=())['groups']}; rep = {k: audit.REPORTED[('Wang 2014', k)] for k in arms}
    from model import compartments as cm
    G8 = cm.geometries(); s_in = np.array([x[1][:x[2]].sum() for x in G8]); rem = np.array([x[3] for x in G8]); a = cm.A_RATIO
    cvp = rng.uniform(CVPS[0], CVPS[-1], draws); g0 = rng.uniform(*NORMAL_GRADIENTS, draws); ref = rng.uniform(*DONOR_FLOW, draws)
    i = rng.integers(0, len(G8), draws); phi = rng.uniform(0, 1, draws); p = rng.uniform(*DOUBLING_MMHG, draws); f = rng.uniform(*FLOW_FIXED_SHARE, draws)
    zr_obs, zr_mod = {}, {}
    for k, r in arms.items():
        g = A[('Wang 2014', k)]
        zF = r.PVF_per_100g * np.exp(g['flow_cv_used'] / np.sqrt(rep[k]['n_flow']) * rng.standard_normal(draws)) / ref
        grad = r.PVP + g['gradient_sd_used'] / np.sqrt(rep[k]['n_grad']) * rng.standard_normal(draws) - cvp
        zr_obs[k] = grad / g0 / zF
        rig = (1 - rem[i]) * (s_in[i] * zF + (1 - s_in[i]) * (zF + phi * a) / (1 + a)) / zF; Q = g0 * zF * rig
        zr_mod[k] = {'rigid network': rig,
                     'cube law': _bisect(lambda G: Q * ((cvp + g0 / 2) / (cvp + G / 2)) ** CUBE - G, np.zeros(draws), 2 * Q + 1) / g0 / zF,
                     'compliance': _bisect(lambda G: Q * np.maximum(1 + (G - g0) / p, 1e-9) ** -2.0 - G, np.zeros(draws), 2 * Q + 1) / g0 / zF,
                     'pressure-flow': (f + (1 - f) * zF * rig) / zF}
    s_, n_ = 'splenectomy', 'no splenectomy'; q = lambda x: [round(float(v), 3) for v in np.percentile(x, [2.5, 50, 97.5])]
    obs = zr_obs[s_] / zr_obs[n_]; out = dict(flow_load=[round(float(arms[s_].zF), 3), round(float(arms[n_].zF), 3)], zR_observed=[round(float(arms[s_].zR), 3), round(float(arms[n_].zR), 3)],
                                               observed_ratio=q(obs), models={})
    for m in zr_mod[s_]:
        pr = zr_mod[s_][m] / zr_mod[n_][m]
        out['models'][m] = dict(ratio=q(pr), share_of_the_log_contrast=round(float(np.median(np.log(pr)) / np.median(np.log(obs))), 2),
                                prob_observed_ratio_below_the_model=round(float((obs < pr).mean()), 3))
    out['note'] = ('the arms were not randomised and splenectomy was indicated by the haemodynamics it changes; the relations see only the flow loads, '
                   'so what they do not reproduce is what the flow load of the arms does not carry')
    return out


def band(zFs=np.round(np.arange(1.0, 5.01, 0.1), 2)):
    """The passive relations in the flow-pressure plane: lowest and highest pressure load at each flow load, over the
    rigid range, a central venous pressure of 3 to 9 mmHg and the parameters, for the relations that keep distending
    and for the one that stops."""
    out = {'zF': [float(x) for x in zFs]}
    for fam in ('keeps distending', 'stops distending'):
        lo, hi = [], []
        for zF in zFs:
            r = rigid_range(float(zF)); z = [distensible(float(zF), zr, c, law, p)['zP'] for law in LAWS if FAMILY[law] == fam for zr in r for c in CVPS for p in PARAMS[law]]; lo.append(min(z)); hi.append(max(z))
        out[fam] = dict(zP_low=[round(float(x), 4) for x in lo], zP_high=[round(float(x), 4) for x in hi])
    out['zP_low'] = out['keeps distending']['zP_low']; out['zP_high'] = out['stops distending']['zP_high']
    return out


def required(zR_obs=0.584, zR_rigid=0.853):
    """With no relation assumed: the widening that closes the distance if a share s of the resistance widens uniformly
    (Poiseuille, R as radius^-4). None where even an unlimited widening of that share is not enough."""
    target = zR_obs / zR_rigid; out = []
    for s in (0.15, 0.32, 0.46, 0.5, 0.75, 1.0):
        x = (target - (1 - s)) / s
        out.append(dict(share_that_widens=s, resistance_of_that_share=None if x <= 0 else round(float(x), 3), radius_factor=None if x <= 0 else round(float(x ** -0.25), 3)))
    return dict(ratio_observed_to_rigid=round(float(target), 3), reduction_missing=round(float(1 - target), 3), smallest_share_that_can_close_it=round(float(1 - target), 3), by_share=out)


def robustness(d, B=None):
    """Do the conclusions survive other admissible ways of distending? Each scenario replaces the two relations that
    keep distending (the pressure-flow relation is a whole-organ measurement and is kept); each group is placed with
    its sampling interval against the span of the relations in that scenario, as in `bracket`.
      published          the whole-liver relations of `bracket`
      site: ...          compartments.site_ranges: the vessels of the network distend at their own mean pressure,
                         only those of the named sites, with the share s that distends fixed by the network
      lambda = x         compartments.pressure_variants: Pd = CVP - Pext + lambda G, Pext 0 or 2 mmHg"""
    from model import compartments as cm
    B = B or bracket(d); gs = groups(d); key = lambda r: (r['study'], r['group'])
    pf = {key(r): r['zR_distensible']['pressure-flow'] for r in B['groups']}; ci = {key(r): r['zR_observed_ci'] for r in B['groups']}
    obs = {key(r): r['zR_observed'] for r in B['groups']}; rig = {key(r): r['zR_rigid'] for r in B['groups']}
    scen = {'published': {key(r): [r['zR_keeps_distending'][0], r['zR_keeps_distending'][1]] for r in B['groups']}}
    S = cm.site_ranges(gs, CVPS, DOUBLING_MMHG); V = cm.pressure_variants(gs, CVPS, DOUBLING_MMHG, rigid_range); shares = {}
    for site in cm.SITES:
        scen[f'site: {site}'] = {key(r): [min(r['by_site'][site][l][0] for l in ('cube law', 'compliance')), max(r['by_site'][site][l][1] for l in ('cube law', 'compliance'))] for r in S}
        shares[site] = S[0]['by_site'][site]['share_that_distends']
    for lam in cm.LAMBDAS:
        scen[f'lambda = {lam:g}'] = {key(r): [min(r['by_lambda'][f'{lam:g}'][l][0] for l in ('cube law', 'compliance')), max(r['by_lambda'][f'{lam:g}'][l][1] for l in ('cube law', 'compliance'))] for r in V}
    out = {}
    for name, kd in scen.items():
        rows = []
        for k, v in kd.items():
            lo, hi = min(v[0], pf[k][0]), max(v[1], pf[k][1]); c = ci[k]
            pos = 'below the relations' if c[1] < lo else ('above the relations' if c[0] > hi else 'within the relations')
            rows.append(dict(study=k[0], group=k[1], zR_observed=obs[k], zR_observed_ci=c, zR_rigid=rig[k], zR_keeps_distending=[round(v[0], 3), round(v[1], 3)],
                             zR_all=[round(lo, 3), round(hi, 3)], position=pos, interval_reaches_a_relation_that_keeps_distending=bool(c[0] <= v[1] and c[1] >= v[0])))
        out[name] = dict(groups=rows, within=int(sum(r['position'] == 'within the relations' for r in rows)),
                         below=[f"{r['study']} / {r['group']}" for r in rows if r['position'] == 'below the relations'],
                         above=[f"{r['study']} / {r['group']}" for r in rows if r['position'] == 'above the relations'],
                         reach_a_relation_that_keeps_distending=int(sum(r['interval_reaches_a_relation_that_keeps_distending'] for r in rows)),
                         keeps_distending_span=[round(min(r['zR_keeps_distending'][0] for r in rows), 3), round(max(r['zR_keeps_distending'][1] for r in rows), 3)],
                         rigid_above_every_group=bool(all(r['zR_observed_ci'][1] < r['zR_rigid'][0] for r in rows)))
    return dict(scenarios=out, share_that_distends=shares, site_ranges=S, pressure_variants=V,
                note='the share that distends is fixed by the network (portal tree, sinusoids, hepatic venous tree), not by the groups')


# Ranges drawn from in `joint`; each is the range used elsewhere in the repository.
DONOR_FLOW = (80.0, 100.0)            # predictions.REF_FLOWS
NORMAL_GRADIENTS = (4.0, 6.0)         # predictions.REF_GRADIENTS


def _bisect(fun, lo, hi, n=70):
    for _ in range(n):
        mid = 0.5 * (lo + hi); pos = fun(mid) > 0; lo = np.where(pos, mid, lo); hi = np.where(pos, hi, mid)
    return 0.5 * (lo + hi)


def joint(d, draws=20000, seed=7):
    """Every source of uncertainty drawn together, group by group.

    measurement  the two group means, each with the dispersion the paper reports (the largest reported where none
                 is), independent; for a gradient derived from portal pressure, the CVP the paper gives or, where it
                 gives none, uniform over 3 to 9 mmHg
    references   the donor portal flow per 100 g (80 to 100, or the same relative range around a value measured in
                 the series) and the normal gradient (4 to 6 mmHg)
    network      one of the eight networks of the rigid range, and the arterial flow per lobule uniform between none
                 and the donor's (the two inflows of P2, whose upper end is below that of the network alone, so the
                 probability of lying below the network is not inflated)
    relations    the doubling pressure uniform over 8 to 13 mmHg, the flow-independent share uniform over 0.26 to 0.35
    The observed zR and the zR of each model are computed from the same draw, so the references they share cancel where
    they should. The references are common to the groups, so these are statements about one group at a time."""
    import audit
    rng = np.random.default_rng(seed); A = {(g['study'], g['group']): g for g in audit.paired(d, correlations=())['groups']}
    from model import compartments as cm
    G8 = cm.geometries(); s_in = np.array([x[1][:x[2]].sum() for x in G8]); rem = np.array([x[3] for x in G8]); a = cm.A_RATIO
    rows = []; U = lambda lo, hi: rng.uniform(lo, hi, draws)
    for _, r in d.dropna(subset=['zF', 'zP']).iterrows():
        g = A[(r.study, r.group)]; rep = audit.REPORTED[(r.study, r.group)]
        stored = r.donor_PVF_per_100g_ref; ref = U(*DONOR_FLOW) * (1.0 if (pd.isna(stored) or stored == indices.DEFAULT_DONOR_REF) else stored / indices.DEFAULT_DONOR_REF)
        flow = r.PVF_per_100g * np.exp(g['flow_cv_used'] / np.sqrt(rep['n_flow']) * rng.standard_normal(draws))
        if pd.notna(r.gradient_mmHg):
            grad = r.gradient_mmHg + g['gradient_sd_used'] / np.sqrt(rep['n_grad']) * rng.standard_normal(draws); cvp = np.full(draws, float(r.PVP - r.gradient_mmHg))
        else:
            cvp = U(CVPS[0], CVPS[-1]) if r.CVP_basis == 'imputed' else np.full(draws, float(r.CVP))
            grad = r.PVP + g['gradient_sd_used'] / np.sqrt(rep['n_grad']) * rng.standard_normal(draws) - cvp
        g0 = U(*NORMAL_GRADIENTS); zF = flow / ref; zR_obs = grad / g0 / zF
        i = rng.integers(0, len(G8), draws); phi = rng.uniform(0, 1, draws)
        zR_rig = (1 - rem[i]) * (s_in[i] * zF + (1 - s_in[i]) * (zF + phi * a) / (1 + a)) / zF; Q = g0 * zF * zR_rig
        cube = _bisect(lambda G: Q * ((cvp + g0 / 2) / (cvp + G / 2)) ** CUBE - G, np.zeros(draws), 2 * Q + 1) / g0 / zF
        p = U(*DOUBLING_MMHG); comp = _bisect(lambda G: Q * np.maximum(1 + (G - g0) / p, 1e-9) ** -2.0 - G, np.zeros(draws), 2 * Q + 1) / g0 / zF
        f = U(*FLOW_FIXED_SHARE); pf = (f + (1 - f) * zF * zR_rig) / zF
        q = lambda x: [round(float(v), 3) for v in np.percentile(x, [2.5, 50, 97.5])]
        models = {'rigid network': zR_rig, 'cube law': cube, 'compliance': comp, 'pressure-flow': pf}
        row = dict(study=r.study, group=r.group, zR_observed=q(zR_obs), models={})
        for name, z in models.items():
            ratio = zR_obs / z
            row['models'][name] = dict(zR=q(z), observed_over_model=q(ratio), prob_observed_below=round(float((ratio < 1).mean()), 4))
        span_lo = np.minimum(np.minimum(cube, comp), pf); span_hi = np.maximum(np.maximum(cube, comp), pf)
        row['prob_observed_below_every_relation'] = round(float((zR_obs < span_lo).mean()), 4)
        row['prob_observed_within_the_relations'] = round(float(((zR_obs >= span_lo) & (zR_obs <= span_hi)).mean()), 4)
        row['prob_observed_above_every_relation'] = round(float((zR_obs > span_hi).mean()), 4)
        rows.append(row)
    rb = [x['models']['rigid network']['prob_observed_below'] for x in rows]
    return dict(groups=rows, draws=draws, prob_below_the_rigid_network=[min(rb), max(rb)],
                ranges=dict(donor_flow=list(DONOR_FLOW), normal_gradient_mmHg=list(NORMAL_GRADIENTS), assumed_cvp_mmHg=[CVPS[0], CVPS[-1]], doubling_mmHg=list(DOUBLING_MMHG),
                            flow_fixed_share=list(FLOW_FIXED_SHARE), arterial_load='uniform between 0 and the donor arterial flow per lobule', networks=len(G8)),
                note='uniform ranges are a choice, not a measurement; the probabilities answer "if every quantity could be anywhere in its plausible range", '
                     'one group at a time. The relations are each a model: no probability is given that one of them is the mechanism')


def summary():
    d = indices.compute(pd.read_csv(f'{ROOT}/data/series.tsv', sep='\t')); B = bracket(d); sh = sinusoidal_share()
    two = [rigid_range(g['zF'], alone=False) for g in groups(d)]      # the network with its two inflows, as in P2
    ends = (float(np.mean([a for a, _ in two])), float(np.mean([b for _, b in two])))      # mean over the groups of the lowest and of the highest configuration
    req = [required(B['mean_observed'], e) for e in ends]
    return dict(
        sinusoids=dict(share_of_the_drop={f'b = {b:g}': round(v, 3) for b, v in sh.items()}, zR_floor_relative_to_rigid=[round(1 - max(sh.values()), 3), round(1 - min(sh.values()), 3)],
                       note='unlimited dilation or recruitment of the sinusoids removes at most their share of the drop'),
        outflow=dict(share_of_the_resistance=OUTLET_SHARE, zR_floor_relative_to_rigid=1 - OUTLET_SHARE, note='the outlet anastomosis of whole_graft.py'),
        rigid_with_two_inflows=dict(mean_of_the_lowest=round(ends[0], 3), mean_of_the_highest=round(ends[1], 3), range=[round(min(a for a, _ in two), 3), round(max(b for _, b in two), 3)]),
        required=dict(reduction_missing=[req[0]['reduction_missing'], req[1]['reduction_missing']],
                      whole_bed_radius_factor=[req[0]['by_share'][-1]['radius_factor'], req[1]['by_share'][-1]['radius_factor']],
                      from_the_lowest_network=req[0], from_the_highest_network=req[1],
                      note='the observed mean against the lowest and the highest configuration of the network, not against the midpoint of its range'),
        robustness=robustness(d, B), joint=joint(d), consistency=consistency(), heterogeneity=heterogeneity(d, B), wang_contrast=wang_contrast(d),
        distension=dict(relations={'cube law': f'R proportional to Pd^-{CUBE:g}, Pd the mean of the pressures on either side of the resistance',
                                   'compliance': f'blood volume linear in intrahepatic pressure, doubling over {DOUBLING_MMHG[0]:g} to {DOUBLING_MMHG[1]:g} mmHg; R proportional to volume^-2',
                                   'pressure-flow': f'gradient linear in flow, a share {FLOW_FIXED_SHARE[0]:g} to {FLOW_FIXED_SHARE[1]:g} of the control gradient independent of flow'},
                        family=FAMILY, published_check=published_check(), bracket=B, band_in_the_plane=band(),
                        flow_doubled_from_the_donor={law: dict(rise_mmHg=[round(x - indices.NORMAL_GRADIENT, 2) for x in over(law, lambda p: distensible(2.0, 1.0, 5.0, law, p)['gradient_mmHg'])]) for law in LAWS},
                        note='the relations were measured on whole livers of cats over the normal range of flow and pressure and are extrapolated; nothing is fitted to the six groups'),
        verdict='sinusoids and outlet cannot close the distance; the passive relations of the normal liver, extrapolated, span the observations of five groups, and one lies below them')


if __name__ == '__main__':
    s = summary(); os.makedirs(f'{ROOT}/results', exist_ok=True); json.dump(s, open(f'{ROOT}/results/model_distension.json', 'w'), indent=1, default=float)
    B = s['distension']['bracket']
    print(f"sinusoids carry {s['sinusoids']['share_of_the_drop']} of the drop; missing reduction {s['required']['reduction_missing']}, whole-bed radius {s['required']['whole_bed_radius_factor']}")
    for k, v in s['robustness']['scenarios'].items(): print(f"  {k:36s} within {v['within']}, below {v['below']}, keeps-distending span {v['keeps_distending_span']}")
    J = s['joint']; print('  joint: P(observed below the rigid network)', J['prob_below_the_rigid_network'])
    for r in J['groups']: print(f"    {r['study'][:13]:13s} below every relation {r['prob_observed_below_every_relation']:.3f}, within {r['prob_observed_within_the_relations']:.3f}, above {r['prob_observed_above_every_relation']:.3f}")
    for r in B['groups']: print(f"  {r['study'][:13]:13s} zF {r['zF']:.2f} observed {r['zR_observed']:.2f} {r['zR_observed_ci']}  rigid {r['zR_rigid']}  {r['zR_distensible']}  {r['position']}, in {r['relations_containing_the_value']}")
    print(f"within {B['within_the_relations']}, below {B['below_the_relations']}, above {B['above_the_relations']}; in a relation that keeps distending {B['value_in_a_relation_that_keeps_distending']}, in the one that stops {B['value_in_the_relation_that_stops_distending']}")
    print('flow doubled:', s['distension']['flow_doubled_from_the_donor'])
    for k, v in s['distension']['published_check']['states'].items(): print(f"  pig, {k}: measured {v['measured_resistance_ratio']}, " + ', '.join(f"{l} {p['resistance_ratio']}" for l, p in v['predicted'].items()))
