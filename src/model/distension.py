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
distension    The venous resistance of the normal liver is passively distensible. Two published relations, both
              measured on whole livers of cats and neither on a graft:
                cube law     resistance proportional to 1 / Pd^3 (index of contractility, IC = R Pd^3), the distending
                             pressure Pd being the mean of the pressures above and below the resistance site (Lautt,
                             Greenway and Legare, Am J Physiol 1991;260:G325; Lautt and Legare, Am J Physiol
                             1992;263:G702; stated in these words by Lautt and Ming, in Portal Hypertension, Humana
                             2005, p. 88, and by Ingles, Legare and Lautt, Can J Physiol Pharmacol 1994;72:687);
                compliance   hepatic blood volume doubles for about 8 mmHg of intrahepatic pressure, which the portal
                             venous pressure estimates, linearly over the physiological range (Lautt and Greenway,
                             Am J Physiol 1976;231:292; Greenway, Seaman and Innes, Am J Physiol 1985;248:H468; both
                             as summarised by Lautt, Hepatol Res 2007;37:891). With the length of the vessels fixed
                             the radius goes as the square root of the volume, and with Poiseuille flow the
                             resistance as volume^-2. This takes the resistance vessels to distend as the capacitance
                             vessels do, which is an assumption and an upper bound.
              The sources give the size of the effect in a normal liver: a doubling of portal flow raises portal
              pressure by about 2 mmHg (Lautt and Ming 2005). They also give its limit: in the normal liver the changes
              of resistance are notable below 10 mmHg and minor above (Ingles et al. 1994), and a graft is at 14 to 22
              mmHg. So the fully distensible bed is a lower bracket, not a prediction.
check         `published_check` sets the two relations against pig livers they were not derived from (Ricciardi et al.,
              J Gastrointest Surg 2001;5:490): in situ and isolated, and the same livers after 2 h of cold ischaemia,
              when resistance was twice as high at every flow.

For a group with flow load zF and central venous pressure CVP the distensible bed is solved self-consistently: the
gradient G satisfies G = 5 zF zR_rigid d(G), with d the relation above and the donor at a gradient of 5 mmHg. Its
inputs are the flow load, the rigid network and the central venous pressure; the portal pressure of the recipient is
the output and is not used. The rigid value runs from the network with the two inflows (0.81) to the network alone
(1.00), because the published relations were measured with the hepatic artery in place and may already contain its part.

What comes out is a bracket: rigid network above, fully distensible bed below. Each group is placed with its 95% sampling
interval: below the bed, at the bed (the interval reaches it) or between the two. `distensible_share` is the fraction of
the bed that would have to follow the published relation, at the pressure the series reports, to give the observed value.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import indices
from model import liver, human

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
DOUBLING_MMHG = 8.0                      # distending pressure that doubles hepatic blood volume (Lautt and Greenway 1976)
CUBE = 3.0                               # exponent of the index of contractility, IC = R Pd^3 (Lautt, Greenway and Legare 1991)
CVPS = (3.0, 5.0, 7.0, 9.0)              # central venous pressure where a series reports none, as in predictions.REF_CVPS
OUTLET_SHARE = 0.05                      # share of the resistance held by the outlet anastomosis (whole_graft.py, r_out)

# The two published relations. Each takes the gradient G across the liver, the central venous pressure and the gradient
# g0 of the reference state, and returns resistance relative to that state.
#   cube law     the distending pressure is the mean of the pressures above and below the resistance site, CVP + G / 2
#   compliance   the distending pressure is the intrahepatic pressure, which the portal venous pressure estimates, CVP + G
LAWS = {
    'cube law': lambda G, cvp, g0: ((cvp + g0 / 2.0) / (cvp + G / 2.0)) ** CUBE,
    'compliance': lambda G, cvp, g0: (1.0 + (G - g0) / DOUBLING_MMHG) ** -2.0,
}


def distensible(zF, zR_rigid, cvp, law, normal_gradient=indices.NORMAL_GRADIENT):
    """Gradient and loads of a bed that follows `law`, at flow load zF: the root of G = g0 zF zR_rigid d(G)."""
    d = LAWS[law]; lo, hi = 0.0, normal_gradient * zF * zR_rigid * 2 + 1
    f = lambda G: normal_gradient * zF * zR_rigid * d(G, cvp, normal_gradient) - G        # decreasing in G above the reference gradient
    for _ in range(80):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    G = 0.5 * (lo + hi); return dict(gradient_mmHg=float(G), zP=float(G / normal_gradient), zR=float(G / normal_gradient / zF), factor=float(d(G, cvp, normal_gradient)))


# Pressure and flow of normal livers of a third species, not used to derive either relation: the pig livers of Ricciardi
# et al. (J Gastrointest Surg 2001;5:490), in situ and isolated on a bypass circuit, and the same livers after 2 h of
# cold ischaemia. Portal flow in mL/min, portal venous pressure in mmHg, resistance (portal minus caval pressure over
# portal flow) in mmHg min/mL, all as printed in the Results of that paper.
RICCIARDI = {
    'in situ': dict(n=9, flow=(400.0, 1100.0), pvp=(6.8, 9.6), hvr=(0.016, 0.0086)),
    'isolated, perfused': dict(n=7, flow=(300.0, 800.0), pvp=(5.14, 8.14), hvr=(0.017, 0.010)),
    'after 2 h of cold ischaemia': dict(n=6, flow=(300.0, 800.0), pvp=(10.2, 15.8), hvr=(0.034, 0.020)),
}


def published_check():
    """The two relations against measurements they were not derived from. For each state of the pig liver: the flow is
    raised by the published factor from the published starting gradient and caval pressure, and the resistance the
    relation gives at the higher flow, relative to the lower, is set against the measured one."""
    out = {}
    for name, r in RICCIARDI.items():
        g0 = r['hvr'][0] * r['flow'][0]; cvp = r['pvp'][0] - g0; zF = r['flow'][1] / r['flow'][0]; pred = {}
        for law in LAWS:
            x = distensible(zF, 1.0, cvp, law, normal_gradient=g0); pred[law] = dict(resistance_ratio=round(x['zR'], 3), portal_pressure_mmHg=round(cvp + x['gradient_mmHg'], 2))
        out[name] = dict(n=r['n'], flow_factor=round(zF, 3), starting_gradient_mmHg=round(g0, 2), caval_pressure_mmHg=round(cvp, 2), measured_resistance_ratio=round(r['hvr'][1] / r['hvr'][0], 3),
                         measured_portal_pressure_mmHg=r['pvp'][1], predicted=pred, resistance_relative_to_the_isolated_liver_at_low_flow=round(r['hvr'][1] / RICCIARDI['isolated, perfused']['hvr'][0], 3))
    return out


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


def bracket(d):
    import audit
    CI = {(x['study'], x['group']): x['zR_ci_sampling'] for x in audit.paired(d, correlations=())['groups']}; rows = []
    for g in groups(d):
        lo, hi = rigid_range(g['zF']); cvps = (g['cvp'],) if g['cvp'] is not None else CVPS; sol = {}
        for law in LAWS:
            z = [distensible(g['zF'], zr, c, law)['zR'] for zr in (lo, hi) for c in cvps]; sol[law] = [round(min(z), 3), round(max(z), 3)]
        low = min(v[0] for v in sol.values()); low_hi = max(v[1] for v in sol.values())
        # the share of the bed that must follow each relation, at the portal pressure the series reports, to give the observed zR
        share = {}
        for law, dlaw in LAWS.items():
            s = []
            for zr in (lo, hi):
                for c in cvps:
                    dd = dlaw(indices.NORMAL_GRADIENT * g['zP'], c, indices.NORMAL_GRADIENT)
                    s.append((1 - g['zR'] / zr) / (1 - dd) if dd < 1 else np.inf)
            share[law] = [round(float(min(s)), 2), round(float(max(s)), 2)]
        # where the observation sits, with its 95% sampling interval: below the bed, at the bed (the interval reaches it), or between bed and rigid network
        ci = CI[(g['study'], g['group'])]
        pos = 'below the distensible bed' if ci[1] < low else ('at the distensible bed' if ci[0] <= low_hi else ('between the two' if g['zR'] < lo else 'at the rigid network'))
        rows.append(dict(study=g['study'], group=g['group'], zF=round(g['zF'], 3), zR_observed=round(g['zR'], 3), cvp_mmHg=g['cvp'], cvp_basis=g['cvp_basis'],
                         zR_rigid=[round(lo, 3), round(hi, 3)], zR_distensible=sol, zR_distensible_all=[round(low, 3), round(low_hi, 3)],
                         zR_observed_ci=[round(float(ci[0]), 3), round(float(ci[1]), 3)], distensible_share=share, position=pos,
                         share_of_the_departure_reached=[round((1 - lo) / (1 - g['zR']), 2), round((1 - low) / (1 - g['zR']), 2)]))
    obs = np.array([r['zR_observed'] for r in rows])
    return dict(groups=rows, mean_observed=round(float(obs.mean()), 3),
                mean_rigid=[round(float(np.mean([r['zR_rigid'][k] for r in rows])), 3) for k in (0, 1)],
                mean_distensible={law: [round(float(np.mean([r['zR_distensible'][law][k] for r in rows])), 3) for k in (0, 1)] for law in LAWS},
                mean_distensible_all=[round(float(np.mean([r['zR_distensible_all'][k] for r in rows])), 3) for k in (0, 1)],
                between_the_two=int(sum(r['position'] == 'between the two' for r in rows)), at_the_distensible_bed=int(sum(r['position'] == 'at the distensible bed' for r in rows)),
                below_the_distensible_bed=[f"{r['study']} / {r['group']}" for r in rows if r['position'] == 'below the distensible bed'],
                above_the_bed_as_point_values=int(sum(r['zR_observed'] > r['zR_distensible_all'][1] for r in rows)))


def band(zFs=np.round(np.arange(1.0, 5.01, 0.1), 2)):
    """The fully distensible bed in the flow-pressure plane: lowest and highest pressure load at each flow load, over the
    two relations, the rigid range and a central venous pressure of 3 to 9 mmHg."""
    lo, hi = [], []
    for zF in zFs:
        r = rigid_range(float(zF)); z = [distensible(float(zF), zr, c, law)['zP'] for law in LAWS for zr in r for c in CVPS]; lo.append(min(z)); hi.append(max(z))
    return dict(zF=[float(x) for x in zFs], zP_low=[round(float(x), 4) for x in lo], zP_high=[round(float(x), 4) for x in hi])


def required(zR_obs=0.584, zR_rigid=0.853):
    """With no relation assumed: the widening that closes the distance if a share s of the resistance widens uniformly
    (Poiseuille, R as radius^-4). None where even an unlimited widening of that share is not enough."""
    target = zR_obs / zR_rigid; out = []
    for s in (0.15, 0.32, 0.46, 0.5, 0.75, 1.0):
        x = (target - (1 - s)) / s
        out.append(dict(share_that_widens=s, resistance_of_that_share=None if x <= 0 else round(float(x), 3), radius_factor=None if x <= 0 else round(float(x ** -0.25), 3)))
    return dict(ratio_observed_to_rigid=round(float(target), 3), reduction_missing=round(float(1 - target), 3), smallest_share_that_can_close_it=round(float(1 - target), 3), by_share=out)


def summary():
    d = indices.compute(pd.read_csv(f'{ROOT}/data/series.tsv', sep='\t')); B = bracket(d); sh = sinusoidal_share()
    two = [rigid_range(g['zF'], alone=False) for g in groups(d)]; mid = float(np.mean([0.5 * (a + b) for a, b in two]))      # the network with its two inflows, as in P2
    return dict(
        sinusoids=dict(share_of_the_drop={f'b = {b:g}': round(v, 3) for b, v in sh.items()}, zR_floor_relative_to_rigid=[round(1 - max(sh.values()), 3), round(1 - min(sh.values()), 3)],
                       note='unlimited dilation or recruitment of the sinusoids removes at most their share of the drop'),
        outflow=dict(share_of_the_resistance=OUTLET_SHARE, zR_floor_relative_to_rigid=1 - OUTLET_SHARE, note='the outlet anastomosis of whole_graft.py'),
        rigid_with_two_inflows=dict(mean=round(mid, 3), range=[round(min(a for a, _ in two), 3), round(max(b for _, b in two), 3)]),
        required=required(B['mean_observed'], mid),
        distension=dict(relations={'cube law': f'R proportional to Pd^-{CUBE:g}, Pd the mean of the pressures above and below the resistance site', 'compliance': f'blood volume doubles per {DOUBLING_MMHG:g} mmHg of intrahepatic pressure, R proportional to volume^-2'},
                        published_check=published_check(),
                        bracket=B, band_in_the_plane=band(), flow_doubled_from_the_donor={law: dict(gradient_mmHg=round(distensible(2.0, 1.0, 5.0, law)['gradient_mmHg'], 2), rise_mmHg=round(distensible(2.0, 1.0, 5.0, law)['gradient_mmHg'] - indices.NORMAL_GRADIENT, 2)) for law in LAWS},
                        note='both relations were measured on whole livers over the normal range of portal pressure and are extrapolated; nothing is fitted to the six groups'),
        verdict='sinusoids and outlet cannot close the distance; a passively distensible bed can: the observations run from the fully distensible bed to short of the rigid network, and one lies below both')


if __name__ == '__main__':
    s = summary(); os.makedirs(f'{ROOT}/results', exist_ok=True); json.dump(s, open(f'{ROOT}/results/model_distension.json', 'w'), indent=1, default=float)
    B = s['distension']['bracket']
    print(f"sinusoids carry {s['sinusoids']['share_of_the_drop']} of the drop: zR floor {s['sinusoids']['zR_floor_relative_to_rigid']} of the rigid value; outlet {s['outflow']['zR_floor_relative_to_rigid']}")
    print(f"missing reduction {s['required']['reduction_missing']}: " + '; '.join(f"share {x['share_that_widens']}: radius x{x['radius_factor']}" for x in s['required']['by_share']))
    for r in B['groups']: print(f"  {r['study'][:13]:13s} zF {r['zF']:.2f} observed {r['zR_observed']:.2f}  rigid {r['zR_rigid']}  distensible {r['zR_distensible']}  share needed {r['distensible_share']}  {r['position']}")
    print(f"means: observed {B['mean_observed']}, rigid {B['mean_rigid']}, distensible {B['mean_distensible']}; between the two {B['between_the_two']}, at the bed {B['at_the_distensible_bed']}, below {len(B['below_the_distensible_bed'])}; flow doubled: {s['distension']['flow_doubled_from_the_donor']}")
    for k, v in s['distension']['published_check'].items(): print(f"  pig liver, {k}: flow x{v['flow_factor']}, resistance measured {v['measured_resistance_ratio']}, relations {[p['resistance_ratio'] for p in v['predicted'].values()]}; relative to the isolated liver at low flow {v['resistance_relative_to_the_isolated_liver_at_low_flow']}")
