"""External check: published hepatectomies, the same liver measured before and after resection.

A resection leaves a fraction g of the network with the whole portal inflow, as a partial graft does, and the liver is
its own reference: flow load, pressure load and resistance load are ratios of the state after to the state before in
the same animals, so the donor references of the clinical series are not needed. The series reach flow loads the six
clinical groups do not (1.3 to 10 against 3.3 to 4.1). None of them was used to build or calibrate the model.

data     data/hepatectomy.tsv, one row per state, with the place in the paper each value comes from.
observed zF from the flow per 100 g of remnant (or the flow and the resected fraction), zP from the portocaval
         gradient: the one reported, or portal pressure minus the outflow pressure reported (central venous or hepatic
         vein); where the outflow pressure is not given, the gradient before resection is taken between 3 and 6 mmHg
         (those reported here are 3.8 to 5.8) and the outflow pressure unchanged. Intervals: the group means drawn with
         their standard errors, before and after as independent (the papers do not give the paired differences).
model    the rigid network of a remnant: the generations of the network at the flow load, the trunk at g times it
         (the trunk stays, with less flow), the artery joining at the sinusoid with the arterial flows measured; the
         three relations of distension.py on that load, with the gradient before resection as reference; where the
         haematocrit before and after is reported, the factor of viscosity.py.

    python src/model/hepatectomy.py     -> results/model_hepatectomy.json
"""
import json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from model import distension as ds, compartments as cm, viscosity

ROOT = ds.ROOT
BASELINE_GRADIENT = (3.0, 6.0)        # mmHg, where the outflow pressure is not reported
MODELS = ('rigid network', 'cube law', 'compliance', 'pressure-flow')


def load():
    return pd.read_csv(os.path.join(ROOT, 'data', 'hepatectomy.tsv'), sep='\t')


def rigid(zF, g, phi, a, geoms=None):
    """zR of a remnant that keeps a fraction g of a rigid network: range over the networks of the rigid range."""
    out = []
    for _, w, k, rem in (geoms or cm.geometries()):
        s_in = float(w[:k].sum()); out.append((1 - rem * (1 - g)) * (s_in * zF + (1 - s_in) * (zF + phi * a) / (1 + a)) / zF)
    return float(min(out)), float(max(out))


def predicted(zF, g, phi, a, g0, cvp, v=1.0, geoms=None):
    lo, hi = rigid(zF, g, phi, a, geoms); lo, hi = lo * v, hi * v; out = {'rigid network': [lo, hi]}
    for law, ps in ds.PARAMS.items():
        z = [ds.distensible(zF, zr, cvp, law, p, normal_gradient=g0)['zR'] for zr in (lo, hi) for p in ps]; out[law] = [min(z), max(z)]
    return out


def state(r, draws=20000, rng=None, geoms=None):
    rng = rng or np.random.default_rng(3); n = float(r.n); se = (lambda x: x) if r.dispersion == 'SEM' else (lambda x: x / np.sqrt(n))
    N = lambda m, s: np.maximum(m + se(s) * rng.standard_normal(draws), 1e-6); g = 1.0 - float(r.resected_fraction)
    scale = 1.0 if r.flow_basis == 'per 100 g' else 1.0 / g
    zF = float(r.flow_after / r.flow_before) * scale; zF_d = N(r.flow_after, r.flow_after_disp) / N(r.flow_before, r.flow_before_disp) * scale
    phi = float(r.HAF_after / r.HAF_before) * scale; a = float(r.HAF_before / r.flow_before)
    if pd.notna(r.gradient_before):
        g0s = [float(r.gradient_before)]; g1s = [float(r.gradient_after)]; cvps = [float(np.mean([r.out_before, r.out_after]))]; basis = 'portocaval gradient reported'
        G0 = N(r.gradient_before, r.gradient_before_disp); G1 = N(r.gradient_after, r.gradient_after_disp)
    elif pd.notna(r.out_before):
        g0s = [float(r.PVP_before - r.out_before)]; g1s = [float(r.PVP_after - r.out_after)]; cvps = [float(np.mean([r.out_before, r.out_after]))]; basis = f'portal pressure minus {r.out_basis}'
        G0 = np.maximum(N(r.PVP_before, r.PVP_before_disp) - N(r.out_before, r.out_before_disp), 1e-6); G1 = np.maximum(N(r.PVP_after, r.PVP_after_disp) - N(r.out_after, r.out_after_disp), 1e-6)
    else:
        g0s = [BASELINE_GRADIENT[0], float(np.mean(BASELINE_GRADIENT)), BASELINE_GRADIENT[1]]; cvps = [float(r.PVP_before) - x for x in g0s]; g1s = [float(r.PVP_after) - c for c in cvps]
        basis = f'outflow pressure not reported: gradient before resection {BASELINE_GRADIENT[0]:g} to {BASELINE_GRADIENT[1]:g} mmHg'
        u = rng.uniform(*BASELINE_GRADIENT, draws); out = r.PVP_before - u; G0 = np.maximum(N(r.PVP_before, r.PVP_before_disp) - out, 1e-6); G1 = np.maximum(N(r.PVP_after, r.PVP_after_disp) - out, 1e-6)
    zR_points = [g1 / g0 / zF for g0, g1 in zip(g0s, g1s)]; mid = len(g0s) // 2; zR = zR_points[mid]; zR_d = G1 / G0 / zF_d
    ci = [float(x) for x in np.percentile(zR_d, [2.5, 97.5])]
    v = viscosity.factor(float(r.hct_after), float(r.hct_before)) if pd.notna(r.hct_before) else 1.0
    regions = {m: [1e9, -1e9] for m in MODELS}
    for g0, c in zip(g0s, cvps):
        for m, (lo, hi) in predicted(zF, g, phi, a, g0, c, v, geoms).items(): regions[m] = [min(regions[m][0], lo), max(regions[m][1], hi)]
    regions['keeps distending'] = [min(regions['cube law'][0], regions['compliance'][0]), max(regions['cube law'][1], regions['compliance'][1])]
    rnd = lambda x: [round(float(y), 3) for y in x]
    return dict(study=r.study, state=r.state, n=int(r.n), shunt=r.shunt, remnant_fraction=round(g, 3), zF=round(zF, 3), zF_ci=rnd(np.percentile(zF_d, [2.5, 97.5])),
                gradient_before_mmHg=rnd(g0s), gradient_after_mmHg=rnd(g1s), gradient_basis=basis, zP=rnd([g1 / g0 for g0, g1 in zip(g0s, g1s)]),
                zR_observed=round(zR, 3), zR_observed_range=rnd([min(zR_points), max(zR_points)]), zR_observed_ci=rnd(ci), zR_observed_median=round(float(np.median(zR_d)), 3),
                arterial_load=round(phi, 3), arterial_to_portal_before=round(a, 3), viscosity_factor=round(float(v), 3),
                haematocrit=[float(r.hct_before), float(r.hct_after)] if pd.notna(r.hct_before) else None,
                model={m: rnd(x) for m, x in regions.items()}, placement={m: ds._place(zR, ci, x) for m, x in regions.items()})


def band(zFs=np.round(np.geomspace(0.25, 12.0, 60), 4), gradients=(3.0, 4.5, 6.0), outflow=(3.0, 6.0), a=0.25, geoms=None):
    """The regions of the model along the flow load for a liver that keeps its portal inflow: the remnant fraction is
    1 / zF above one, and below one the whole liver receives less (a shunt). Arterial flow per lobule unchanged; range
    over the gradients before resection and the outflow pressures of the series."""
    geoms = geoms or cm.geometries(); out = {m: [] for m in MODELS}
    for zF in zFs:
        reg = {m: [1e9, -1e9] for m in MODELS}
        for g0 in gradients:
            for c in outflow:
                for m, (lo, hi) in predicted(float(zF), min(1.0, 1.0 / float(zF)), 1.0, a, g0, c, 1.0, geoms).items(): reg[m] = [min(reg[m][0], lo), max(reg[m][1], hi)]
        for m in MODELS: out[m].append([round(reg[m][0], 4), round(reg[m][1], 4)])
    return dict(zF=[float(z) for z in zFs], **out, gradients_mmHg=list(gradients), outflow_mmHg=list(outflow), arterial_to_portal=a)


def summary(draws=20000):
    d = load(); rng = np.random.default_rng(3); geoms = cm.geometries(); rows = [state(r, draws, rng, geoms) for _, r in d.iterrows()]
    res = [x for x in rows if x['shunt'] == 'none']; zF = np.array([x['zF'] for x in res]); zR = np.array([x['zR_observed'] for x in res])
    from scipy import stats
    rho = stats.spearmanr(zF, zR); mid = lambda x, m: float(np.sqrt(x['model'][m][0] * x['model'][m][1]))
    err = {m: round(float(np.sqrt(np.mean([np.log(x['zR_observed'] / mid(x, m)) ** 2 for x in res]))), 3) for m in MODELS}
    cnt = {m: dict(point_inside=sum(x['placement'][m]['point'] == 'inside' for x in res), interval_reaches=sum(x['placement'][m]['interval'] == 'reaches' for x in res),
                   interval_below=sum(x['placement'][m]['interval'] == 'below' for x in res), interval_above=sum(x['placement'][m]['interval'] == 'above' for x in res)) for m in MODELS + ('keeps distending',)}
    low = [x for x in res if x['zF'] < 2]; high = [x for x in res if x['zF'] >= 2]
    return dict(states=rows, resections_without_shunt=len(res), flow_load_range=[float(zF.min()), float(zF.max())],
                trend=dict(spearman_rho=round(float(rho.statistic), 3), p_value=round(float(rho.pvalue), 4), n=len(res),
                           mean_zR_below_twice_the_flow=round(float(np.mean([x['zR_observed'] for x in low])), 3) if low else None,
                           mean_zR_at_twice_the_flow_or_more=round(float(np.mean([x['zR_observed'] for x in high])), 3)),
                rms_log_error_of_the_midrange=err, observed_over_midrange={m: [round(float(x['zR_observed'] / mid(x, m)), 3) for x in res] for m in MODELS}, counts=cnt, band=band(geoms=geoms), baseline_gradient_where_not_reported_mmHg=list(BASELINE_GRADIENT),
                note='group means of published series; before and after drawn as independent, which widens the intervals; states with a portocaval shunt are listed and not counted')


if __name__ == '__main__':
    s = summary(); json.dump(s, open(os.path.join(ROOT, 'results', 'model_hepatectomy.json'), 'w'), indent=1)
    for x in s['states']:
        print(f"{x['study']:13s} {x['state']:46s} zF {x['zF']:5.2f}  zR {x['zR_observed']:.2f} [{x['zR_observed_ci'][0]:.2f}, {x['zR_observed_ci'][1]:.2f}]  " + '  '.join(f"{m} {x['model'][m][0]:.2f}-{x['model'][m][1]:.2f}" for m in MODELS))
    print('trend', s['trend']); print('rms log error', s['rms_log_error_of_the_midrange']); print('counts', s['counts'])
