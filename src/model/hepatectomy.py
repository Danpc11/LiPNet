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
         their standard errors. The papers do not give the paired differences, so the correlation between before and
         after is a scenario (0, 0.5, 0.8; the main intervals take 0), and the draw of the state before is shared by
         the states that use the same measurement before (column baseline_group). Which animals a state was measured
         in is a column of the table (animal_group), not inferred from the values.
trend    seven states from five groups of animals: rho of the seven is a description with no P value; the test is on
         one state per group of animals, with the exact permutation P.
errors   of every admissible configuration of a model, and to its range; the mid-range of a model is not a prediction.
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
PRE_POST = (0.0, 0.5, 0.8)          # correlation between the group means before and after, where the papers give no paired difference


def load():
    return pd.read_csv(os.path.join(ROOT, 'data', 'hepatectomy.tsv'), sep='\t')


def rigid_each(zF, g, phi, a, geoms=None):
    """zR of a remnant that keeps a fraction g of a rigid network, for each network of the rigid range: the trunk stays
    and carries g zF times its flow, the other generations zF, and the artery joins at the sinusoid."""
    out = []
    for _, w, k, rem in (geoms or cm.geometries()):
        s_in = float(w[:k].sum()); out.append(float((1 - rem * (1 - g)) * (s_in * zF + (1 - s_in) * (zF + phi * a) / (1 + a)) / zF))
    return out


def rigid(zF, g, phi, a, geoms=None):
    """Range of `rigid_each` over the networks."""
    z = rigid_each(zF, g, phi, a, geoms); return min(z), max(z)


def predicted(zF, g, phi, a, g0, cvp, v=1.0, geoms=None):
    lo, hi = rigid(zF, g, phi, a, geoms); lo, hi = lo * v, hi * v; out = {'rigid network': [lo, hi]}
    for law, ps in ds.PARAMS.items():
        z = [ds.distensible(zF, zr, cvp, law, p, normal_gradient=g0)['zR'] for zr in (lo, hi) for p in ps]; out[law] = [min(z), max(z)]
    return out


def group_of(r):
    """The animals a state was measured in, as the table says (`animal_group`): the states of a stepwise resection are
    the same animals. It is not inferred from the numbers: two groups of animals can share the state before."""
    return str(r.animal_group)


def baseline_of(r):
    """The measurement before the procedure that a state uses (`baseline_group`): the same animals in a stepwise
    resection, and one pooled value for the three groups of Golriz 2022. Its draws are shared by those states."""
    return str(r.baseline_group)


def baseline(rng, draws):
    """Standard draws for the state before resection of one group of animals, shared by every state measured in them."""
    return dict(flow=rng.standard_normal(draws), PVP=rng.standard_normal(draws), out=rng.standard_normal(draws), grad=rng.standard_normal(draws), u=rng.uniform(0, 1, draws))


def configurations(zF, g, phi, a, g0, cvp, v=1.0, geoms=None):
    """zR of each model in every admissible configuration, in a fixed order: network by network and, for the relations,
    parameter by parameter. The range of `predicted` is the span of these; none of them is a central prediction."""
    out = {m: [] for m in MODELS}
    for z in rigid_each(zF, g, phi, a, geoms):
        zr = z * v; out['rigid network'].append(float(zr))
        for law, ps in ds.PARAMS.items(): out[law] += [float(ds.distensible(zF, zr, cvp, law, p, normal_gradient=g0)['zR']) for p in ps]
    return out


def state(r, draws=20000, rng=None, geoms=None, corr=0.0, base=None, private=False):
    """One state against the model. `corr` is the correlation between the group means before and after (0: independent,
    which is the widest interval for a positive correlation); `base` holds the draws of the state before resection, shared
    by the states measured in the same animals. With `private`, also the draws and the configurations of each model."""
    rng = rng or np.random.default_rng(3); base = base or baseline(rng, draws); n = float(r.n); se = (lambda x: x) if r.dispersion == 'SEM' else (lambda x: x / np.sqrt(n)); w = np.sqrt(1 - corr ** 2)
    Nb = lambda m, s, key: np.maximum(m + se(s) * base[key], 1e-6); Na = lambda m, s, key: np.maximum(m + se(s) * (corr * base[key] + w * rng.standard_normal(draws)), 1e-6); g = 1.0 - float(r.resected_fraction)
    scale = 1.0 if r.flow_basis == 'per 100 g' else 1.0 / g
    zF = float(r.flow_after / r.flow_before) * scale; zF_d = Na(r.flow_after, r.flow_after_disp, 'flow') / Nb(r.flow_before, r.flow_before_disp, 'flow') * scale
    phi = float(r.HAF_after / r.HAF_before) * scale; a = float(r.HAF_before / r.flow_before)
    if pd.notna(r.gradient_before):
        g0s = [float(r.gradient_before)]; g1s = [float(r.gradient_after)]; cvps = [float(np.mean([r.out_before, r.out_after]))]; basis = 'portocaval gradient reported'
        G0 = Nb(r.gradient_before, r.gradient_before_disp, 'grad'); G1 = Na(r.gradient_after, r.gradient_after_disp, 'grad')
    elif pd.notna(r.out_before):
        g0s = [float(r.PVP_before - r.out_before)]; g1s = [float(r.PVP_after - r.out_after)]; cvps = [float(np.mean([r.out_before, r.out_after]))]; basis = f'portal pressure minus {r.out_basis}'
        G0 = np.maximum(Nb(r.PVP_before, r.PVP_before_disp, 'PVP') - Nb(r.out_before, r.out_before_disp, 'out'), 1e-6); G1 = np.maximum(Na(r.PVP_after, r.PVP_after_disp, 'PVP') - Na(r.out_after, r.out_after_disp, 'out'), 1e-6)
    else:
        g0s = [BASELINE_GRADIENT[0], float(np.mean(BASELINE_GRADIENT)), BASELINE_GRADIENT[1]]; cvps = [float(r.PVP_before) - x for x in g0s]; g1s = [float(r.PVP_after) - c for c in cvps]
        basis = f'outflow pressure not reported: gradient before resection {BASELINE_GRADIENT[0]:g} to {BASELINE_GRADIENT[1]:g} mmHg'
        u = BASELINE_GRADIENT[0] + (BASELINE_GRADIENT[1] - BASELINE_GRADIENT[0]) * base['u']; out = r.PVP_before - u
        G0 = np.maximum(Nb(r.PVP_before, r.PVP_before_disp, 'PVP') - out, 1e-6); G1 = np.maximum(Na(r.PVP_after, r.PVP_after_disp, 'PVP') - out, 1e-6)
    zR_points = [g1 / g0 / zF for g0, g1 in zip(g0s, g1s)]; mid = len(g0s) // 2; zR = zR_points[mid]; zR_d = G1 / G0 / zF_d
    ci = [float(x) for x in np.percentile(zR_d, [2.5, 97.5])]
    v = viscosity.factor(float(r.hct_after), float(r.hct_before)) if pd.notna(r.hct_before) else 1.0
    conf = [configurations(zF, g, phi, a, g0, c, v, geoms) for g0, c in zip(g0s, cvps)]
    regions = {m: [min(min(k[m]) for k in conf), max(max(k[m]) for k in conf)] for m in MODELS}
    regions['keeps distending'] = [min(regions['cube law'][0], regions['compliance'][0]), max(regions['cube law'][1], regions['compliance'][1])]
    rnd = lambda x: [round(float(y), 3) for y in x]
    pub = dict(study=r.study, state=r.state, n=int(r.n), shunt=r.shunt, remnant_fraction=round(g, 3), zF=round(zF, 3), zF_ci=rnd(np.percentile(zF_d, [2.5, 97.5])),
               gradient_before_mmHg=rnd(g0s), gradient_after_mmHg=rnd(g1s), gradient_basis=basis, zP=rnd([g1 / g0 for g0, g1 in zip(g0s, g1s)]),
               zR_observed=round(zR, 3), zR_observed_range=rnd([min(zR_points), max(zR_points)]), zR_observed_ci=rnd(ci), zR_observed_median=round(float(np.median(zR_d)), 3),
               arterial_load=round(phi, 3), arterial_to_portal_before=round(a, 3), viscosity_factor=round(float(v), 3),
               haematocrit=[float(r.hct_before), float(r.hct_after)] if pd.notna(r.hct_before) else None,
               model={m: rnd(x) for m, x in regions.items()}, observed_over_model={m: rnd([zR / x[1], zR / x[0]]) for m, x in regions.items()},
               placement={m: ds._place(zR, ci, x) for m, x in regions.items()})
    return (pub, dict(zF_d=zF_d, zR_d=zR_d, conf=conf, zR_points=zR_points, regions=regions, zR=zR)) if private else pub


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


def _rho(x, y):
    from scipy import stats
    return float(stats.spearmanr(x, y).statistic)


def exact_p(x, y):
    """Spearman rho with its exact permutation P values (two-sided, and one-sided for a fall): every ordering of y."""
    from itertools import permutations
    r = _rho(x, y); all_r = np.array([_rho(x, q) for q in permutations(y)])
    return dict(rho=round(r, 3), p_two_sided=round(float(np.mean(np.abs(all_r) >= abs(r) - 1e-12)), 4), p_one_sided=round(float(np.mean(all_r <= r + 1e-12)), 4), n=len(x))


def _rho_draws(zF_d, zR_d):
    """Spearman rho in every draw: columns are states."""
    from scipy import stats
    a = stats.rankdata(zF_d, axis=1); b = stats.rankdata(zR_d, axis=1); a -= a.mean(1, keepdims=True); b -= b.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


def run(d, draws, geoms, corr):
    """Every state at one correlation between before and after; the states that use one measurement before share its draws."""
    rng = np.random.default_rng(3); bases = {}; out = []
    for _, r in d.iterrows():
        key = baseline_of(r)
        if key not in bases: bases[key] = baseline(rng, draws)
        out.append(state(r, draws, rng, geoms, corr=corr, base=bases[key], private=True))
    return out


def trend(d, runs):
    """How the resistance load moves with the flow load across the resections. Seven states come from five groups of
    animals, so no P value is given for the seven: rho is a description there. The test is on one state per group of
    animals, for every choice of that state, with the exact permutation P. The draws carry the sampling error of the
    group means (state before shared inside a group) into rho."""
    from itertools import product
    keep = [i for i, (_, r) in enumerate(d.iterrows()) if r.shunt == 'none']; rows = [r for _, r in d.iterrows()]; pub = [runs[0.0][i][0] for i in keep]
    zF = np.array([x['zF'] for x in pub]); zR = np.array([x['zR_observed'] for x in pub]); groups = {}
    for j, i in enumerate(keep): groups.setdefault(group_of(rows[i]), []).append(j)
    sel = []
    for choice in product(*groups.values()):
        c = list(choice); e = exact_p(zF[c], zR[c]); e['states'] = [f"{pub[j]['study']}, {pub[j]['state']}" for j in c]; e['index'] = c; sel.append(e)
    mc = {}
    for corr, run_ in runs.items():
        F = np.column_stack([run_[i][1]['zF_d'] for i in keep]); Z = np.column_stack([run_[i][1]['zR_d'] for i in keep]); r_all = _rho_draws(F, Z)
        r_sel = np.concatenate([_rho_draws(F[:, e['index']], Z[:, e['index']]) for e in sel])
        q = lambda x: [round(float(y), 3) for y in np.percentile(x, [2.5, 50, 97.5])]
        mc[f'{corr:g}'] = dict(all_states=dict(rho=q(r_all), share_of_draws_negative=round(float((r_all < 0).mean()), 4)),
                               one_state_per_group=dict(rho=q(r_sel), share_of_draws_negative=round(float((r_sel < 0).mean()), 4)))
    low = [x for x in pub if x['zF'] < 2]; high = [x for x in pub if x['zF'] >= 2]
    for e in sel: del e['index']
    return dict(spearman_rho=round(_rho(zF, zR), 3), n=len(pub), groups_of_animals=len(groups), p_value=None,
                one_state_per_group=sel, rho_one_state_per_group=[min(e['rho'] for e in sel), max(e['rho'] for e in sel)],
                p_two_sided_one_state_per_group=[min(e['p_two_sided'] for e in sel), max(e['p_two_sided'] for e in sel)],
                p_one_sided_one_state_per_group=[min(e['p_one_sided'] for e in sel), max(e['p_one_sided'] for e in sel)],
                draws=mc, mean_zR_below_twice_the_flow=round(float(np.mean([x['zR_observed'] for x in low])), 3) if low else None,
                mean_zR_at_twice_the_flow_or_more=round(float(np.mean([x['zR_observed'] for x in high])), 3),
                note='no P value for the seven states: three of them are the same animals. rho of the seven is a description')


def errors(priv):
    """How far the observations are from each model, without a central prediction. by_configuration: the root-mean-square
    logarithmic error of every admissible configuration of the model (network, parameter and, where the outflow pressure
    is not reported, the gradient assumed before resection), as its smallest and largest value. to_the_range: the same
    error with the distance to the range of the model, zero for an observation inside it."""
    by, rng_ = {}, {}
    for m in MODELS:
        e = []
        for k in range(3):
            for c in range(len(priv[0]['conf'][0][m])):
                e.append(np.sqrt(np.mean([np.log(p['zR_points'][min(k, len(p['zR_points']) - 1)] / p['conf'][min(k, len(p['conf']) - 1)][m][c]) ** 2 for p in priv])))
        by[m] = [round(float(min(e)), 3), round(float(max(e)), 3)]
    for m in MODELS + ('keeps distending',):
        dist = [0.0 if p['regions'][m][0] <= p['zR'] <= p['regions'][m][1] else float(np.log(p['zR'] / p['regions'][m][1]) if p['zR'] > p['regions'][m][1] else np.log(p['regions'][m][0] / p['zR'])) for p in priv]
        rng_[m] = round(float(np.sqrt(np.mean(np.square(dist)))), 3)
    return dict(by_configuration=by, to_the_range=rng_, note='the mid-range of a model is not a prediction and is not used')


def summary(draws=20000):
    d = load(); geoms = cm.geometries(); runs = {c: run(d, draws, geoms, c) for c in PRE_POST}; rows = [x[0] for x in runs[0.0]]
    idx = [i for i, x in enumerate(rows) if x['shunt'] == 'none']; res = [rows[i] for i in idx]; zF = np.array([x['zF'] for x in res])
    count = lambda rs: {m: dict(point_inside=sum(x['placement'][m]['point'] == 'inside' for x in rs), interval_reaches=sum(x['placement'][m]['interval'] == 'reaches' for x in rs),
                                interval_below=sum(x['placement'][m]['interval'] == 'below' for x in rs), interval_above=sum(x['placement'][m]['interval'] == 'above' for x in rs)) for m in MODELS + ('keeps distending',)}
    cov = {f'{c:g}': dict(zR_observed_ci=[runs[c][i][0]['zR_observed_ci'] for i in range(len(rows))], counts=count([runs[c][i][0] for i in idx]),
                          widest_interval_ratio=round(float(max(runs[c][i][0]['zR_observed_ci'][1] / runs[c][i][0]['zR_observed_ci'][0] for i in idx)), 2)) for c in PRE_POST}
    changes = [dict(correlation=c, state=f"{rows[i]['study']}, {rows[i]['state']}", model=m, at_zero=rows[i]['placement'][m]['interval'], now=runs[c][i][0]['placement'][m]['interval'])
               for c in PRE_POST[1:] for i in idx for m in MODELS + ('keeps distending',) if runs[c][i][0]['placement'][m]['interval'] != rows[i]['placement'][m]['interval']]
    fam = lambda c: {m: cov[f'{c:g}']['counts'][m] for m in ('rigid network', 'keeps distending', 'pressure-flow')}
    animals = int(d[d.shunt == 'none'].drop_duplicates('animal_group').n.sum())                         # the states of a stepwise resection are the same animals
    return dict(states=rows, resections_without_shunt=len(res), animals_resected=animals, flow_load_range=[float(zF.min()), float(zF.max())],
                trend=trend(d, runs), errors=errors([runs[0.0][i][1] for i in idx]), counts=count(res),
                covariance_sensitivity=dict(correlations=list(PRE_POST), by_correlation=cov, placements_that_change=changes, unchanged_against_the_rigid_network_and_the_two_families=bool(all(fam(c) == fam(0.0) for c in PRE_POST)),
                                            note='correlation between the group means before and after; the main intervals take 0. The state before is shared by the states of the same animals'),
                band=band(geoms=geoms), baseline_gradient_where_not_reported_mmHg=list(BASELINE_GRADIENT),
                note='group means of published series; states with a portocaval shunt are listed and not counted')


if __name__ == '__main__':
    s = summary(); json.dump(s, open(os.path.join(ROOT, 'results', 'model_hepatectomy.json'), 'w'), indent=1)
    for x in s['states']:
        print(f"{x['study']:13s} {x['state']:46s} zF {x['zF']:5.2f}  zR {x['zR_observed']:.2f} [{x['zR_observed_ci'][0]:.2f}, {x['zR_observed_ci'][1]:.2f}]  " + '  '.join(f"{m} {x['model'][m][0]:.2f}-{x['model'][m][1]:.2f}" for m in MODELS))
    T = s['trend']; print('rho of the seven states', T['spearman_rho'], '| one state per group of animals: rho', T['rho_one_state_per_group'], 'exact two-sided P', T['p_two_sided_one_state_per_group'], '| draws', T['draws']['0'])
    print('errors', s['errors']); print('covariance', {k: (v['counts']['rigid network'], v['counts']['keeps distending'], v['widest_interval_ratio']) for k, v in s['covariance_sensitivity']['by_correlation'].items()}, s['covariance_sensitivity']['placements_that_change'])
