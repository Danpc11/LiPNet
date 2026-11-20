"""Test the predictions of the perfusion-network model against independent data.

    python src/predictions.py        -> results/predictions.json, printed summary

The model (src/model/) describes the liver as a dissipation-optimal perfusion network at fixed vascular maintenance
cost. It fixes the scale of the design tau0 = sqrt(D*/C) (a wall shear stress only at b = 1), the scaling of the optimum with the number of lobules,
and the local relation dP = f * R, from which the two indices follow. It says nothing about recipient severity,
donor age or steatosis. Four predictions follow, and each is tested here against data the model was not fitted to:

  P1  Comparative anatomy. With portal pressure invariant across mammals, the scale of the design tau0 must be nearly
      invariant with body mass and the vascular mass must scale close to the observed exponent of hepatic blood volume (0.86).
  P2  zP = zF * zR is an identity and proves nothing. The prediction is the value of zR: the network gives a partial
      graft 0.8 to 1 from its geometry and the two inflows, with no pressure of a recipient in it, and that number is
      compared with the zR each cohort reports. Whether zR differs between surgical configurations is a third,
      separate question, asked of the mechanistic classes of the table.
  P3  Within a cohort, the outcome increases with zP, with a common slope and a centre-specific level. The contrast is
      between groups of one cohort: an ecological association, not the effect of pressure on a patient.
  P4  Exploratory comparison of haemodynamic thresholds across liver diseases and surgical settings. The thresholds
      are each about 10 mmHg, so they sit at zP = 2 by the choice of a 5 mmHg reference: arithmetic, not a prediction.
      HVPG and the portocaval gradient of a graft are different measurements and no equivalence has been validated.

Cached fits are reused only while the table and the model code are the ones that produced them (results/bayes_inputs.json);
LIPNET_REBUILD=1 ignores every cache.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import indices
from model import scaling as sc

ROOT = os.path.join(os.path.dirname(__file__), '..'); OUT = f'{ROOT}/results'


def p1_allometry():
    """Exponents of tau0 and of the vascular mass under an invariant portal pressure, for the measured flow and lobule
    scalings (Kruepunga 2019: portal flow ~ M^0.77, portocentral distance ~ M^0.10, portal pressure 6-11 mmHg)."""
    gs = np.arange(8, 19); rows = []
    for f, flab in ((0.77, 'portal flow ~ M^0.77'), (0.88, 'total hepatic flow ~ M^0.88')):
        for l, n in ((0.0, 0.89), (0.10, 0.59)):
            for b in (1.0, 0.9, 0.85, 0.8, 0.75, 2 / 3, 0.5):
                R = [sc.canopy_optimum(g, b, 3) for g in gs]; N = np.array([r['N'] for r in R])
                eD = np.polyfit(np.log(N), np.log([r['D'] for r in R]), 1)[0]
                eV = np.polyfit(np.log(N), np.log([r['V'] for r in R]), 1)[0]
                c = b * (f + eD * n + 3 * l) / 2
                rows.append(dict(flow=flab, lobule_exponent=l, b=b,
                                 tau0_exponent=(2 * f - (2 + b) * c / b + eD * n + 3 * l) / 2,
                                 vascular_mass_exponent=c / b + eV * n))
    t = pd.DataFrame(rows); t.to_csv(f'{OUT}/prediction_allometry.tsv', sep='\t', index=False, float_format='%.4f')
    # reported over the whole prespecified range of b, with no value selected for agreeing with the observation
    obs_lo, obs_hi = 0.80, 0.92          # hepatic blood volume exponent with its reported spread
    return dict(b_range=[float(t.b.min()), float(t.b.max())],
                tau0_exponent_range=[float(t.tau0_exponent.min()), float(t.tau0_exponent.max())],
                tau0_exponent_iqr=[float(t.tau0_exponent.quantile(.25)), float(t.tau0_exponent.quantile(.75))],
                tau0_invariant=bool(t.tau0_exponent.abs().max() < 0.1),
                vascular_mass_exponent_range=[float(t.vascular_mass_exponent.min()), float(t.vascular_mass_exponent.max())],
                observed_vascular_mass_exponent=0.86, observed_range=[obs_lo, obs_hi],
                fraction_compatible=float(((t.vascular_mass_exponent >= obs_lo) & (t.vascular_mass_exponent <= obs_hi)).mean()),
                note='no value of b was selected for fitting the observation; the whole prespecified range is reported')


def p2_discordance(d):
    """What the six groups with both loads say, in three levels that are not to be mixed.

    identity      zR = zP / zF by definition: it cannot be evidence of a mechanism and is not tested.
    description   zR of each group with the interval of its ratio (audit.paired: the dispersions each paper reports,
                  the largest reported one where none is given, and the interval repeated for a correlation between
                  the sampling errors of the two means). The unit of the test is the cohort: two groups of Wang 2014
                  are one cohort. The sign test is exact and its direction is the one the network gives (zR <= 1).
    configuration the mechanistic class of each group (table column `configuration`, taken from the methods of each
                  paper): outflow reconstruction, inflow modulation, portosystemic shunt, combined, none. The only
                  contrast inside one cohort is splenectomy against none in Wang 2014.
    Outcome is carried along and is not evidence for any of this."""
    from scipy.stats import binomtest
    import audit
    A = {(g['study'], g['group']): g for g in audit.paired(d)['groups']}; b = d.dropna(subset=['zP', 'zF', 'zR']); detail = []
    for _, r in b.iterrows():
        g = A[(r.study, r.group)]
        detail.append(dict(study=r.study, group=r.group, cohort_id=r.cohort_id, configuration=r.configuration, outflow=r.outflow, n=int(r.n),
                           zP=round(float(r.zP), 3), zF=round(float(r.zF), 3), ratio=round(float(r.zR), 3), ratio_ci=[round(x, 3) for x in g['zR_ci_sampling']],
                           ratio_ci_by_correlation={k: [round(x, 3) for x in v] for k, v in g['zR_ci_by_correlation'].items()},
                           dispersion_reported=bool(g['flow_cv_reported'] and g['gradient_sd_reported']),
                           interval_below_one=bool(max(g['zR_ci_sampling'][1], *(v[1] for v in g['zR_ci_by_correlation'].values())) < 1),
                           outcome_pct=float(r.pct) if pd.notna(r.pct) else None))
    cohorts = {}
    for x in detail: cohorts.setdefault(x['cohort_id'], []).append(x)
    below = [all(x['ratio'] < 1 for x in v) for v in cohorts.values()]; K, k = len(cohorts), int(sum(below))
    conf = {}
    for x in detail: conf.setdefault(x['configuration'], []).append(x['ratio'])
    rec, none = conf.get('outflow reconstruction', []), conf.get('none', [])
    within = []
    for c, v in cohorts.items():
        if len({x['configuration'] for x in v}) > 1:
            v = sorted(v, key=lambda x: x['configuration'] == 'none'); u, w = v[0], v[-1]; g1, g2 = A[(u['study'], u['group'])], A[(w['study'], w['group'])]
            se = float(np.hypot(g1['zR_log_se'], g2['zR_log_se'])); q = u['ratio'] / w['ratio']
            within.append(dict(cohort=c, contrast=f"{u['configuration']} against {w['configuration']}", zR=[u['ratio'], w['ratio']], zF=[u['zF'], w['zF']], zP=[u['zP'], w['zP']],
                               ratio_of_zR=round(q, 3), ratio_ci=[round(q * float(np.exp(-1.96 * se)), 3), round(q * float(np.exp(1.96 * se)), 3)],
                               note='groups of one cohort that were not randomised: the manoeuvre was chosen for the recipients who needed it'))
    one = lambda a, n: float(binomtest(a, n, 0.5, alternative='greater').pvalue) if n else None
    centres = {}
    for x, c in zip(detail, b.centre_id): centres.setdefault(c, []).append(x['ratio'] < 1)
    arms = [x for x in detail if len(cohorts[x['cohort_id']]) > 1]
    sens = dict(unit_cohort=dict(below_one=k, of=K, one_sided_p=one(k, K)),
                unit_centre=dict(below_one=int(sum(all(v) for v in centres.values())), of=len(centres), one_sided_p=one(int(sum(all(v) for v in centres.values())), len(centres))),
                leave_one_cohort_out=[dict(left_out=c, below_one=int(sum(all(x['ratio'] < 1 for x in v) for cc, v in cohorts.items() if cc != c)), of=K - 1,
                                           one_sided_p=one(int(sum(all(x['ratio'] < 1 for x in v) for cc, v in cohorts.items() if cc != c)), K - 1)) for c in cohorts],
                one_group_per_cohort=[dict(group_dropped=f"{x['study']} / {x['group']}",
                                           below_one=int(sum(all(y['ratio'] < 1 for y in v if y is not x) for v in cohorts.values())), of=K,
                                           one_sided_p=one(int(sum(all(y['ratio'] < 1 for y in v if y is not x) for v in cohorts.values())), K)) for x in arms],
                smallest_one_sided_p_attainable=float(0.5 ** K),
                note='the six groups are six descriptions and five independent units; with five units the exact sign test cannot give less than 0.031, and leaving any cohort out gives 0.0625')
    return dict(groups_with_both=len(detail), cohorts=K, cohorts_below_one=k,
                sign_test_cohorts_one_sided_p=float(binomtest(k, K, 0.5, alternative='greater').pvalue), sign_test_cohorts_two_sided_p=float(binomtest(k, K, 0.5).pvalue),
                sign_test_sensitivity=sens,
                groups_below_one=int(sum(x['ratio'] < 1 for x in detail)), groups_with_interval_below_one=int(sum(x['interval_below_one'] for x in detail)),
                ratio_all=[x['ratio'] for x in detail], by_configuration={c: sorted(v) for c, v in conf.items()},
                no_manoeuvre_inside_the_range_of_outflow_reconstruction=bool(rec and none and min(rec) <= min(none) and max(none) <= max(rec)),
                within_cohort_contrasts=within,
                outcome_pct={f"{x['study']} / {x['group']}": x['outcome_pct'] for x in detail},
                note='the cohort is the unit of the sign test; groups of one cohort are not independent. Configuration was classified from the '
                     'methods of each paper; outcome is shown separately and is not evidence for the mechanism',
                detail=detail)


def p2_model_against_observed(d):
    """The prediction that does not use the pressure of any recipient: zR of a partial graft from the geometry of the
    network (the trunks a lobar graft leaves behind) and from the two inflows (the artery joins at the sinusoid, and the
    load is read from portal flow), over b in [0.5, 1] and an arterial flow per lobule between none and the donor's.
    Its only input from a series is the flow load. Compared here with what each group reports."""
    import audit
    from model import liver, human
    A = {(g['study'], g['group']): g for g in audit.paired(d)['groups']}; p = d.dropna(subset=['zF', 'zP']); rows = []
    lob = [liver.lobar_graft(liver.HUMAN_GENERATIONS, b, 1) for b in (1.0, 0.75, 0.5)]; cal = [human.liver(b) for b in (0.6, 2 / 3, 0.75, 0.85, 0.93)]
    measured = human.two_inflows_measured()
    for _, r in p.iterrows():
        z = [liver.two_inflows(r.zF, l['share_inlet'], l['zR'], arterial_load=phi)['zR'] for l in lob for phi in (0.0, 1.0)]
        z += [liver.two_inflows(r.zF, c['share_inlet'], 1 - c['trunks_removed'], arterial_load=phi)['zR'] for c in cal for phi in (0.0, 1.0)]
        lo, hi = float(min(z)), float(max(z)); g = A[(r.study, r.group)]; obs = float(r.zR); ci = g['zR_ci_sampling']
        row = dict(study=r.study, group=r.group, configuration=r.configuration, zF=round(float(r.zF), 3), zR_predicted=[round(lo, 3), round(hi, 3)], zR_observed=round(obs, 3),
                   zR_observed_ci=[round(x, 3) for x in ci], observed_minus_predicted=[round(obs - hi, 3), round(obs - lo, 3)],
                   share_of_the_departure_from_one_predicted=[round((1 - hi) / (1 - obs), 2), round((1 - lo) / (1 - obs), 2)],
                   interval_reaches_the_prediction=bool(ci[1] >= lo))
        if r.study in measured and 'zR_model' in measured[r.study]: row['zR_predicted_with_the_arterial_flow_measured_in_that_series'] = [round(x, 3) for x in measured[r.study]['zR_model']]
        rows.append(row)
    lo_all = min(x['zR_predicted'][0] for x in rows); obs = np.array([x['zR_observed'] for x in rows])
    return dict(groups=rows, lowest_predicted=round(lo_all, 3), mean_observed=round(float(obs.mean()), 3),
                mean_predicted_range=[round(float(np.mean([x['zR_predicted'][0] for x in rows])), 3), round(float(np.mean([x['zR_predicted'][1] for x in rows])), 3)],
                groups_below_the_prediction=int(sum(x['zR_observed'] < x['zR_predicted'][0] for x in rows)),
                groups_whose_interval_is_below_the_prediction=int(sum(not x['interval_reaches_the_prediction'] for x in rows)),
                share_of_the_departure_predicted=[min(x['share_of_the_departure_from_one_predicted'][0] for x in rows), max(x['share_of_the_departure_from_one_predicted'][1] for x in rows)],
                verdict='the network gives the direction of the departure of zR from one and under-predicts its size in every group',
                note='inputs: geometry of the network, arterial share of the donor, flow load of the group. No pressure of a recipient enters')



CENTRAL = dict(b=0.75, arterial_load=1.0)     # the network of the human cast (human.B_HUMAN) with the donor's arterial flow per lobule
ARTERIAL_GRID = tuple(np.round(np.linspace(0.0, 1.0, 11), 2))


def prediction_metrics(observed, central, admissible, interval, observed_ci):
    """Four separate statements about a prediction, none of which uses the midpoint of a range of configurations.
    observed       zR of each group
    central        zR of one specified configuration, per group
    admissible     zR of every admissible configuration, shape (configurations, groups): one configuration is used for
                   all groups at once, as a model would be
    interval       (lowest, highest) prediction per group
    observed_ci    95% interval of each observation (sampling and references)
    Returns the error against the specified configuration, the smallest error any single admissible configuration
    attains (and which), and the coverage of the predicted interval by the observations."""
    obs = np.asarray(observed, float); c = np.asarray(central, float); adm = np.asarray(admissible, float)
    e = obs - c; rm = np.sqrt(np.mean((obs[None, :] - adm) ** 2, axis=1)); i = int(np.argmin(rm))
    inside = [bool(lo <= o <= hi) for o, (lo, hi) in zip(obs, interval)]
    overlap = [bool(ci[0] <= hi and ci[1] >= lo) for ci, (lo, hi) in zip(observed_ci, interval)]
    return dict(central=dict(mean_error=round(float(e.mean()), 3), rmse=round(float(np.sqrt(np.mean(e ** 2))), 3), errors=[round(float(x), 3) for x in e]),
                best_admissible=dict(rmse=round(float(rm[i]), 3), mean_error=round(float((obs - adm[i]).mean()), 3), index=i),
                coverage=dict(observed_inside_the_predicted_interval=int(sum(inside)), observed_interval_overlaps_it=int(sum(overlap)), of=len(obs)))


def p2_performance(d, model, uncertainty):
    """How well does the rigid network predict zR? Four statements, kept apart:
    central      the error against one specified configuration (CENTRAL: the network of the human cast at b = 0.75, each
                 lobule keeping the donor's arterial flow), chosen before this comparison and not adjusted to it
    best         the smallest root-mean-square error that any single admissible configuration (eight networks x eleven
                 arterial loads from none to the donor's) attains on all six groups at once: the best the rigid
                 network can do
    coverage     how many observations fall inside the predicted interval, and how many observed 95% intervals
                 (sampling and references, p2_uncertainty) reach it
    references   the error of the central configuration when the donor flow, the normal gradient and the assumed CVP
                 are varied together over the 60 combinations of p2_reference_sweep: the observed zR and the flow load
                 at which the network is evaluated both change"""
    from model import liver, human, compartments as cm
    p = d.dropna(subset=['zF', 'zP']); G8 = [g for g in cm.geometries()]; zF = p.zF.to_numpy(float); obs = p.zR.to_numpy(float)
    cen = [g for g in G8 if g[0] == f"human b={CENTRAL['b']:.2f}"][0]
    rig = lambda zf, geom, phi: cm.rigid(zf, geom, phi)
    central = np.array([rig(z, cen, CENTRAL['arterial_load']) for z in zF])
    configs = [(g, phi) for g in G8 for phi in ARTERIAL_GRID]; adm = np.array([[rig(z, g, phi) for z in zF] for g, phi in configs])
    U = {(x['study'], x['group']): x for x in uncertainty['groups']}; M = {(x['study'], x['group']): x for x in model['groups']}
    keys = list(zip(p.study, p.group)); interval = [M[k]['zR_predicted'] for k in keys]; oci = [U[k]['ci_total'] for k in keys]
    met = prediction_metrics(obs, central, adm, interval, oci); g, phi = configs[met['best_admissible'].pop('index')]
    met['best_admissible'].update(network=g[0], arterial_load=float(phi))
    # the references varied together
    sweep = []
    for f in REF_FLOWS:
        for ng in REF_GRADIENTS:
            for c in REF_CVPS:
                ref = np.where(p.donor_PVF_per_100g_ref.isna() | (p.donor_PVF_per_100g_ref == indices.DEFAULT_DONOR_REF), f, p.donor_PVF_per_100g_ref)
                grad = np.where(p.gradient_mmHg.notna(), p.gradient_mmHg, p.PVP - np.where(p.CVP_basis.eq('imputed'), c, p.CVP))
                zf = (p.PVF_per_100g / ref).to_numpy(float); zr = (grad / ng) / zf
                e = zr - np.array([rig(z, cen, CENTRAL['arterial_load']) for z in zf]); sweep.append((float(np.sqrt(np.mean(e ** 2))), float(e.mean())))
    met['references'] = dict(rmse=[round(min(x[0] for x in sweep), 3), round(max(x[0] for x in sweep), 3)], mean_error=[round(min(x[1] for x in sweep), 3), round(max(x[1] for x in sweep), 3)], combinations=len(sweep))
    met.update(central_configuration=dict(network=cen[0], arterial_load=CENTRAL['arterial_load'], zR=[round(float(x), 3) for x in central]),
               groups=[f'{a} / {b}' for a, b in keys],
               note='no statement uses the midpoint of the predicted range, which is not a configuration of the network')
    return met


REF_FLOWS = (80.0, 85.0, 90.0, 95.0, 100.0)      # donor PVF per 100 g, default 90
REF_GRADIENTS = (4.0, 5.0, 6.0)                   # normal portocaval gradient, default 5 mmHg
REF_CVPS = (3.0, 5.0, 7.0, 9.0)                   # assumed CVP where a series reports none, default 5 mmHg


def p2_reference_sweep(d, model_lowest=0.81):
    """Could the normalisation alone put every paired group below one? Both denominators of zR were fixed from donor
    physiology before the extraction, so they are swept here: the default donor reference flow (used where a series
    reports no donor value of its own; a stored reference equal to the default 90 is the default, only Troisi and Chan
    report their own) over 80-100 mL/min/100 g, the normal gradient over 4-6 mmHg, and, separately, the CVP assumed
    where none is reported over 3-9 mmHg. For each combination the largest zR among the groups that
    report both loads is recorded; the combination that overturns the result, if any, is named."""
    raw = d[d.in_main_analysis & d.PVF_per_100g.notna() & (d.PVP.notna() | d.gradient_mmHg.notna())].copy()

    def zr(ref_flow, normal_gradient, cvp):
        stored = raw.donor_PVF_per_100g_ref
        ref = np.where(stored.isna() | (stored == indices.DEFAULT_DONOR_REF), ref_flow, stored)
        grad = np.where(raw.gradient_mmHg.notna(), raw.gradient_mmHg,
                        raw.PVP - np.where(raw.CVP_basis.eq('imputed'), cvp, raw.CVP))
        return np.asarray((grad / normal_gradient) / (raw.PVF_per_100g / ref), float)

    def zr_max(ref_flow, normal_gradient, cvp):
        z = zr(ref_flow, normal_gradient, cvp); i = int(np.argmax(z))
        return float(z[i]), f"{raw.study.iloc[i]} / {raw.group.iloc[i]}"

    # the three references varied together, and where each conclusion stops holding
    names = [f"{a} / {b}" for a, b in zip(raw.study, raw.group)]; joint = []; per_group = {n: dict(below_one=0, below_model=0) for n in names}
    for f in REF_FLOWS:
        for g in REF_GRADIENTS:
            for c in REF_CVPS:
                z = zr(f, g, c); joint.append(dict(donor_ref_flow=f, normal_gradient=g, assumed_cvp=c, max_zR=round(float(z.max()), 3),
                                                   groups_below_one=int((z < 1).sum()), groups_below_model=int((z < model_lowest).sum())))
                for n, v in zip(names, z): per_group[n]['below_one'] += int(v < 1); per_group[n]['below_model'] += int(v < model_lowest)

    grid = []
    for f in REF_FLOWS:
        for g in REF_GRADIENTS:
            m, who = zr_max(f, g, 5.0)
            grid.append(dict(donor_ref_flow=f, normal_gradient=g, assumed_cvp=5.0, max_zR=round(m, 3), max_group=who,
                             all_below_one=bool(m < 1)))
    cvp = []
    for c in REF_CVPS:
        m, who = zr_max(90.0, 5.0, c)
        cvp.append(dict(donor_ref_flow=90.0, normal_gradient=5.0, assumed_cvp=c, max_zR=round(m, 3), max_group=who,
                        all_below_one=bool(m < 1)))
    corner_m, corner_who = zr_max(100.0, 4.0, 3.0)
    return dict(groups=int(len(raw)), grid=grid, combinations=len(grid),
                combinations_all_below_one=int(sum(x['all_below_one'] for x in grid)),
                exceptions=[x for x in grid if not x['all_below_one']],
                cvp_sweep=cvp, cvp_sweep_all_below_one=bool(all(x['all_below_one'] for x in cvp)),
                worst_corner=dict(donor_ref_flow=100.0, normal_gradient=4.0, assumed_cvp=3.0, max_zR=round(corner_m, 3),
                                  max_group=corner_who),
                joint=dict(combinations=len(joint), all_below_one=int(sum(x['groups_below_one'] == len(names) for x in joint)),
                           all_below_model=int(sum(x['groups_below_model'] == len(names) for x in joint)), model_lowest=float(model_lowest),
                           per_group={n: dict(share_below_one=round(v['below_one'] / len(joint), 3), share_below_model=round(v['below_model'] / len(joint), 3)) for n, v in per_group.items()},
                           not_all_below_one=[x for x in joint if x['groups_below_one'] < len(names)],
                           by_normal_gradient={f'{g:g} mmHg': dict(all_below_one=int(sum(x['groups_below_one'] == len(names) for x in joint if x['normal_gradient'] == g)),
                                                                    all_below_model=int(sum(x['groups_below_model'] == len(names) for x in joint if x['normal_gradient'] == g)),
                                                                    of=int(sum(1 for x in joint if x['normal_gradient'] == g))) for g in REF_GRADIENTS},
                           note='60 combinations of donor reference flow, normal gradient and assumed CVP, varied together'),
                note='flow and gradient references crossed on a grid; the assumed CVP swept separately at the default '
                     'references; the worst corner crosses all three and is reported for completeness')


FLOW_CVS = (0.25, 0.35, 0.45, 0.55, 0.70)          # SD / mean of portal flow per 100 g where a series reports none (reported: 0.32 to 0.56)
GRADIENT_SDS = (2.5, 3.5, 4.5, 6.0)               # SD of the gradient in mmHg where a series reports none (reported: 2.6 to 4.5)
CORRELATIONS = (-0.5, 0.0, 0.5)                   # between the sampling errors of the two group means; no series reports it


def p2_uncertainty(d, model, draws=40000, seed=0):
    """How uncertain is each zR, and from what? Three sources, kept apart.

    sampling     the error of the two group means. With the dispersions each paper reports and, where it reports none,
                 the largest reported one, at zero correlation: `ci_sampling`. Where a dispersion is not reported it is
                 then swept over FLOW_CVS and GRADIENT_SDS, and the correlation between the two means over CORRELATIONS
                 for every group: `ci_sampling_widest` is the widest interval of the sweep, and the sweep says in how
                 many of its combinations the upper limit stays below one and below the prediction of the network.
    references   the denominators of the two loads, which are not data of the group: the donor portal flow per 100 g
                 (80 to 100 where the default of 90 is used; the same relative range around the value a series measured
                 in its own donors), the normal gradient (4 to 6 mmHg) and the CVP assumed where none is reported
                 (3 to 9 mmHg), each drawn uniformly over the range of p2_reference_sweep: `ci_references`.
    both         `ci_total` (reported or largest dispersions, zero correlation) and `ci_total_widest` (the widest of the
                 sweep). The references are common to the groups that share them, so these intervals are marginal, one
                 group at a time; what the references do to all groups at once is p2_reference_sweep.
    The uniform ranges are a choice, not a measurement: the interval answers "if the references could be anywhere in
    their plausible range", not a 95% statement about the donors."""
    import audit
    A = {(g['study'], g['group']): g for g in audit.paired(d, correlations=())['groups']}; pred = {(g['study'], g['group']): g['zR_predicted'][0] for g in model['groups']}
    rng = np.random.default_rng(seed); u = rng.uniform(size=(3, draws)); e = rng.standard_normal(draws); rows = []
    span = lambda x: [round(float(v), 3) for v in np.percentile(x, [2.5, 97.5])]
    for _, r in d.dropna(subset=['zF', 'zP']).iterrows():
        g = A[(r.study, r.group)]; rep = audit.REPORTED[(r.study, r.group)]; grad = 5 * r.zP; lo = pred[(r.study, r.group)]
        se = lambda cv, sd, rho: audit.ratio_log_se(cv, rep['n_flow'], sd, rep['n_grad'], grad, rho)
        cvs = (g['flow_cv_used'],) if g['flow_cv_reported'] else FLOW_CVS; sds = (g['gradient_sd_used'],) if g['gradient_sd_reported'] else GRADIENT_SDS
        sweep = [(cv, sd, rho, se(cv, sd, rho)) for cv in cvs for sd in sds for rho in CORRELATIONS]; worst = max(sweep, key=lambda x: x[3])
        upper = np.array([r.zR * np.exp(1.96 * x[3]) for x in sweep]); se0 = g['zR_log_se']; ci = lambda x: [round(float(r.zR * np.exp(-1.96 * x)), 3), round(float(r.zR * np.exp(1.96 * x)), 3)]
        # the references, drawn over their ranges; a measured gradient has no CVP in it
        f_rel = (REF_FLOWS[0] + u[0] * (REF_FLOWS[-1] - REF_FLOWS[0])) / indices.DEFAULT_DONOR_REF
        ng = REF_GRADIENTS[0] + u[1] * (REF_GRADIENTS[-1] - REF_GRADIENTS[0]); cvp = REF_CVPS[0] + u[2] * (REF_CVPS[-1] - REF_CVPS[0])
        gradient = np.full(draws, float(r.gradient_mmHg)) if pd.notna(r.gradient_mmHg) else r.PVP - (cvp if r.CVP_basis == 'imputed' else r.CVP)
        z_ref = (gradient / ng) / (r.zF / f_rel); tot = z_ref * np.exp(se0 * e); tot_w = z_ref * np.exp(worst[3] * e)
        v_s, v_r = se0 ** 2, float(np.var(np.log(z_ref)))
        ups = dict(ci_sampling=r.zR * np.exp(1.96 * se0), ci_sampling_widest=r.zR * np.exp(1.96 * worst[3]), ci_references=np.percentile(z_ref, 97.5),
                   ci_total=np.percentile(tot, 97.5), ci_total_widest=np.percentile(tot_w, 97.5))
        rows.append(dict(study=r.study, group=r.group, cohort_id=r.cohort_id, zR=round(float(r.zR), 3), zR_predicted_lowest=lo, upper_limits={k: float(v) for k, v in ups.items()},
                         flow_cv_reported=bool(g['flow_cv_reported']), gradient_sd_reported=bool(g['gradient_sd_reported']),
                         cvp_assumed=bool(r.CVP_basis == 'imputed'), donor_flow_measured_in_series=bool(g['donor_flow_reference_measured_in_series']),
                         ci_sampling=ci(se0), ci_sampling_widest=ci(worst[3]), widest_at=dict(flow_cv=float(worst[0]), gradient_sd=float(worst[1]), correlation=float(worst[2])),
                         sweep=dict(combinations=len(sweep), upper_limit_below_one=int((upper < 1).sum()), upper_limit_below_prediction=int((upper < lo).sum()),
                                    upper_limit_range=[round(float(upper.min()), 3), round(float(upper.max()), 3)]),
                         ci_references=span(z_ref), ci_total=span(tot), ci_total_widest=span(tot_w),
                         prob_below_one=round(float((tot < 1).mean()), 3), prob_below_prediction=round(float((tot < lo).mean()), 3),
                         share_of_variance_from_references=round(v_r / (v_r + v_s), 2)))
    under = lambda x, key, bound: x['upper_limits'][key] < (bound if bound else x['zR_predicted_lowest'])        # on the unrounded limit
    count = lambda key, bound: int(sum(under(x, key, bound) for x in rows))
    kinds = ('ci_sampling', 'ci_sampling_widest', 'ci_references', 'ci_total', 'ci_total_widest')
    cohorts = sorted({x['cohort_id'] for x in rows}); coh = lambda key, bound: int(sum(all(under(x, key, bound) for x in rows if x['cohort_id'] == c) for c in cohorts))
    return dict(groups=rows, of=len(rows), cohorts=len(cohorts),
                groups_with_interval_below_one={k: count(k, 1.0) for k in kinds}, groups_with_interval_below_prediction={k: count(k, None) for k in kinds},
                cohorts_with_interval_below_one={k: coh(k, 1.0) for k in kinds}, cohorts_with_interval_below_prediction={k: coh(k, None) for k in kinds},
                grids=dict(flow_cv=list(FLOW_CVS), gradient_sd_mmHg=list(GRADIENT_SDS), correlation=list(CORRELATIONS), donor_flow=[REF_FLOWS[0], REF_FLOWS[-1]],
                           normal_gradient_mmHg=[REF_GRADIENTS[0], REF_GRADIENTS[-1]], assumed_cvp_mmHg=[REF_CVPS[0], REF_CVPS[-1]], draws=draws),
                note='sampling: error of the two group means. references: donor flow, normal gradient and assumed CVP drawn uniformly over their ranges. '
                     'The intervals with references are marginal for each group; the references are shared between groups')


def same_cohort_contrast(d):
    """The only place where flow and pressure can be set against each other on the same recipients: a cohort whose
    groups report both loads and an outcome. Fisher's exact test of the outcome between its two groups, with the
    direction in which each load differs. One cohort, not randomised: a description, not a test of the two indices."""
    from scipy.stats import fisher_exact
    out = []
    for c, g in d[d.zF.notna() & d.zP.notna() & d.events.notna()].groupby('cohort_id'):
        if len(g) != 2: continue
        g = g.sort_values('pct'); lo, hi = g.iloc[0], g.iloc[1]; tab = [[int(hi.events), int(hi.n - hi.events)], [int(lo.events), int(lo.n - lo.events)]]; odds, pv = fisher_exact(tab)
        out.append(dict(cohort=c, study=str(hi.study), groups=[str(hi.group), str(lo.group)], outcome_pct=[float(hi.pct), float(lo.pct)], events=[int(hi.events), int(lo.events)], n=[int(hi.n), int(lo.n)],
                        odds_ratio=float(odds), fisher_two_sided_p=float(pv), zP=[float(hi.zP), float(lo.zP)], zF=[float(hi.zF), float(lo.zF)],
                        worse_group_has_higher_pressure_load=bool(hi.zP > lo.zP), worse_group_has_higher_flow_load=bool(hi.zF > lo.zF)))
    return out


def p3_within_cohort(d):
    """Primary analysis: hierarchical binomial model with the pressure index centred within each study, fitted by
    NUTS. Reported in the prespecified order (primary outcome, secondary outcome, all outcomes combined), with the
    sensitivity to the prior on tau, the overlap analysis for the Kyoto series, and the conditional
    internal-external validation, which scores the slope without fitting anything on the held-out centre."""
    import bayes, cache
    c = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]
    FITS = ('bayes_main', 'bayes_iecv', 'spec_curve', 'bayes_priors', 'bayes_measurement', 'bayes_recovery', 'bayes_extra')
    old = cache.stale(FITS)
    if old:                                          # another table, other model code, or LIPNET_REBUILD: fit again, all of them
        print('cached fits not computed from the current table and model code: ' + ', '.join(old) + '; running bayes.py', flush=True); bayes.main()
        if cache.stale(FITS) and not os.environ.get('LIPNET_REBUILD'): raise RuntimeError('bayes.py did not write ' + ', '.join(cache.stale(FITS)))
    read = lambda name, make=None: json.load(open(f'{OUT}/{name}.json'))
    main = read('bayes_main'); iecv_ = read('bayes_iecv'); spec = read('spec_curve'); prio = read('bayes_priors'); meas = read('bayes_measurement'); recov = read('bayes_recovery')
    pairs = []
    for st, g in c.groupby(indices.strata(c).values):
        g = g.sort_values('zP')
        if len(g) >= 2: pairs.append(dict(stratum=st, higher_zP_worse=bool(g.pct.iloc[-1] >= g.pct.iloc[0])))
    extra = read('bayes_extra')
    acc = dict(series=int(d.study.nunique()), centres=int(d.centre_id.nunique()), groups=int(len(d)),
               groups_with_zP=int(d.zP.notna().sum()), groups_with_zF=int(d.zF.notna().sum()),
               series_with_a_within_study_contrast=len(set(x.split(' / ')[0] for x in main['exploratory (all outcomes)']['strata'])),
               series_with_one_group_only=int(sum(1 for _, g in d[d.zP.notna()].groupby('study') if len(g) == 1)))
    mus = [x['mu_beta'] for x in spec if x.get('mu_beta') is not None]
    return dict(primary=main['primary (SFSS or early dysfunction)'], secondary=main['secondary (mortality or graft loss)'],
                exploratory=main['exploratory (all outcomes)'],
                every_pair_in_predicted_direction=bool(all(p['higher_zP_worse'] for p in pairs)), pairs=pairs,
                prior_sensitivity=[dict(tau_prior=r['tau_prior'], mu_beta=r['mu_beta'], ci=r['mu_beta_ci'], tau_beta=r['tau_beta'], prob_positive=r['prob_mu_positive'],
                                        slope_in_a_new_centre=r['prediction_interval'], prob_new_centre_positive=r['prob_new_study_positive']) for r in prio],
                interval_excludes_zero_under_every_prior=bool(all(r['mu_beta_ci'][0] > 0 for r in prio)),
                slope_in_a_new_centre={k: dict(interval=v['prediction_interval'], prob_positive=v['prob_new_study_positive']) for k, v in main.items()},
                same_cohort_contrast=same_cohort_contrast(d),
                reading='an association between groups of one cohort (ecological), with groups often defined by the pressure itself or by a manoeuvre chosen '
                        'because of it: not the effect of pressure on a patient and not a causal effect',
                conditional_iecv=[dict(centre=v['held_out'], mu_from_others=v['mu_from_others'],
                                       log_score_gain=v['total_gain'], directions_correct=v['all_directions_correct']) for v in iecv_],
                accounting=acc, index_contrast=extra['index_contrast'], with_cutoff_groups=extra['with_cutoff_groups'],
                specification_curve=spec, mu_range_across_specifications=[float(min(mus)), float(max(mus))],
                all_specifications_positive=bool(all(x['mu_beta'] > 0 for x in spec if x.get('mu_beta') is not None)),
                measurement_monte_carlo=meas, design_recovery=recov)

def p4_thresholds():
    """Exploratory comparison of haemodynamic thresholds across liver diseases and surgical settings. Not a statistical
    result, not a physiological threshold and not an independent prediction: the thresholds are each a gradient of about
    10 mmHg before any normalisation, so dividing by the 5 mmHg reference gives 2 by construction, and the resection
    limit is the threshold of cirrhosis applied to surgery. HVPG is wedged minus free hepatic venous pressure and
    reads the sinusoidal gradient; it equals the portocaval gradient only when the resistance is sinusoidal and the
    hepatic veins and cava add nothing. PVP minus CVP in a graft includes the outflow. The two are placed on one axis
    without any validated equivalence. The landing at 2 depends on the 5 mmHg reference, so 4 and 6 mmHg are shown."""
    z = indices.zP
    sens = {f'{ref:g} mmHg reference': {k: round(v / ref, 2) for k, v in
            dict(cirrhosis_HVPG_10=10, graft_PVP15_CVP5=10, resection_HVPG_10=10, variceal_HVPG_12=12).items()}
            for ref in (4, 5, 6)}
    return dict(title='Exploratory comparison of haemodynamic thresholds across liver diseases and surgical settings',
                reading='comparison of scales, exploratory; the thresholds are each about 10 mmHg, so they sit at 2 by the choice of a 5 mmHg reference; '
                        'HVPG and the portocaval gradient of a graft are different measurements',
                agreement_is_by_construction=True, thresholds_mmHg=dict(cirrhosis_HVPG=10, graft_PVP_minus_CVP=10, resection_HVPG=10, variceal_HVPG=12),
                reference_sensitivity=sens, clinically_significant_portal_hypertension=dict(source='HVPG >= 10 mmHg (Baveno)', zP=z(gradient=10)),
                graft_threshold=dict(source='PVP 15 mmHg at CVP 5 (ILTS 2023)', zP=z(15, 5)),
                variceal_bleeding=dict(source='HVPG >= 12 mmHg', zP=z(gradient=12)),
                resection_limit=dict(source='HVPG 10 mmHg, any resection puts the remnant above', zP=z(gradient=10)),
                all_at_two=bool(abs(z(gradient=10) - 2) < 1e-9 and abs(z(15, 5) - 2) < 1e-9))


def main():
    os.makedirs(OUT, exist_ok=True)
    d = indices.compute(pd.read_csv(f'{ROOT}/data/series.tsv', sep='\t'))
    res = dict(P1_allometry=p1_allometry(), P2_discordance=p2_discordance(d),
               P3_within_cohort=p3_within_cohort(d), P4_thresholds=p4_thresholds())
    res['P2_discordance']['model_against_observed'] = p2_model_against_observed(d)
    res['P2_discordance']['uncertainty'] = p2_uncertainty(d, res['P2_discordance']['model_against_observed'])
    res['P2_discordance']['reference_sweep'] = p2_reference_sweep(d, model_lowest=res['P2_discordance']['model_against_observed']['lowest_predicted'])
    res['P2_discordance']['performance'] = p2_performance(d, res['P2_discordance']['model_against_observed'], res['P2_discordance']['uncertainty'])
    json.dump(res, open(f'{OUT}/predictions.json', 'w'), indent=1)
    a, b2, c3, d4 = res['P1_allometry'], res['P2_discordance'], res['P3_within_cohort'], res['P4_thresholds']
    print('P1  comparative anatomy')
    print(f"    over b = {a['b_range'][0]:.2f}-{a['b_range'][1]:.2f}: exponent of the design scale tau0 stays in "
          f"[{a['tau0_exponent_range'][0]:+.3f}, {a['tau0_exponent_range'][1]:+.3f}] -> invariant: {a['tau0_invariant']}")
    print(f"    vascular mass exponent spans {a['vascular_mass_exponent_range'][0]:.2f}-{a['vascular_mass_exponent_range'][1]:.2f}; "
          f"{100*a['fraction_compatible']:.0f}% of the range is compatible with the observed {a['observed_range']}")
    print('P2  the resistance load of a partial graft: prediction of the network against what the cohorts report')
    m = b2['model_against_observed']
    print(f"    {b2['groups_with_both']} groups from {b2['cohorts']} cohorts report both loads: zR {sorted(b2['ratio_all'])}")
    print(f"    cohorts below one: {b2['cohorts_below_one']} of {b2['cohorts']} (exact sign test, one-sided p = {b2['sign_test_cohorts_one_sided_p']:.3f}, two-sided {b2['sign_test_cohorts_two_sided_p']:.3f}); "
          f"groups whose interval stays below one at any correlation of the two means: {b2['groups_with_interval_below_one']} of {b2['groups_with_both']}")
    print(f"    network, from geometry and the two inflows: lowest value {m['lowest_predicted']}; mean predicted range {m['mean_predicted_range']} against mean observed {m['mean_observed']}"
          f"; {m['groups_below_the_prediction']} of {len(m['groups'])} groups below it, "
          f"{m['groups_whose_interval_is_below_the_prediction']} with the whole interval below; share of the departure predicted {m['share_of_the_departure_predicted']}")
    pf = b2['performance']
    print(f"    performance: central configuration ({pf['central_configuration']['network']}, arterial load {pf['central_configuration']['arterial_load']:g}) mean error {pf['central']['mean_error']:+.3f}, rmse {pf['central']['rmse']:.3f}; "
          f"best single admissible configuration rmse {pf['best_admissible']['rmse']:.3f} ({pf['best_admissible']['network']}, {pf['best_admissible']['arterial_load']:g}); "
          f"observations inside the predicted interval {pf['coverage']['observed_inside_the_predicted_interval']} of {pf['coverage']['of']}, observed intervals reaching it {pf['coverage']['observed_interval_overlaps_it']}; "
          f"rmse over the references {pf['references']['rmse']}")
    print(f"    by configuration: {b2['by_configuration']}; the group with no manoeuvre is inside the range of outflow reconstruction: {b2['no_manoeuvre_inside_the_range_of_outflow_reconstruction']}")
    for w in b2['within_cohort_contrasts']: print(f"    inside {w['cohort']}: {w['contrast']}, zR {w['zR']}, ratio {w['ratio_of_zR']} ({w['ratio_ci'][0]} to {w['ratio_ci'][1]})")
    sw = b2['reference_sweep']
    print(f"    reference sweep: every paired group below 1 in {sw['combinations_all_below_one']} of {sw['combinations']} "
          f"flow x gradient combinations; exceptions {[(x['donor_ref_flow'], x['normal_gradient'], x['max_zR']) for x in sw['exceptions']]}; "
          f"CVP 3-9 all below 1: {sw['cvp_sweep_all_below_one']}; worst corner (100, 4, CVP 3): {sw['worst_corner']['max_zR']}")
    j = sw['joint']
    print(f"    the three references together ({j['combinations']} combinations): every group below one in {j['all_below_one']}, every group below the network in {j['all_below_model']}; by normal gradient {j['by_normal_gradient']}")
    print('P3  within-cohort ordering (primary analysis: hierarchical, zP centred within study)')
    a = c3['accounting']
    print(f"    material: {a['series']} series from {a['centres']} centres, {a['groups']} groups; {a['groups_with_zP']} with a pressure index "
          f"and {a['groups_with_zF']} with a flow index. A slope needs two groups of the same outcome at different gradients, which "
          f"{a['series_with_a_within_study_contrast']} series provide for pressure and 4 for flow; "
          f"{a['series_with_one_group_only']} report a single group and inform the level, not the slope.")
    ic = c3['index_contrast']
    print('    exploratory, the same model on each index, both restricted to the primary outcome (the pressure arm is the primary fit):')
    for k in ('zP, pressure per lobule', 'zF, flow per lobule'):
        r = ic[k]
        print(f"      {k:24s} {r['groups']:2d} groups / {len(r['strata'])} series: mu {r['mu_beta']:+.2f} "
              f"({r['mu_beta_ci'][0]:+.2f} to {r['mu_beta_ci'][1]:+.2f}), P(mu>0) = {r['prob_positive']:.3f}")
    print(f"      the two arms are different cohorts (shared: {ic['shared_series']}); the flow arm rests on {ic['zF, flow per lobule']['events']} events "
          f"in groups of median {ic['zF, flow per lobule']['median_group_size']:.0f}. The two slopes are not compared: a joint model needs cohorts with both loads and an outcome, and there is one")
    for w in c3['same_cohort_contrast']:
        print(f"    inside one cohort ({w['study']}): outcome {w['outcome_pct'][0]:.1f}% against {w['outcome_pct'][1]:.1f}% (Fisher p = {w['fisher_two_sided_p']:.3f}); the worse group has the "
              f"{'higher' if w['worse_group_has_higher_pressure_load'] else 'lower'} pressure load ({w['zP'][0]:.2f} against {w['zP'][1]:.2f}) and the "
              f"{'higher' if w['worse_group_has_higher_flow_load'] else 'lower'} flow load ({w['zF'][0]:.2f} against {w['zF'][1]:.2f})")
    for k, v in c3['slope_in_a_new_centre'].items(): print(f"    slope in a new centre, {k}: {v['interval'][0]:+.2f} to {v['interval'][1]:+.2f}, P>0 {v['prob_positive']:.3f}")
    print(f"    interval of the mean slope excludes zero under every prior on tau: {c3['interval_excludes_zero_under_every_prior']}")
    cg = c3['with_cutoff_groups']
    print(f"    adding the series whose groups are defined by a gradient cut-off: {cg['groups']} groups / {len(cg['strata'])} strata, "
          f"mu {cg['mu_beta']:+.2f} ({cg['mu_beta_ci'][0]:+.2f} to {cg['mu_beta_ci'][1]:+.2f})")
    for k in ('primary', 'secondary', 'exploratory'):
        r = c3[k]
        print(f"    {k:12s} {r['groups']:2d} groups / {len(r['strata'])} strata: mu {r['mu_beta']:+.2f} "
              f"({r['mu_beta_ci'][0]:+.2f} to {r['mu_beta_ci'][1]:+.2f}), OR {r['OR_per_zP']:.2f}, "
              f"P(mu>0) = {r['prob_mu_positive']:.3f}, tau {r['tau_beta']:.2f}; R-hat {r['rhat']['max']:.3f}, "
              f"ESS {r['ess']['min']:.0f}, divergences {r['divergences']}, PPC {r['ppc_inside']}/{r['ppc_total']}")
    print(f"    every pair in the predicted direction: {c3['every_pair_in_predicted_direction']}")
    print('    prior on tau: ' + '; '.join(f"HN({p['tau_prior']}) mu {p['mu_beta']:+.2f} ({p['ci'][0]:+.2f},{p['ci'][1]:+.2f}) P>0 {p['prob_positive']:.3f}" for p in c3['prior_sensitivity']))
    print('    conditional validation (slope only, no intercept fitted on the held-out centre):')
    for v in c3['conditional_iecv']:
        print(f"      {v['centre']:10s} mu from the others {v['mu_from_others']:+.2f}, log-score gain over no association {v['log_score_gain']:+.2f}, directions correct {v['directions_correct']}")
    print(f"    specification curve: mu between {c3['mu_range_across_specifications'][0]:+.2f} and "
          f"{c3['mu_range_across_specifications'][1]:+.2f} across {len(c3['specification_curve'])} specifications; all positive: {c3['all_specifications_positive']}")
    if c3['measurement_monte_carlo']:
        mm = c3['measurement_monte_carlo']
        print(f"    with measurement uncertainty propagated: mu {mm['mu_beta']:+.2f} ({mm['mu_beta_ci'][0]:+.2f} to {mm['mu_beta_ci'][1]:+.2f}), P>0 {mm['prob_positive']:.2f}")
    if c3['design_recovery']:
        print('    design recovery: ' + '; '.join(f"true {r['mu_true']:.1f} -> {r['mean_estimate']:.2f} (coverage {r['coverage']:.2f}, {r['coverage_ci'][0]:.2f} to {r['coverage_ci'][1]:.2f}; power {r['power_to_exclude_zero']:.2f}, {r['power_ci'][0]:.2f} to {r['power_ci'][1]:.2f}; {r['reps']} refits)" for r in c3['design_recovery'])
              + f"; heterogeneity set to {c3['design_recovery'][0]['tau_true']}, one level per stratum")
    print('P4  exploratory comparison of haemodynamic thresholds across liver diseases and surgical settings')
    print('    each threshold is a gradient of about 10 mmHg, so 2 follows from the 5 mmHg reference; HVPG and the portocaval gradient of a graft are different measurements')
    for k, v in d4.items():
        if isinstance(v, dict) and 'source' in v: print(f"    {v['source']:45s} zP = {v['zP']:.2f}")
    print('    with a reference gradient of 4, 5 or 6 mmHg the same thresholds give:')
    for k, v in d4['reference_sensitivity'].items(): print(f"      {k:18s} {v}")


if __name__ == '__main__':
    main()
