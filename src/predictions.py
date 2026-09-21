"""Test the predictions of the perfusion-network model against independent data.

    python src/predictions.py        -> results/predictions.json, printed summary

The model (src/model/) describes the liver as a dissipation-optimal perfusion network at fixed vascular maintenance
cost. It fixes the wall shear set-point tau0 = sqrt(D*/C), the scaling of the optimum with the number of lobules,
and the local relation dP = f * R, from which the two indices follow. It says nothing about recipient severity,
donor age or steatosis. Four predictions follow, and each is tested here against data the model was not fitted to:

  P1  Comparative anatomy. With portal pressure invariant across mammals, the set-point shear must be invariant with
      body mass and the vascular mass must scale close to the observed exponent of hepatic blood volume (0.86).
  P2  Flow and pressure are the same variable only at donor outflow resistance. Grafts whose outflow was
      reconstructed must sit below the diagonal (zP/zF < 1) and tolerate high flow.
  P3  Within a cohort, the outcome must increase with zP, with a common slope and a centre-specific level.
  P4  The clinical thresholds of three different fields must coincide once normalised: HVPG 10 mmHg in cirrhosis,
      PVP 15 mmHg at CVP 5 in a graft, and the limit of safe resection all fall at zP = 2.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import indices
from model import scaling as sc

ROOT = os.path.join(os.path.dirname(__file__), '..'); OUT = f'{ROOT}/results'


def p1_allometry():
    """Set-point and vascular-mass exponents under an invariant portal pressure, for the measured flow and lobule
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
    best = t[(t.flow.str.startswith('portal')) & (t.lobule_exponent == 0.10)]
    fit = best.iloc[(best.vascular_mass_exponent - 0.86).abs().argsort()[:3]]
    return dict(max_abs_tau0_exponent=float(t.tau0_exponent.abs().max()),
                tau0_invariant=bool(t.tau0_exponent.abs().max() < 0.1),
                observed_vascular_mass_exponent=0.86,
                best_b=[float(x) for x in fit.b.tolist()],
                vascular_mass_exponent_at_best_b=[float(x) for x in fit.vascular_mass_exponent.tolist()])


def p2_discordance(d):
    """zP/zF by outflow reconstruction: the model predicts < 1 when the outflow was enlarged, ~1 otherwise."""
    b = d.dropna(subset=['zP', 'zF']).copy(); b['ratio'] = b.zP / b.zF
    rec = b[b.outflow.str.startswith('reconstructed')]; std = b[b.outflow == 'standard']
    return dict(groups_with_both=int(len(b)), reconstructed=int(len(rec)), standard=int(len(std)),
                ratio_reconstructed=[float(x) for x in rec.ratio.round(3)],
                ratio_reconstructed_mean=float(rec.ratio.mean()) if len(rec) else None,
                ratio_standard=[float(x) for x in std.ratio.round(3)],
                all_below_one=bool((rec.ratio < 1).all()) if len(rec) else None,
                outcome_reconstructed_pct=[float(x) for x in rec.pct], 
                detail=[dict(study=r.study, group=r.group, zP=float(r.zP), zF=round(float(r.zF), 2),
                             ratio=round(float(r.ratio), 2), outcome_pct=float(r.pct), outflow=r.outflow)
                        for _, r in b.iterrows()])


def p3_within_cohort(d):
    """Common slope with a centre-specific level, and its internal-external validation by centre."""
    c = indices.contributing(d[d.in_main_analysis & d.zP.notna() & d.events.notna()])
    W = indices.within_study_fit(c, boot=2000)
    H = indices.bayes_hierarchical(c, primary_only=False)
    V = indices.iecv(c)
    pairs = []
    for st, g in c.groupby(indices.strata(c).values):
        g = g.sort_values('zP')
        pairs.append(dict(stratum=st, higher_zP_worse=bool(g.pct.iloc[-1] >= g.pct.iloc[0])))
    return dict(slope=W['beta'], slope_ci=W['beta_ci'], OR_per_zP=W['OR_per_zP'], OR_per_zP_ci=W['OR_per_zP_ci'],
                OR_per_mmHg=W['OR_per_mmHg'], groups=W['groups'], events=W['events'],
                intercept_range=[float(min(W['alphas'].values())), float(max(W['alphas'].values()))],
                hierarchical_mu=H['mu_beta'], hierarchical_ci=H['mu_beta_ci'], prob_positive=H['prob_mu_positive'],
                every_pair_in_predicted_direction=bool(all(p['higher_zP_worse'] for p in pairs)), pairs=pairs,
                iecv=[dict(centre=v['held_out'], slope_from_others=v['beta_from_others'],
                           slope_in_centre=v['beta_in_held_out'], OE=v['OE']) for v in V])


def p4_thresholds():
    """Three thresholds from three fields, on the model's scale."""
    z = indices.zP
    return dict(clinically_significant_portal_hypertension=dict(source='HVPG >= 10 mmHg (Baveno)', zP=z(gradient=10)),
                graft_threshold=dict(source='PVP 15 mmHg at CVP 5 (ILTS 2023)', zP=z(15, 5)),
                variceal_bleeding=dict(source='HVPG >= 12 mmHg', zP=z(gradient=12)),
                resection_limit=dict(source='HVPG 10 mmHg, any resection puts the remnant above', zP=z(gradient=10)),
                all_at_two=bool(abs(z(gradient=10) - 2) < 1e-9 and abs(z(15, 5) - 2) < 1e-9))


def main():
    os.makedirs(OUT, exist_ok=True)
    d = indices.compute(pd.read_csv(f'{ROOT}/data/series.tsv', sep='\t'))
    res = dict(P1_allometry=p1_allometry(), P2_discordance=p2_discordance(d),
               P3_within_cohort=p3_within_cohort(d), P4_thresholds=p4_thresholds())
    json.dump(res, open(f'{OUT}/predictions.json', 'w'), indent=1)
    a, b2, c3, d4 = res['P1_allometry'], res['P2_discordance'], res['P3_within_cohort'], res['P4_thresholds']
    print('P1  comparative anatomy')
    print(f"    shear set-point exponent with body mass: |max| = {a['max_abs_tau0_exponent']:.3f} -> invariant: {a['tau0_invariant']}")
    print(f"    vascular mass exponent closest to the observed 0.86 at b = {a['best_b']} ({[round(x,2) for x in a['vascular_mass_exponent_at_best_b']]})")
    print('P2  flow-pressure discordance')
    print(f"    {b2['reconstructed']} groups with reconstructed outflow, zP/zF = {[round(x,2) for x in b2['ratio_reconstructed']]} "
          f"(all below 1: {b2['all_below_one']}), outcomes {b2['outcome_reconstructed_pct']}%")
    print('P3  within-cohort ordering')
    print(f"    slope {c3['slope']:.2f} (95% CI {c3['slope_ci'][0]:.2f}-{c3['slope_ci'][1]:.2f}), OR {c3['OR_per_zP']:.2f} per unit zP, "
          f"{c3['OR_per_mmHg']:.2f} per mmHg; hierarchical {c3['hierarchical_mu']:.2f} "
          f"({c3['hierarchical_ci'][0]:.2f} to {c3['hierarchical_ci'][1]:.2f}), P(>0) = {c3['prob_positive']:.2f}")
    print(f"    every pair in the predicted direction: {c3['every_pair_in_predicted_direction']}; "
          f"intercepts span {c3['intercept_range'][0]:.2f} to {c3['intercept_range'][1]:.2f} logits")
    for v in c3['iecv']: print(f"    leave out {v['centre']:10s} slope {v['slope_from_others']:.2f} vs {v['slope_in_centre']:.2f} in that centre, O/E {v['OE']:.2f}")
    print('P4  thresholds of three fields on one scale')
    for k, v in d4.items():
        if isinstance(v, dict): print(f"    {v['source']:45s} zP = {v['zP']:.2f}")


if __name__ == '__main__':
    main()
