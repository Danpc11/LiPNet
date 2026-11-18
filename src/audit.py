"""What the extraction table and the fits support, stated without the model's narrative.

    python src/audit.py      -> results/audit.json

Four blocks, each a set of facts about data/series.tsv and the cached fits:

  table      arithmetic and provenance of the table itself
  paired     the five groups that report both loads: what each ingredient of zR is, how uncertain the ratio is with
             the dispersions the papers report, and which reference values would move it across one
  loads      how much of the pressure load is readable from the flow load, with no reference value involved
  slope      what each stratum of the outcome model contributes: how the groups are defined, how far apart they are
             in the predictor, and the slope each gives on its own
  model      the model's own interval for zR of a partial graft (src/model/liver.py) next to the observed values
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy.stats import norm, spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indices

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

# Dispersions the papers report for the five paired groups, with where each comes from. cv_flow is SD / mean of the
# portal flow per 100 g; sd_grad the SD of the portocaval gradient in mmHg; n_flow and n_grad the recipients measured.
REPORTED = {
    ('Chan 2011', 'right lobe with MHV, no modulation'): dict(
        cv_flow=None, flow_range=(127.0, 754.0), sd_grad=None, grad_range=(-1.0, 11.0), n_flow=46, n_grad=46,
        basis='ranges of Table 3 (data/patient_level/chan2011_groups.tsv), SD by the range estimator of Wan et al.'),
    ('Yagi 2006', 'PVP kept <20'): dict(
        cv_flow=592.0 / 1828.0, sd_grad=None, n_flow=27, n_grad=27,
        basis='PVF 1828 +- 592 mL/min (notes); no dispersion of the portal pressure in the table'),
    ('Wang 2014', 'splenectomy'): dict(cv_flow=None, sd_grad=None, n_flow=154, n_grad=154, basis='no dispersion in the table'),
    ('Wang 2014', 'no splenectomy'): dict(cv_flow=None, sd_grad=None, n_flow=122, n_grad=122, basis='no dispersion in the table'),
    ('Ishizaki 2012', 'left lobe, no inflow modulation (gradient measured)'): dict(
        cv_flow=167.0 / 301.0, sd_grad=4.4, n_flow=38, n_grad=31,
        basis='PVF 301 +- 167 (n = 38), gradient 12.4 +- 4.4 (n = 26 to 31 of 42) (notes)'),
}


def sd_from_range(lo, hi, n):
    """Wan et al., BMC Med Res Methodol 2014: SD from the range of n observations."""
    return (hi - lo) / (2 * norm.ppf((n - 0.375) / (n + 0.25)))


def table_checks(d):
    pct = 100 * d.events / d.n
    t = d[d.PVP.notna() & d.CVP.notna() & d.gradient_mmHg.notna()]
    return dict(rows=int(len(d)), series=int(d.study.nunique()), centres=int(d.centre_id.nunique()),
                pct_disagrees_with_events_over_n=int(((d.pct - pct).abs() > 0.5).sum()),
                gradient_disagrees_with_PVP_minus_CVP=int(((t.PVP - t.CVP - t.gradient_mmHg).abs() > 0.05).sum()),
                groups_with_pressure_load=int(d.zP.notna().sum()), groups_with_flow_load=int(d.zF.notna().sum()),
                groups_with_both=int((d.zP.notna() & d.zF.notna()).sum()),
                CVP_basis=d[d.zP.notna()].CVP_basis.value_counts().to_dict(),
                donor_flow_reference_measured={s: float(v) for s, v in d[['study', 'donor_PVF_per_100g_ref']].drop_duplicates().dropna().values if v != indices.DEFAULT_DONOR_REF},
                donor_flow_reference_default=float(indices.DEFAULT_DONOR_REF),
                series_with_a_measured_donor_gradient=0, normal_gradient_assumed_mmHg=float(indices.NORMAL_GRADIENT),
                columns_recording_when_flow_was_measured=[c for c in d.columns if 'flow' in c.lower() and 'timing' in c.lower()],
                gradient_timing_of_paired_groups=d.dropna(subset=['zF', 'zP']).set_index(d.dropna(subset=['zF', 'zP']).study + ' / ' + d.dropna(subset=['zF', 'zP']).group).gradient_timing.to_dict())


def paired(d):
    p = d.dropna(subset=['zF', 'zP']).copy(); rows = []
    cvs = [v['cv_flow'] for v in REPORTED.values() if v['cv_flow']] + [sd_from_range(*REPORTED[('Chan 2011', 'right lobe with MHV, no modulation')]['flow_range'], 46) / 318.0]
    sds = [v['sd_grad'] for v in REPORTED.values() if v['sd_grad']] + [sd_from_range(*REPORTED[('Chan 2011', 'right lobe with MHV, no modulation')]['grad_range'], 46)]
    cv_max, sd_max = max(cvs), max(sds)                       # used where a paper reports no dispersion
    measured_refs = sorted(v for v in d.donor_PVF_per_100g_ref.dropna().unique() if v != indices.DEFAULT_DONOR_REF)
    for _, r in p.iterrows():
        rep = REPORTED[(r.study, r.group)]; grad = 5 * r.zP
        cv = rep['cv_flow'] or (sd_from_range(*rep['flow_range'], rep['n_flow']) / r.PVF_per_100g if rep.get('flow_range') else cv_max)
        sd = rep['sd_grad'] or (sd_from_range(*rep['grad_range'], rep['n_grad']) if rep.get('grad_range') else sd_max)
        se = np.sqrt((cv / np.sqrt(rep['n_flow'])) ** 2 + (sd / np.sqrt(rep['n_grad']) / grad) ** 2)      # relative SE of the ratio of the two means
        own_ref = bool(pd.notna(r.donor_PVF_per_100g_ref) and r.donor_PVF_per_100g_ref != indices.DEFAULT_DONOR_REF)
        ref = r.donor_PVF_per_100g_ref if pd.notna(r.donor_PVF_per_100g_ref) else indices.DEFAULT_DONOR_REF
        rows.append(dict(
            study=r.study, group=r.group, n=int(r.n), zF=float(r.zF), zP=float(r.zP), zR=float(r.zR),
            flow_basis=r.PVF_basis, donor_flow_reference=float(ref), donor_flow_reference_measured_in_series=own_ref,
            pressure_basis=r.PVP_basis, CVP_basis=r.CVP_basis, gradient_basis=r.gradient_basis, gradient_timing=r.gradient_timing,
            flow_cv_used=float(cv), flow_cv_reported=bool(rep['cv_flow'] or rep.get('flow_range')),
            gradient_sd_used=float(sd), gradient_sd_reported=bool(rep['sd_grad'] or rep.get('grad_range')), dispersion_basis=rep['basis'],
            zR_ci_sampling=[float(r.zR * np.exp(-1.96 * se)), float(r.zR * np.exp(1.96 * se))],
            # reference values that put the group at zR = 1, one at a time
            donor_gradient_at_which_zR_is_one_mmHg=float(5 * r.zR),
            donor_flow_reference_at_which_zR_is_one=(None if own_ref else float(ref / r.zR)),
            zR_at_largest_measured_donor_reference=(None if own_ref else float(r.zR * max(measured_refs) / ref))))
    return dict(groups=rows, flow_cv_where_unreported=float(cv_max), gradient_sd_where_unreported=float(sd_max),
                measured_donor_flow_references=[float(x) for x in measured_refs],
                note='zR_ci_sampling is the sampling error of the two group means only; it holds the donor references fixed')


def loads(d):
    p = d.dropna(subset=['zF', 'zP']); lf, lp, lr = np.log(p.zF), np.log(p.zP), np.log(p.zR)
    vF, vR, cov, vP = lf.var(ddof=1), lr.var(ddof=1), np.cov(lf, lr, ddof=1)[0, 1], lp.var(ddof=1)
    rho, pval = spearmanr(p.zF, p.zP); w = p[p.study == 'Wang 2014'].set_index('group')
    return dict(groups=int(len(p)), zF_range=[float(p.zF.min()), float(p.zF.max())], zF_fold=float(p.zF.max() / p.zF.min()),
                zP_range=[float(p.zP.min()), float(p.zP.max())], zP_fold=float(p.zP.max() / p.zP.min()),
                zR_range=[float(p.zR.min()), float(p.zR.max())], zR_fold=float(p.zR.max() / p.zR.min()),
                spearman_zF_zP=float(rho), spearman_p=float(pval),
                share_of_var_log_zP=dict(flow=float(vF / vP), resistance=float(vR / vP), covariance=float(2 * cov / vP)),
                within_cohort_wang_2014=dict(flow_change_pct=float(100 * (w.PVF_per_100g['splenectomy'] / w.PVF_per_100g['no splenectomy'] - 1)),
                                             portal_pressure_change_pct=float(100 * (w.PVP['splenectomy'] / w.PVP['no splenectomy'] - 1)),
                                             outcome_pct=[float(w.pct['no splenectomy']), float(w.pct['splenectomy'])]),
                note='ratios and ranks between groups do not depend on the normal gradient; they depend on the donor flow reference only through the series that report their own')


def slope(d):
    m = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]; c = indices.contributing(m); meta = indices.meta_slope(c)
    own = {r['stratum']: r for r in meta['strata']}; rows = []
    for st, g in c.groupby(indices.strata(c).values):
        sd = g.gradient_sd.dropna()
        rows.append(dict(stratum=st, groups=int(len(g)), defined_by=g.group.tolist(), events=int(g.events.sum()), recipients=int(g.n.sum()),
                         span_mmHg=float(5 * (g.zP.max() - g.zP.min())), within_group_sd_mmHg=([float(sd.min()), float(sd.max())] if len(sd) else None),
                         span_over_within_group_sd=(float(5 * (g.zP.max() - g.zP.min()) / sd.mean()) if len(sd) else None),
                         own_slope=float(own[st]['beta']), own_slope_se=float(own[st]['se']),
                         groups_defined_by_pressure=bool(g.group.str.contains('PVP|PVPG|PSG').all())))
    out = dict(strata=rows)
    f = os.path.join(ROOT, 'results', 'bayes_overlap.json')
    if os.path.exists(f):
        for r in json.load(open(f)):
            if sorted(s.split(' / ')[0] for s in r['strata']) == ['Osman 2017', 'Wang 2014']:
                out['primary_without_uemura'] = dict(mu_beta=r['mu_beta'], mu_beta_ci=r['mu_beta_ci'], events=r['events'], recipients=r['patients']); break
    f = os.path.join(ROOT, 'results', 'bayes_main.json')
    if os.path.exists(f):
        p = json.load(open(f))['primary (SFSS or early dysfunction)']
        out['primary'] = dict(mu_beta=p['mu_beta'], mu_beta_ci=p['mu_beta_ci'], events=p['events'], recipients=p['patients'], prediction_interval_new_study=p['prediction_interval'])
    f = os.path.join(ROOT, 'results', 'bayes_extra.json')
    if os.path.exists(f):
        ic = json.load(open(f))['index_contrast']
        out['pressure_against_flow'] = {k: dict(groups=v['groups'], series=len(v['strata']), events=v['events'], median_group_size=v['median_group_size'],
                                                mu_beta=v['mu_beta'], mu_beta_ci=v['mu_beta_ci']) for k, v in ic.items() if isinstance(v, dict)}
        out['pressure_against_flow']['prob_pressure_slope_exceeds_flow_slope'] = ic['prob_pressure_slope_exceeds_flow_slope']
        out['pressure_against_flow']['series_in_both_arms'] = ic['shared_series']
    return out


def model_against_data(d):
    from model import liver
    p = d.dropna(subset=['zF', 'zP']); out = dict(observed_zR=[float(x) for x in p.zR])
    lob = [liver.lobar_graft(liver.HUMAN_GENERATIONS, b, 1) for b in (1.0, 0.75, 0.5)]
    out['network_lobar_graft'] = [float(min(l['zR'] for l in lob)), float(max(l['zR'] for l in lob))]
    two = [liver.two_inflows(zF, l['share_inlet'], l['zR'], arterial_load=phi)['zR'] for l in lob for zF in (p.zF.min(), p.zF.max()) for phi in (0.0, 1.0)]
    out['network_with_arterial_inflow'] = [float(min(two)), float(max(two))]
    out['network_with_arterial_inflow_any_split'] = [float(liver.two_inflows(float(p.zF.max()), 0.0, arterial_load=0.0)['zR']), 1.0]
    f = os.path.join(ROOT, 'results', 'model_liver.json')
    if os.path.exists(f) and 'graft_on_explicit_graph_range' in json.load(open(f)): out['explicit_graphs'] = json.load(open(f))['graft_on_explicit_graph_range']
    f = os.path.join(ROOT, 'results', 'fem_network_link.json')
    if os.path.exists(f): out['with_shear_feedback_of_the_lobule_layer'] = [float(g['zR_predicted']) for g in json.load(open(f))['against_the_published_groups']['groups']]
    lo = min(out['network_lobar_graft'][0], out['network_with_arterial_inflow'][0], *(out.get('explicit_graphs') or [1.0]))
    out['lowest_model_value'] = float(lo); out['observed_below_lowest_model_value'] = int(sum(z < lo for z in out['observed_zR']))
    return out


def main():
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    out = dict(table=table_checks(d), paired=paired(d), loads=loads(d), slope=slope(d), model=model_against_data(d))
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True); json.dump(out, open(os.path.join(ROOT, 'results', 'audit.json'), 'w'), indent=1)
    t = out['table']
    print(f"table: {t['rows']} groups, {t['series']} series; pct against events/n: {t['pct_disagrees_with_events_over_n']} disagreements; gradient against PVP - CVP: {t['gradient_disagrees_with_PVP_minus_CVP']}")
    print(f"   donor flow reference measured in {len(t['donor_flow_reference_measured'])} series {t['donor_flow_reference_measured']}, default {t['donor_flow_reference_default']:.0f} elsewhere; "
          f"donor gradient measured in {t['series_with_a_measured_donor_gradient']} series (assumed {t['normal_gradient_assumed_mmHg']:.0f} mmHg); timing of the flow measurement recorded: {bool(t['columns_recording_when_flow_was_measured'])}")
    print('paired groups: zR, sampling interval of the ratio, and the reference that would put it at one')
    for g in out['paired']['groups']:
        print(f"   {g['study']:14s} {g['group'][:16]:16s} zR {g['zR']:.2f} ({g['zR_ci_sampling'][0]:.2f} to {g['zR_ci_sampling'][1]:.2f}; flow CV {g['flow_cv_used']:.2f}{'' if g['flow_cv_reported'] else '*'}, gradient SD {g['gradient_sd_used']:.1f}{'' if g['gradient_sd_reported'] else '*'}) "
              f"| at one if donor gradient <= {g['donor_gradient_at_which_zR_is_one_mmHg']:.1f} mmHg" + ('' if g['donor_flow_reference_measured_in_series'] else f", or donor flow >= {g['donor_flow_reference_at_which_zR_is_one']:.0f}; at reference {max(out['paired']['measured_donor_flow_references']):.0f}: {g['zR_at_largest_measured_donor_reference']:.2f}"))
    print('   * not reported by the paper: the largest reported value is used')
    l = out['loads']; s = l['share_of_var_log_zP']
    print(f"loads: zF {l['zF_fold']:.2f}-fold, zP {l['zP_fold']:.2f}-fold, zR {l['zR_fold']:.2f}-fold; Spearman zF-zP {l['spearman_zF_zP']:+.2f}; var(log zP): flow {100*s['flow']:.0f}%, resistance {100*s['resistance']:.0f}%, covariance {100*s['covariance']:+.0f}%")
    print('slope model, stratum by stratum')
    for r in out['slope']['strata']:
        print(f"   {r['stratum']:42s} span {r['span_mmHg']:.2f} mmHg, own slope {r['own_slope']:+.2f} (SE {r['own_slope_se']:.2f}), {r['events']}/{r['recipients']}, groups defined by pressure: {r['groups_defined_by_pressure']}")
    if 'primary_without_uemura' in out['slope']:
        a, b = out['slope']['primary'], out['slope']['primary_without_uemura']
        print(f"   primary {a['mu_beta']:+.2f} ({a['mu_beta_ci'][0]:+.2f} to {a['mu_beta_ci'][1]:+.2f}), {a['events']}/{a['recipients']}; without Uemura 2016 {b['mu_beta']:+.2f} ({b['mu_beta_ci'][0]:+.2f} to {b['mu_beta_ci'][1]:+.2f}), {b['events']}/{b['recipients']}")
    m = out['model']
    print(f"model against data: network {m['network_lobar_graft'][0]:.2f} to {m['network_lobar_graft'][1]:.2f}; with the arterial inflow {m['network_with_arterial_inflow'][0]:.2f} to {m['network_with_arterial_inflow'][1]:.2f}"
          + (f"; explicit graphs {m['explicit_graphs'][0]:.2f} to {m['explicit_graphs'][1]:.2f}" if 'explicit_graphs' in m else '')
          + f"; observed {min(m['observed_zR']):.2f} to {max(m['observed_zR']):.2f}; {m['observed_below_lowest_model_value']} of {len(m['observed_zR'])} below the lowest model value {m['lowest_model_value']:.2f}")


if __name__ == '__main__':
    main()
