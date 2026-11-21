"""Checks that matter: the indices are what we say they are, the fitted effect is the published one,
the data cannot drift out of step with the provenance labels, and the calculator matches the analysis.

    python -m pytest tests -q        or        python tests/test_indices.py
"""
import json, os, re, subprocess, sys
import numpy as np, pandas as pd
ROOT = os.path.join(os.path.dirname(__file__), '..'); sys.path.insert(0, os.path.join(ROOT, 'src'))
import indices

def test_sources_parse_on_older_python():
    """Every source file must parse under the oldest Python the CI uses (3.9 grammar), so nested-quote f-strings
    and other 3.12-only syntax cannot slip in."""
    import ast, glob
    files = glob.glob(os.path.join(ROOT, 'src', '**', '*.py'), recursive=True) + [__file__]
    for f in files:
        try:
            ast.parse(open(f, encoding='utf-8').read(), feature_version=(3, 9))
        except SyntaxError as e:
            raise AssertionError(f'{os.path.relpath(f, ROOT)}:{e.lineno} is not valid on Python 3.9: {e.msg}')

def test_load_identity():
    """zR is defined from the resistances, zR = R_graft / R_donor with R = dP / Q per 100 g. With resistances
    estimated from the measured gradient and flow this must reproduce zP / zF, so that zP = zF * zR is an exact
    decomposition of the pressure load (a check of the definitions, not of the data)."""
    import pandas as pd, numpy as np
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    both = d.dropna(subset=['zP', 'zF', 'zR'])
    assert len(both) == 6
    assert np.allclose(both.zR, both.R_graft / both.R_donor), 'zR must be the resistance ratio'
    assert np.allclose(both.zR, both.zP / both.zF), 'the resistance ratio must equal zP / zF'
    assert np.allclose(both.zF * both.zR, both.zP), 'zP = zF * zR must hold exactly'
    Rd = indices.R(indices.NORMAL_GRADIENT, indices.DEFAULT_DONOR_REF)
    assert abs(indices.zR(Rd, Rd) - 1.0) < 1e-12, 'a graft at donor resistance has zR = 1'
    assert all(both.zR < 1) and both.cohort_id.nunique() == 5, 'six groups from five cohorts: two are Wang 2014'
    assert both.configuration.value_counts().to_dict() == {'outflow reconstruction': 3, 'combined': 1, 'inflow modulation': 1, 'none': 1}
    w = d[(d.study == 'Wang 2014') & (d.group == 'splenectomy')].iloc[0]; assert w.outflow == 'standard' and w.configuration == 'inflow modulation', 'a splenectomy is not an outflow reconstruction'
    assert set(d.configuration) == {'outflow reconstruction', 'inflow modulation', 'portosystemic shunt', 'combined', 'none', 'not separable'}

def test_indices_and_baseline():
    assert indices.zP(15, 5) == 2.0                    # PVP 15 at CVP 5 is twice the normal gradient
    assert indices.zP(gradient=10) == 2.0              # a reported gradient is used directly
    assert abs(indices.zF(250, 120) - 2.083) < 1e-3
    a = indices.baseline_from_rate(0.10, 1.6, 1.66)
    assert abs(1 / (1 + np.exp(-(a + 1.66 * 1.6))) - 0.10) < 1e-9
    rng = np.random.default_rng(0); z = rng.normal(1.6, 0.5, 400)
    r = indices.recalibrate(rng.binomial(1, 1 / (1 + np.exp(-(a + 1.66 * z)))), z, 1.66)
    assert abs(r['alpha'] - a) < 3 * r['se']

def test_data_and_provenance():
    d = pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t')
    c = indices.compute(d)
    for col in ('PVP', 'CVP'):                          # a basis must not claim a value that is not there
        assert not (d[col].isna() & ~d[f'{col}_basis'].isin(['not applicable', 'imputed'])).any(), col
    assert len(c) == 32 and c.zP.notna().sum() == 22 and c.zF.notna().sum() == 15
    assert set(c[~c.in_main_analysis].study) == {'Vasavada 2014', 'Yao 2018', 'Kanetkar 2017', 'Wang 2014'}
    m = indices.contributing(c[c.in_main_analysis & c.zP.notna() & c.events.notna()])
    assert m[m.study == 'Wang 2014'].n.sum() == 276, 'a cohort partitioned twice must count once'
    for col, bad in (('CVP_basis', 'reported'), ('gradient_basis', 'not applicable')):
        e = d.copy(); i = e.index[e.study.eq('Osman 2017' if col == 'CVP_basis' else 'Botha 2010')][0]
        e.loc[i, 'CVP' if col == 'CVP_basis' else col] = float('nan') if col == 'CVP_basis' else bad
        if col == 'CVP_basis': e.loc[i, col] = bad
        try:
            indices.compute(e); assert False, f'a label contradicting the data must raise ({col})'
        except ValueError:
            pass

def test_fitted_effect():
    for name in ('cvp_scenarios', 'centred', 'meta_slope', 'within_study_fit'):
        assert hasattr(indices, name), f'src/indices.py is out of date: no {name}()'
    indices.main()                                      # full run, so the JSONs match the committed calculator
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    M = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))
    assert abs(W['beta'] - 1.96) < 0.05 and W['groups'] == 11 and W['events'] == 104
    assert all(v > 1 for v in list(W['leave_one_study_out'].values()) + list(W['by_outcome'].values()))
    assert abs(M['meta']['beta'] - 1.90) < 0.05 and M['meta']['I2'] < 5 and M['meta']['Q_p'] > 0.5
    assert M['meta']['beta_ci'][1] - M['meta']['beta_ci'][0] > M['meta']['beta_ci_normal'][1] - M['meta']['beta_ci_normal'][0]  # Hartung-Knapp
    assert abs(M['centred']['beta'] - 1.66) < 0.05 and M['centred']['r2_weighted'] > 0.8   # descriptive view
    assert M['overdispersion']['se_scale'] == 1.0 and M['attenuation']['loss_pct'] < 10
    assert 'beta_individual_implied' not in M['attenuation']          # no per-patient extrapolation

def test_cvp_scenarios_and_single_source():
    """The assumed CVP must not move the slope, and there must be exactly one hierarchical implementation."""
    import json
    f = os.path.join(ROOT, 'results', 'cvp_scenarios.json')
    if not os.path.exists(f): indices.main()
    betas = {round(x['beta'], 3) for x in json.load(open(f))['cvp_scenarios']}
    assert len(betas) == 1, 'a constant CVP shift within a study must be absorbed by its intercept'
    assert not hasattr(indices, 'bayes_hierarchical'), 'the superseded sampler must not coexist with bayes.fit'
    assert not hasattr(indices, 'iecv'), 'the O/E validation was replaced by bayes.conditional_iecv'

def test_calculator_javascript_parses():
    """The calculator is one HTML file with inline JavaScript: a stray apostrophe breaks the whole page silently.
    Parse it with node when available."""
    import re, shutil, subprocess, tempfile
    html = open(os.path.join(ROOT, 'sfss_calculator.html'), encoding='utf-8').read()
    js = re.search(r'<script>(.*)</script>', html, re.S)
    assert js, 'the calculator must contain its script'
    node = shutil.which('node')
    if not node:
        return                                              # nothing to check with
    with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False, encoding='utf-8') as f:
        f.write(js.group(1)); path = f.name
    r = subprocess.run([node, '--check', path], capture_output=True, text=True)
    assert r.returncode == 0, 'the calculator JavaScript does not parse:\n' + r.stderr[-600:]

def test_the_calculator_gives_the_resistance_load_of_the_model():
    """The resistance load the page shows for a flow load and a CVP is the one src/model gives: the rigid network with
    two inflows (P2) and the three passive relations, for the six groups that report both loads. Run with node."""
    import re, shutil, tempfile
    node = shutil.which('node')
    if not node: return
    assert subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True).returncode == 0
    html = open(os.path.join(ROOT, 'sfss_calculator.html'), encoding='utf-8').read()
    model = html[html.index('// MODEL BEGIN'):html.index('// MODEL END')]
    from model import distension as ds
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t')); B = ds.bracket(d)
    G = {(g['study'], g['group']): g for g in ds.groups(d)}
    cases = [dict(zF=r['zF'], cvps=[G[(r['study'], r['group'])]['cvp']] if G[(r['study'], r['group'])]['cvp'] is not None else list(ds.CVPS)) for r in B['groups']]
    with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as f:
        f.write(model + '\nconst C=' + json.dumps([dict(zF=G[(r['study'], r['group'])]['zF'], cvps=c['cvps']) for r, c in zip(B['groups'], cases)]) +
                ';\nconsole.log(JSON.stringify(C.map(c=>({rigid:rigidRange(c.zF,false), alone:rigidRange(c.zF,true), p:passiveRanges(c.zF,c.cvps)}))));')
    out = json.loads(subprocess.run([node, f.name], capture_output=True, text=True, check=True).stdout)
    for r, o in zip(B['groups'], out):
        lo, hi = ds.rigid_range(G[(r['study'], r['group'])]['zF'])
        assert abs(o['alone'][0] - lo) < 1e-9 and abs(o['alone'][1] - hi) < 1e-9
        assert all(abs(a - b) < 6e-4 for a, b in zip(o['p']['keep'], r['zR_keeps_distending'])) and all(abs(a - b) < 6e-4 for a, b in zip(o['p']['stop'], r['zR_stops_distending'])), (r['study'], o['p'])
        lo2, hi2 = ds.rigid_range(G[(r['study'], r['group'])]['zF'], alone=False); assert abs(o['rigid'][0] - lo2) < 1e-9 and abs(o['rigid'][1] - hi2) < 1e-9
    # with a haematocrit: the factor of model/viscosity.py and the relations solved on the network with that blood
    from model import viscosity as vs
    with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as f:
        f.write(model + '\nconst C=' + json.dumps([dict(zF=G[(r['study'], r['group'])]['zF'], cvps=c['cvps']) for r, c in zip(B['groups'], cases)]) +
                ';\nconsole.log(JSON.stringify({v:[0.25,0.30,0.35].map(viscFactor), p:C.map(c=>passiveRanges(c.zF,c.cvps,viscFactor(0.30))), thr:NET.thr}));')
    out = json.loads(subprocess.run([node, f.name], capture_output=True, text=True, check=True).stdout)
    for H, v in zip((0.25, 0.30, 0.35), out['v']):
        ref = [vs.factor(H, Hr, b) for Hr in vs.HCT_REFERENCE for b in vs.B_HUMAN]; assert abs(v[0] - min(ref)) < 1e-9 and abs(v[1] - max(ref)) < 1e-9
    v30 = out['v'][1]
    for r, o in zip(B['groups'], out['p']):
        e = vs._regions(G[(r['study'], r['group'])], v30); assert all(abs(a - b) < 1e-3 for a, b in zip(o['keep'], e['keeps_distending'])) and all(abs(a - b) < 1e-3 for a, b in zip(o['stop'], e['stops_distending'])), (r['study'], o, e)
    assert abs(out['thr']['zR'] - vs.thresholds()['zR_where_they_agree']) < 1e-3 and out['thr']['flow'] == 250 and out['thr']['gradient'] == 10


def test_plots_and_calculator():
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    names = ['nomogram', 'within_study', 'forest', 'centred', 'risk_change', 'series_pressure', 'series_flow']
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'plots.py'), *names], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert all(os.path.exists(os.path.join(ROOT, 'results', f'{n}.png')) for n in names)
    assert subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True).returncode == 0
    html = open(os.path.join(ROOT, 'sfss_calculator.html')).read()
    B = json.load(open(os.path.join(ROOT, 'results', 'bayes_main.json')))['primary (SFSS or early dysfunction)']
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert '__' not in html.replace('__proto__', '')
    assert f"const BETA={B['mu_beta']:.4f}" in html, 'the calculator must use the primary analysis of the paper'
    opts_ = re.findall(r'<option value="([^"]+)"', html)
    assert all('(mortality' not in o for o in opts_), 'reference series must match the outcome the slope was fitted on'
    assert f"BLO={B['mu_beta_ci'][0]:.4f}" in html, 'and its credible interval'
    assert 'id="pvp2"' in html and 'id="ref"' in html and 'id="basegrad"' in html
    assert 'id="copyrow"' in html and 'caseRow' in html, 'the demo needs the case-collection card'
    assert 'Nothing left this page' in html, 'and must say the data stay in the browser'
    assert 'cohort_template.tsv' in html and 'caseRow' in html, 'the export must follow the cohort template'
    assert 'let LAST' in html and 'calc();' in html.split('function caseRow')[1][:120], \
        'the export must reuse the values the calculator validated, not re-read the fields'
    cols = open(os.path.join(ROOT, 'data', 'cohort_template.tsv')).readline().rstrip('\n').split('\t')
    for c in ('SFSS_ILTS2023', 'graft_loss_90d', 'death_90d'):
        assert c in cols and f"// {c}" in html, f'{c} must be exported in its own column'
    body = html.split('function caseRow')[1]
    assert body.count("flag('o_") == 3, 'each outcome must be flagged independently'
    assert 'const ZGRID=' in html and '"band"' in html and 'slope only' in html   # joint band for published strata

def test_whole_graft_model():
    """The whole-graft limit must reduce to the partial-graft formula and separate the two anastomoses."""
    from model.whole_graft import whole_graft as wg
    n = wg()
    assert abs(n['zF'] - 1) < 1e-9 and abs(n['zP'] - 1) < 1e-9, 'an intact graft at donor inflow must give 1, 1'
    for g, h in ((0.4, 1.5), (0.6, 1.0), (1.0, 2.0)):
        assert abs(wg(h=h / g, r_coll=1e9)['zF'] - h / g) < 1e-9, 'zF must be h/g'
    for d in (0.8, 0.6, 0.5):
        o, i = wg(d_out=d, h=1.5), wg(d_in=d, h=1.5)
        assert o['zP_over_zF'] > 1.05, 'outlet stenosis must raise pressure above flow'
        assert abs(i['zP_over_zF'] - 1) < 1e-9, 'inlet stenosis is upstream of the lobule: the ratio stays at 1'
        assert i['zP'] < 1.5 and i['portal_trunk_pressure'] > 1.5, 'and it lowers graft pressure while raising the trunk'
        assert abs(o['zP_over_zF'] - (0.9 + 0.05 * d ** -4) / 0.95) < 1e-9, 'the ratio is the resistance ratio'
    assert wg(d_in=0.5, h=1.5, r_coll=1e9)['collateral_steal'] < 1e-9, 'no collaterals, no steal'
    assert wg(d_in=0.5, h=1.5)['collateral_steal'] > 0.02, 'with collaterals, an inlet stenosis diverts flow'

def test_modules_are_in_step():
    """Every function the analysis calls must exist, so a half-updated checkout fails here with a clear message
    instead of crashing inside a subprocess."""
    import bayes
    for f in ('fit', 'index_contrast', 'with_cutoff_groups', 'conditional_iecv', 'spec_curve', 'prior_sensitivity',
              'monte_carlo_measurement', 'recovery', 'overlap_sets', 'main'):
        assert hasattr(bayes, f), f'src/bayes.py is out of date: no bayes.{f}()'
    for f in ('zP', 'zF', 'zR', 'compute', 'recalibrate', 'baseline_from_rate'):
        assert hasattr(indices, f), f'src/indices.py is out of date: no indices.{f}()'
    for name in ('bayes_main', 'bayes_extra', 'spec_curve', 'bayes_iecv', 'bayes_priors'):
        assert os.path.exists(os.path.join(ROOT, 'results', f'{name}.json')), f'results/{name}.json is missing'

def test_model_predictions():
    """The four predictions must hold with the committed data."""
    import json, subprocess
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'predictions.py')], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    P = json.load(open(os.path.join(ROOT, 'results', 'predictions.json')))
    assert P['P1_allometry']['tau0_invariant'], 'P1: the scale of the design must not scale with body mass'
    assert 'best_b' not in P['P1_allometry'], 'P1 must not select a value of b that fits the observation'
    p2 = P['P2_discordance']
    # the cohort is the unit: five cohorts, exact sign test in the direction the network gives
    assert (p2['groups_with_both'], p2['cohorts'], p2['cohorts_below_one']) == (6, 5, 5)
    assert abs(p2['sign_test_cohorts_one_sided_p'] - 0.5 ** 5) < 1e-12 and abs(p2['sign_test_cohorts_two_sided_p'] - 2 * 0.5 ** 5) < 1e-12
    assert p2['groups_with_interval_below_one'] == 6, 'at any correlation between the two means from -0.5 to +0.5'
    assert all(x['ratio_ci_by_correlation']['+0.5'][1] < x['ratio_ci'][1] < x['ratio_ci_by_correlation']['-0.5'][1] for x in p2['detail']), 'a positive correlation narrows the ratio'
    # the six groups are five units: no arm of Wang 2014 and no single cohort carries the sign test
    ss = p2['sign_test_sensitivity']
    assert (ss['unit_centre']['below_one'], ss['unit_centre']['of']) == (5, 5) and all(x['below_one'] == x['of'] == 5 for x in ss['one_group_per_cohort']) and len(ss['one_group_per_cohort']) == 2
    assert len(ss['leave_one_cohort_out']) == 5 and all(abs(x['one_sided_p'] - 0.0625) < 1e-12 for x in ss['leave_one_cohort_out']) and ss['smallest_one_sided_p_attainable'] == 0.5 ** 5
    # the uncertainty of each zR by source: sampling (dispersions and correlation swept) apart from the references
    u = p2['uncertainty']; G = {(g['study'], g['group']): g for g in u['groups']}
    assert u['grids']['flow_cv'][0] < 0.33 and u['grids']['flow_cv'][-1] > 0.56 and u['grids']['gradient_sd_mmHg'][-1] > 4.5, 'the sweep must go beyond the reported dispersions on both sides'
    assert sum(g['sweep']['combinations'] for g in u['groups']) == 141 and all(g['sweep']['upper_limit_below_one'] == g['sweep']['combinations'] for g in u['groups']), 'sampling error alone brings no group to one'
    assert all(g['ci_sampling_widest'][0] <= g['ci_sampling'][0] and g['ci_sampling'][1] <= g['ci_sampling_widest'][1] and g['ci_total'][1] >= g['ci_sampling'][1] and g['ci_total'][1] >= g['ci_references'][1] for g in u['groups'])
    assert [g['sweep']['combinations'] for g in u['groups'] if g['flow_cv_reported'] and g['gradient_sd_reported']] == [3, 3, 3], 'a reported dispersion is not swept, only the correlation'
    assert u['groups_with_interval_below_one'] == dict(ci_sampling=6, ci_sampling_widest=6, ci_references=6, ci_total=5, ci_total_widest=5)
    assert u['cohorts_with_interval_below_one']['ci_total'] == 4 and u['groups_with_interval_below_prediction']['ci_sampling'] == 5 and u['groups_with_interval_below_prediction']['ci_total'] == 3
    assert G[('Ishizaki 2012', 'left lobe, no inflow modulation (gradient measured)')]['ci_total'][1] > 1 > G[('Ishizaki 2012', 'left lobe, no inflow modulation (gradient measured)')]['ci_sampling_widest'][1]
    assert all(g['share_of_variance_from_references'] > 0.5 for g in u['groups']) and max(u['groups'], key=lambda g: g['share_of_variance_from_references'])['cvp_assumed'], 'the references are the larger source, most where the CVP is assumed'
    assert p2['by_configuration']['none'] == [0.716] and p2['no_manoeuvre_inside_the_range_of_outflow_reconstruction'], 'configuration does not order the groups'
    wc = p2['within_cohort_contrasts']; assert len(wc) == 1 and wc[0]['cohort'] == 'Wang 2014' and wc[0]['ratio_ci'][1] < 1
    # the prediction of the network uses no pressure of a recipient, and it is too high in every group
    m = p2['model_against_observed']
    assert 0.80 < m['lowest_predicted'] < 0.82 and m['groups_below_the_prediction'] == 6 and m['groups_whose_interval_is_below_the_prediction'] == 5
    assert m['share_of_the_departure_predicted'][1] < 0.8 and 'rmse' not in m, 'no error is computed against the midpoint of the range'
    pf = p2['performance']                                      # four statements about the error, none against the midpoint
    assert pf['central_configuration']['network'] == 'human b=0.75' and pf['central_configuration']['arterial_load'] == 1.0
    assert pf['central']['mean_error'] < -0.25 and 0.30 < pf['central']['rmse'] < 0.36 and all(e < 0 for e in pf['central']['errors'])
    assert pf['best_admissible']['rmse'] < pf['central']['rmse'] and pf['best_admissible']['mean_error'] < -0.2, 'no single configuration of the rigid network reaches the groups'
    assert pf['coverage'] == dict(observed_inside_the_predicted_interval=0, observed_interval_overlaps_it=3, of=6)
    assert pf['references']['combinations'] == 60 and pf['references']['mean_error'][1] < 0
    assert all(x['zR_predicted'][0] > x['zR_observed'] for x in m['groups'])
    sw = p2['reference_sweep']                                  # the normalisation does not manufacture zR < 1
    assert sw['groups'] == 6 and sw['combinations'] == 15
    j = sw['joint']                                             # the three references together: where each statement stops holding
    assert (j['combinations'], j['all_below_one'], j['all_below_model']) == (60, 54, 34)
    assert j['by_normal_gradient']['4 mmHg']['all_below_model'] == 0 and j['by_normal_gradient']['6 mmHg']['all_below_model'] == 20, 'being below the network depends on the donor gradient'
    assert sum(v['share_below_model'] == 1.0 for v in j['per_group'].values()) == 2 and all(x['normal_gradient'] == 4.0 for x in j['not_all_below_one'])
    assert sw['combinations_all_below_one'] == 14 and len(sw['exceptions']) == 1
    assert (sw['exceptions'][0]['donor_ref_flow'], sw['exceptions'][0]['normal_gradient']) == (100.0, 4.0)
    assert 1.0 < sw['exceptions'][0]['max_zR'] < 1.05, 'the one exception is the extreme corner, barely above one'
    assert sw['cvp_sweep_all_below_one'], 'the assumed CVP alone never lifts a group above one'
    p3 = P['P3_within_cohort']
    assert not p3['interval_excludes_zero_under_every_prior'], 'under HalfNormal(1.0) on tau the interval of the mean slope reaches zero'
    assert p3['slope_in_a_new_centre']['primary (SFSS or early dysfunction)']['interval'][0] < 0 < p3['slope_in_a_new_centre']['exploratory (all outcomes)']['interval'][0]
    sc = p3['same_cohort_contrast']; assert len(sc) == 1 and sc[0]['worse_group_has_higher_pressure_load'] and not sc[0]['worse_group_has_higher_flow_load'] and sc[0]['fisher_two_sided_p'] < 0.05
    assert 'pressure_carries_the_association' not in p3['index_contrast'], 'two arms from different cohorts are not a comparison of the indices'
    assert 'ecological' in p3['reading'] and 'exploratory' in P['P4_thresholds']['reading']
    p4 = P['P4_thresholds']; assert p4['title'].startswith('Exploratory comparison') and p4['agreement_is_by_construction'] and set(p4['thresholds_mmHg'].values()) == {10, 12}, 'P4: each threshold is about 10 mmHg before any normalisation'
    rec = p3['design_recovery']; assert all(r['reps'] >= 20 and r['coverage_ci'][0] <= r['coverage'] <= r['coverage_ci'][1] and r['power_ci'][0] <= r['power_to_exclude_zero'] <= r['power_ci'][1] and 'stratum' in r['intercepts'] for r in rec), 'recovery: Monte Carlo error reported, one level per stratum'
    assert p3['every_pair_in_predicted_direction'], 'P3: every within-study pair in the predicted direction'
    for k in ('primary', 'secondary', 'exploratory'):           # the prespecified order
        r = p3[k]
        assert r['mu_beta'] > 0 and r['prob_mu_positive'] > 0.95, f'P3 {k}'
        assert r['rhat']['max'] < 1.01 and r['divergences'] == 0, f'P3 {k}: sampler diagnostics'
        assert r['ess']['min'] > 400, f'P3 {k}: effective sample size'
        assert r['ppc_inside'] == r['ppc_total'], f'P3 {k}: posterior predictive check'
    assert p3['all_specifications_positive'] and p3['mu_range_across_specifications'][0] > 0, 'P3: specification curve'
    ic = p3['index_contrast']                                   # same model on each index, matched on outcome
    assert ic['matched_on_outcome'], 'the two arms must be restricted to the same outcome definition'
    assert ic['zP, pressure per lobule']['mu_beta_ci'][0] > 0, 'the pressure interval must exclude zero'
    assert 'prob_pressure_slope_exceeds_flow_slope' not in ic and ic['slopes_compared'] is False, 'arms from different cohorts that share one series: no probability that one slope exceeds the other'
    assert ic['zP, pressure per lobule']['is_the_primary_fit'] and ic['zP, pressure per lobule']['mu_beta'] == p3['primary']['mu_beta'] and ic['zP, pressure per lobule']['mu_beta_ci'] == p3['primary']['mu_beta_ci'], 'one primary estimate, not two runs of the sampler'
    assert ic['zF, flow per lobule']['events'] < ic['zP, pressure per lobule']['events'], 'and its weakness stated'
    a = p3['accounting']                                        # the material must be reported in full
    assert a['series'] == 18 and a['groups'] == 32 and a['series_with_a_within_study_contrast'] == 5
    assert p3['with_cutoff_groups']['mu_beta_ci'][0] > 0, 'P3: adding cut-off groups must not overturn it'
    assert all(v['directions_correct'] and v['log_score_gain'] > 0 for v in p3['conditional_iecv']), 'P3: conditional validation'
    assert all(pr['prob_positive'] > 0.9 for pr in p3['prior_sensitivity']), 'P3: prior sensitivity'
    assert P['P4_thresholds']['all_at_two'], 'P4'

BLOCKED_MACROS = ('operatorname', 'newcommand', 'def', 'require', 'label', 'includegraphics')


def test_documents_are_the_right_documents():
    """Three documents, each with its own content: the README runs the code, THEORY.md states the model, the changelog
    lists the versions; docs/ holds only the release manifest."""
    read = lambda f: open(os.path.join(ROOT, f), encoding='utf-8').read()
    assert read('README.md').startswith('# LiPNet') and '## Running' in read('README.md') and 'THEORY.md' in read('README.md')
    assert read('THEORY.md').startswith('# LiPNet: theory') and 'z_P = z_F\\,z_R' in read('THEORY.md')
    assert read('CHANGELOG.md').startswith('# Changelog')
    version = re.search(r'^version: (\S+)', read('CITATION.cff'), re.M).group(1)
    assert f'## [{version}]' in read('CHANGELOG.md') and f'version-{version}-' in read('README.md'), 'version in CITATION.cff, README and changelog'
    assert sorted(os.listdir(os.path.join(ROOT, 'docs'))) == ['MANIFEST-release.json']
    for f in ('src', 'tests', 'run_all.sh', 'README.md', 'THEORY.md'):
        p = os.path.join(ROOT, f); paths = [os.path.join(r, n) for r, _, ns in os.walk(p) for n in ns if n.endswith(('.py', '.md', '.sh'))] if os.path.isdir(p) else [p]
        for q in paths:
            assert not re.search(r'docs/(?!MANIFEST-release)\w+\.md', open(q, encoding='utf-8').read()), f'{q} points to a document that is gone'


def test_docs_math_renders_on_github():
    """GitHub renders README maths with a restricted KaTeX: a blocked macro silently kills the whole expression."""
    docs = {f: open(os.path.join(ROOT, f), encoding='utf-8').read() for f in ('README.md', 'THEORY.md')}
    for name, s in docs.items():
        exprs = re.findall(r'\$\$(.+?)\$\$', s, re.S) + re.findall(r'(?<!\$)\$([^$\n]+)\$(?!\$)', s)
        for e in exprs:
            for m in BLOCKED_MACROS:
                assert '\\' + m not in e, f'{name}: \\{m} is not allowed by GitHub maths: {e.strip()[:60]}'
            assert '\\\\' not in e, f'{name}: a double backslash breaks the expression: {e.strip()[:60]}'


if __name__ == '__main__':
    for t in (test_sources_parse_on_older_python, test_load_identity, test_indices_and_baseline, test_data_and_provenance, test_fitted_effect, test_cvp_scenarios_and_single_source, test_calculator_javascript_parses, test_plots_and_calculator, test_whole_graft_model, test_modules_are_in_step, test_model_predictions, test_documents_are_the_right_documents, test_docs_math_renders_on_github):
        t(); print('ok', t.__name__)
