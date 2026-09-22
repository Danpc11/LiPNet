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
    """The three LiPNet loads and their identity: pressure load = flow load x resistance load."""
    import pandas as pd, numpy as np
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    both = d.dropna(subset=['zP', 'zF', 'zR'])
    assert len(both) == 5 and np.allclose(both.zF * both.zR, both.zP), 'zP = zF * zR must hold exactly'
    assert abs(indices.zR(2.0, 2.0) - 1.0) < 1e-12, 'a donor-resistance graft has zR = 1'
    assert all(both.zR < 1), 'every graft measured on both scales here had its outflow enlarged'

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
    assert len(c) == 31 and c.zP.notna().sum() == 21 and c.zF.notna().sum() == 14
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

def test_plots_and_calculator():
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    names = ['nomogram', 'within_study', 'forest', 'centred', 'risk_change', 'series_pressure', 'series_flow']
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'plots.py'), *names], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert all(os.path.exists(os.path.join(ROOT, 'results', f'{n}.png')) for n in names)
    assert subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True).returncode == 0
    html = open(os.path.join(ROOT, 'sfss_calculator.html')).read()
    B = json.load(open(os.path.join(ROOT, 'results', 'bayes_main.json')))['exploratory (all outcomes)']
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert '__' not in html.replace('__proto__', '')
    assert f"const BETA={B['mu_beta']:.4f}" in html, 'the calculator must use the primary hierarchical estimate'
    assert f"BLO={B['mu_beta_ci'][0]:.4f}" in html, 'and its credible interval'
    assert f"{W['alphas'][W['studies'][0]]:.4f}" in html, 'with the per-series levels for the reference option'
    assert 'id="pvp2"' in html and 'id="ref"' in html and 'id="basegrad"' in html
    assert 'id="copyrow"' in html and 'caseRow' in html, 'the demo needs the case-collection card'
    assert 'never left this page' in html, 'and must say the data stay in the browser'
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
    assert P['P1_allometry']['tau0_invariant'], 'P1: the shear set-point must not scale with body mass'
    assert 'best_b' not in P['P1_allometry'], 'P1 must not select a value of b that fits the observation'
    p2 = P['P2_discordance']
    assert p2['all_below_one'] and p2['reconstructed'] >= 4, 'P2: every reconstructed outflow below one'
    assert p2['min_prob_below_one'] > 0.9, 'P2: with measurement uncertainty'
    assert max(p2['outcome_reconstructed_pct']) <= 10, 'P2: those grafts did well (reported separately)'
    p3 = P['P3_within_cohort']
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
    assert 'prob_pressure_slope_exceeds_flow_slope' in ic, 'the comparison must be reported as a probability'
    assert ic['zF, flow per lobule']['events'] < ic['zP, pressure per lobule']['events'], 'and its weakness stated'
    a = p3['accounting']                                        # the material must be reported in full
    assert a['series'] == 17 and a['groups'] == 31 and a['series_with_a_within_study_contrast'] == 5
    assert p3['with_cutoff_groups']['mu_beta_ci'][0] > 0, 'P3: adding cut-off groups must not overturn it'
    assert all(v['directions_correct'] and v['log_score_gain'] > 0 for v in p3['conditional_iecv']), 'P3: conditional validation'
    assert all(pr['prob_positive'] > 0.9 for pr in p3['prior_sensitivity']), 'P3: prior sensitivity'
    assert P['P4_thresholds']['all_at_two'], 'P4'

def test_readme_matches_results():
    sens = pd.read_csv(os.path.join(ROOT, 'results', 'sensitivity.tsv'), sep='\t')
    M = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))
    readme = open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read()
    for _, r in sens[sens.analysis.str.startswith('main analysis')].iterrows():
        flat = re.sub(r'\s|\\\\[,;: ]|\\\\rho|ρ', '', readme)     # the README may write these in prose or in LaTeX
        assert f"{int(r.groups)}groups" in flat and f"={r.spearman_rho:.2f}" in flat and f"p={r.p:.3f}" in flat
    assert 'P1.' in readme and 'P2.' in readme and 'P3.' in readme and 'P4.' in readme, 'the README must state the predictions'
    E = json.load(open(os.path.join(ROOT, 'results', 'bayes_extra.json')))['index_contrast']
    assert f"{E['prob_pressure_slope_exceeds_flow_slope']:.2f}" in readme, 'the README must report the index comparison'
    A = json.load(open(os.path.join(ROOT, 'results', 'predictions.json')))['P1_allometry']
    assert f"{100*A['fraction_compatible']:.0f}%" in readme, 'the README must report P1 over the whole range of b'
    assert '0.89–0.90 at b' not in readme, 'stale cherry-picked P1 claim'
    M = json.load(open(os.path.join(ROOT, 'results', 'bayes_main.json')))
    p = M['primary (SFSS or early dysfunction)']
    assert f"{p['mu_beta']:+.2f}" in readme and f"{p['OR_per_zP']:.2f}" in readme, 'README must quote the primary analysis'
    assert f"{p['prob_mu_positive']:.3f}" in readme
    assert 'assets/graphical_abstract.png' in readme
    assert 'LiPNet' in readme, 'the README must carry the model name'
    identity = re.sub(r'\s|\\(?:[,;: ]|cdot|times)', '', readme)   # LaTeX spacing/product must not matter
    assert 'z_P=z_Fz_R' in identity, 'the README must carry the identity zP = zF zR'
    assert 'liver_pressure_index' not in readme, 'stale repository name'

if __name__ == '__main__':
    for t in (test_sources_parse_on_older_python, test_load_identity, test_indices_and_baseline, test_data_and_provenance, test_fitted_effect, test_cvp_scenarios_and_single_source, test_calculator_javascript_parses, test_plots_and_calculator, test_whole_graft_model, test_modules_are_in_step, test_model_predictions, test_readme_matches_results):
        t(); print('ok', t.__name__)
