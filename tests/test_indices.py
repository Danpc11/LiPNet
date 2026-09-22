"""Checks that matter: the indices are what we say they are, the fitted effect is the published one,
the data cannot drift out of step with the provenance labels, and the calculator matches the analysis.

    python -m pytest tests -q        or        python tests/test_indices.py
"""
import json, os, subprocess, sys
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
    for name in ('bayes_hierarchical', 'iecv', 'cvp_scenarios', 'centred', 'meta_slope'):
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

def test_hierarchical_and_validation():
    """The hierarchical model, the leave-one-centre-out validation and the CVP scenarios must run and agree in sign."""
    import json
    assert hasattr(indices, 'bayes_hierarchical'), 'src/indices.py is out of date: no bayes_hierarchical()'
    f = os.path.join(ROOT, 'results', 'hierarchical.json')
    if not os.path.exists(f): indices.main()            # tests must not depend on the order they run in
    H = json.load(open(f))
    b, ba = H['bayes']['primary'], H['bayes']['all']
    assert b['groups'] == 7 and len(b['studies']) == 3 and ba['groups'] == 11 and len(ba['studies']) == 5
    assert 0.5 < b['mu_beta'] < 1.2 and b['mu_beta_ci'][0] < 0 < b['mu_beta_ci'][1]    # honest uncertainty
    assert ba['prob_mu_positive'] > 0.95 and ba['tau_beta'] > 0
    centres = {r['held_out'] for r in H['iecv']}
    assert centres == {'Cairo', 'Fukuoka', 'Kyoto'} and all(r['beta_from_others'] > 1 for r in H['iecv'])
    assert all(abs(r['OE'] - 1) < 0.05 for r in H['iecv'])                              # intercept recalibration works
    betas = {round(x['beta'], 3) for x in H['cvp_scenarios']}
    assert len(betas) == 1, 'a constant CVP shift must be absorbed by the intercept'

def test_plots_and_calculator():
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    names = ['nomogram', 'within_study', 'forest', 'centred', 'risk_change', 'series_pressure', 'series_flow']
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'plots.py'), *names], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert all(os.path.exists(os.path.join(ROOT, 'results', f'{n}.png')) for n in names)
    assert subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True).returncode == 0
    html = open(os.path.join(ROOT, 'sfss_calculator.html')).read()
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert '__' not in html.replace('__proto__', '')
    assert f"const BETA={W['beta']:.4f}" in html, 'slope and intercepts must come from the same model'
    assert f"BLO={W['beta_ci'][0]:.4f}" in html and f"{W['alphas'][W['studies'][0]]:.4f}" in html
    assert 'id="pvp2"' in html and 'id="ref"' in html and 'id="basegrad"' in html
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

def test_model_predictions():
    """The four predictions must hold with the committed data."""
    import json, subprocess
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'predictions.py')], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    P = json.load(open(os.path.join(ROOT, 'results', 'predictions.json')))
    assert P['P1_allometry']['tau0_invariant'], 'P1: the shear set-point must not scale with body mass'
    assert P['P2_discordance']['all_below_one'] and P['P2_discordance']['reconstructed'] >= 4, 'P2'
    assert max(P['P2_discordance']['outcome_reconstructed_pct']) <= 10, 'P2: those grafts did well'
    assert P['P3_within_cohort']['every_pair_in_predicted_direction'], 'P3'
    assert P['P3_within_cohort']['prob_positive'] > 0.95
    assert P['P4_thresholds']['all_at_two'], 'P4'

def test_readme_matches_results():
    sens = pd.read_csv(os.path.join(ROOT, 'results', 'sensitivity.tsv'), sep='\t')
    M = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))
    readme = open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read()
    for _, r in sens[sens.analysis.str.startswith('main analysis')].iterrows():
        assert f"{int(r.groups)} groups" in readme and f"ρ = {r.spearman_rho:.2f}" in readme and f"p = {r.p:.3f}" in readme
    assert 'P1.' in readme and 'P2.' in readme and 'P3.' in readme and 'P4.' in readme, 'the README must state the predictions'
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert f"{W['OR_per_zP']:.2f}" in readme and f"{W['OR_per_mmHg']:.2f}" in readme
    assert 'assets/graphical_abstract.png' in readme

if __name__ == '__main__':
    for t in (test_sources_parse_on_older_python, test_indices_and_baseline, test_data_and_provenance, test_fitted_effect, test_hierarchical_and_validation,
              test_plots_and_calculator, test_whole_graft_model, test_model_predictions, test_readme_matches_results):
        t(); print('ok', t.__name__)
