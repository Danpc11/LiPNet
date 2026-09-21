"""Checks that matter: the indices are what we say they are, the fitted effect is the published one,
the data cannot drift out of step with the provenance labels, and the calculator matches the analysis.

    python -m pytest tests -q        or        python tests/test_indices.py
"""
import json, os, subprocess, sys
import numpy as np, pandas as pd
ROOT = os.path.join(os.path.dirname(__file__), '..'); sys.path.insert(0, os.path.join(ROOT, 'src'))
import indices

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
    indices.main()                                      # full run, so the JSONs match the committed calculator
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    M = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))
    assert abs(W['beta'] - 1.96) < 0.05 and W['groups'] == 11 and W['events'] == 104
    assert all(v > 1 for v in list(W['leave_one_study_out'].values()) + list(W['by_outcome'].values()))
    assert abs(M['meta']['beta'] - 1.90) < 0.05 and M['meta']['I2'] < 5 and M['meta']['Q_p'] > 0.5
    assert abs(M['centred']['beta'] - 1.66) < 0.05 and M['centred']['r2_weighted'] > 0.8
    assert M['overdispersion']['se_scale'] == 1.0 and M['attenuation']['loss_pct'] < 10

def test_plots_and_calculator():
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    names = ['nomogram', 'within_study', 'forest', 'centred', 'risk_change', 'series_pressure', 'series_flow']
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'plots.py'), *names], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert all(os.path.exists(os.path.join(ROOT, 'results', f'{n}.png')) for n in names)
    assert subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True).returncode == 0
    html = open(os.path.join(ROOT, 'sfss_calculator.html')).read()
    C = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))['centred']
    assert '__' not in html.replace('__proto__', '')
    assert f"const BETA={C['beta']:.4f}" in html, 'the calculator must use the intercept-free centred slope'
    assert 'id="pvp2"' in html and 'id="ref"' in html and 'id="basegrad"' in html

def test_readme_matches_results():
    sens = pd.read_csv(os.path.join(ROOT, 'results', 'sensitivity.tsv'), sep='\t')
    M = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))
    readme = open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read()
    for _, r in sens[sens.analysis.str.startswith('main analysis')].iterrows():
        assert f"{int(r.groups)} groups" in readme and f"ρ = {r.spearman_rho:.2f}" in readme and f"p = {r.p:.3f}" in readme
    assert f"{M['centred']['OR_per_zP']:.2f}" in readme and f"{M['centred']['OR_per_mmHg']:.2f}" in readme
    assert 'assets/graphical_abstract.png' in readme

if __name__ == '__main__':
    for t in (test_indices_and_baseline, test_data_and_provenance, test_fitted_effect,
              test_plots_and_calculator, test_readme_matches_results):
        t(); print('ok', t.__name__)
