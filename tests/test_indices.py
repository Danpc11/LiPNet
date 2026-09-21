"""Minimal checks: index definitions, reproducibility of the pooled fit, and that every plot renders.
    python -m pytest tests -q          (or)          python tests/test_indices.py
"""
import os, sys, json, subprocess
ROOT = os.path.join(os.path.dirname(__file__), '..'); sys.path.insert(0, os.path.join(ROOT, 'src'))
import indices

def test_series_indices_match_stored_results():
    import pandas as pd
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    assert d.zP.notna().sum() == 13 and d.zF.notna().sum() == 13
    assert set(d.CVP_basis.unique()) <= {'reported', 'derived', 'imputed', 'not applicable'}

def test_incomplete_edit_is_caught():
    import pandas as pd
    d = pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'); i = d.index[d.study.eq('Yamada 2008')][0]
    e = d.copy(); e.loc[i, 'CVP'] = float('nan')
    try:
        indices.compute(e); assert False, 'a blank CVP labelled as reported must raise'
    except ValueError:
        pass
    e.loc[i, 'CVP_basis'] = 'imputed'; c = indices.compute(e)
    assert bool(c.loc[i, 'CVP_filled']) and bool(c.loc[i, 'any_imputed'])

def test_cut_off_groups_excluded_from_main_analysis():
    import pandas as pd
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    assert (~d.in_main_analysis).sum() == 2 and set(d[~d.in_main_analysis].study) == {'Vasavada 2014'}

def test_index_definitions():
    assert indices.zP(15, 5) == 2.0            # consensus threshold: PVP 15 at CVP 5
    assert indices.zP(10, 5) == 1.0            # normal gradient
    assert abs(indices.zF(250, 120) - 2.083) < 1e-3   # 250 mL/min/100 g with the donor reference of Troisi 2005
    assert indices.zF(90, 90) == 1.0

def test_pooled_fit_reproduces_published_coefficients():
    indices.main()                              # full run (seed 0), so results/logit_zP.json is byte-identical to the committed one
    L = json.load(open(os.path.join(ROOT, 'results', 'logit_zP.json')))
    assert abs(L['b0'] - (-6.40)) < 0.02 and abs(L['b1'] - 1.88) < 0.02
    assert L['normal_gradient_mmHg'] == 5 and L['zP_min'] > 1.4 and L['zP_max'] < 2.7
    assert L['groups'] == 5 and L['patients'] == 362 and L['events'] == 46

def test_fast_plots_render():
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'plots.py'), 'nomogram', 'risk_curve', 'series_pressure', 'series_flow', 'plane', 'resection', 'load_curve', 'interventions'], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for n in ('nomogram', 'risk_curve', 'series_pressure', 'series_flow'):
        assert os.path.exists(os.path.join(ROOT, 'results', f'{n}.png'))

def test_readme_numbers_match_results():
    """The README quotes the main correlations and the fit; they must match results/ so that the two cannot drift apart."""
    import pandas as pd, json
    indices.main()                              # full run: never leave a short-bootstrap JSON behind for build_app
    sens = pd.read_csv(os.path.join(ROOT, 'results', 'sensitivity.tsv'), sep='\t'); L = json.load(open(os.path.join(ROOT, 'results', 'logit_zP.json')))
    readme = open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read()
    main = sens[sens.analysis.str.startswith('main analysis')]
    for _, r in main.iterrows():
        assert f"{int(r.groups)} groups" in readme and f"ρ = {r.spearman_rho:.2f}" in readme, f"README lacks main {r['index']} result"
        assert (f"p = {r.p:.3f}" in readme) or (f"p = {r.p:.2f}" in readme), f"README lacks p for {r['index']}"
    assert f"{L['b0']:.2f}".replace('-', '−') in readme and f"+ {L['b1']:.2f}" in readme, 'README intercept/slope differ from the fit'
    assert '−6.41' not in readme, 'stale intercept in README'

def test_app_builds():
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    html = open(os.path.join(ROOT, 'sfss_calculator.html')).read()
    assert '__CURVE__' not in html and 'const B0=' in html
    assert 'id="gnorm"' not in html            # the normal gradient is fixed by the fitted model, not user-editable
    assert 'cvpRaw===\'\'' in html or "cvpRaw===''" in html   # blank CVP is distinguished from an explicit zero
    assert 'of the 5 groups had no reported central venous pressure' in html
    assert 'function tagP(' in html and 'function tagF(' in html   # pressure and flow have separate labels
    assert 'must be a positive number' in html         # flow and weight inputs are validated

if __name__ == '__main__':
    for t in (test_series_indices_match_stored_results, test_incomplete_edit_is_caught, test_cut_off_groups_excluded_from_main_analysis, test_index_definitions, test_pooled_fit_reproduces_published_coefficients, test_fast_plots_render, test_readme_numbers_match_results, test_app_builds): t(); print('ok', t.__name__)
