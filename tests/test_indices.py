"""Minimal checks: index definitions, reproducibility of the pooled fit, and that every plot renders.
    python -m pytest tests -q          (or)          python tests/test_indices.py
"""
import os, sys, json, subprocess
ROOT = os.path.join(os.path.dirname(__file__), '..'); sys.path.insert(0, os.path.join(ROOT, 'src'))
import indices

def test_series_indices_match_stored_results():
    import pandas as pd
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    assert d.zP.notna().sum() == 21 and d.zF.notna().sum() == 14 and len(d) == 31
    assert (d.gradient_mmHg.notna() & d.gradient_basis.eq('not applicable')).sum() == 0
    assert set(d.CVP_basis.unique()) <= {'reported', 'derived', 'imputed', 'threshold', 'not applicable'}

def test_incomplete_edit_is_caught():
    """A provenance label that contradicts the data must raise, not silently impute."""
    import pandas as pd
    d = pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t')
    i = d.index[d.study.eq('Osman 2017')][0]                   # a row whose zP comes from PVP and CVP
    e = d.copy(); e.loc[i, 'CVP'] = float('nan'); e.loc[i, 'CVP_basis'] = 'reported'
    try:
        indices.compute(e); assert False, 'a blank CVP labelled as reported must raise'
    except ValueError:
        pass
    e.loc[i, 'CVP_basis'] = 'imputed'; c = indices.compute(e)
    assert bool(c.loc[i, 'CVP_filled']) and bool(c.loc[i, 'zP_imputed'])
    j = d.index[d.study.eq('Botha 2010')][0]                   # a row whose zP comes from a reported gradient
    e2 = d.copy(); e2.loc[j, 'gradient_basis'] = 'not applicable'
    try:
        indices.compute(e2); assert False, 'a gradient labelled not applicable must raise'
    except ValueError:
        pass

def test_cut_off_groups_excluded_from_main_analysis():
    import pandas as pd
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    assert set(d[~d.in_main_analysis].study) == {'Vasavada 2014', 'Yao 2018', 'Kanetkar 2017', 'Wang 2014'}
    # a cohort partitioned twice enters the model once
    m = indices.contributing(d[d.in_main_analysis & d.zP.notna() & d.events.notna()])
    assert m[m.study == 'Wang 2014'].n.sum() == 276, 'the Wang cohort must not be counted twice'

def test_index_definitions():
    assert indices.zP(15, 5) == 2.0
    assert indices.zP(gradient=10) == 2.0            # a reported gradient is used directly            # consensus threshold: PVP 15 at CVP 5
    assert indices.zP(10, 5) == 1.0            # normal gradient
    assert abs(indices.zF(250, 120) - 2.083) < 1e-3   # 250 mL/min/100 g with the donor reference of Troisi 2005
    assert indices.zF(90, 90) == 1.0

def test_within_study_fit_reproduces_published_effect():
    indices.main()                              # full run (seed 0), so the JSONs match the committed calculator
    L = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert abs(L['beta'] - 1.96) < 0.03, L['beta']
    assert abs(L['OR_per_mmHg'] - 1.48) < 0.03
    assert L['groups'] == 11 and len(L['studies']) == 5 and L['patients'] == 734 and L['events'] == 104
    assert all(v > 1.0 for v in L['leave_one_study_out'].values()), 'slope must stay positive dropping any study'
    assert all(v > 1.0 for v in L['by_outcome'].values()), 'slope must hold for each outcome definition'
    assert L['bootstrap_failures'] == 0, L['bootstrap_failures']
    lo, hi = L['profile_ci']                                   # profile likelihood must agree with the bootstrap
    assert abs(lo - L['beta_ci'][0]) < 0.3 and abs(hi - L['beta_ci'][1]) < 0.4
    assert all(' / ' in s for s in L['studies']), 'strata are study x outcome'

def test_risk_after_is_an_odds_shift():
    assert abs(indices.risk_after(0.10, -1.0, 1.96) - 0.0154) < 0.001
    assert indices.risk_after(0.10, 0.0, 1.96) == 0.10

def test_meta_analysis_and_diagnostics():
    """Random-effects pooling, overdispersion and regression dilution must agree with the stratified fit."""
    import json
    M = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))
    meta, od, att = M['meta'], M['overdispersion'], M['attenuation']
    assert abs(meta['beta'] - 1.90) < 0.05 and abs(meta['OR_per_zP'] - 6.69) < 0.3
    assert meta['tau2'] == 0 and meta['I2'] < 5 and meta['Q_p'] > 0.5     # the strata estimate one common effect
    assert len(meta['strata']) == 5 and all(r['beta'] > 0 for r in meta['strata'])
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert abs(meta['beta'] - W['beta']) < 0.2, 'meta-analysis and stratified fit must agree'
    assert od['ratio'] < 1.5 and od['se_scale'] == 1.0                    # no extra-binomial spread
    assert att['loss_pct'] < 10 and att['lambda_individual'] < 1 and att['beta_individual_implied'] > W['beta']

def test_centring_removes_the_cohort_level():
    """Centring within stratum must cancel the intercepts exactly and recover the same slope."""
    import json, pandas as pd, numpy as np
    C = json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))['centred']
    P = pd.DataFrame(C['points'])
    for _, g in P.groupby('stratum'):                      # each stratum is centred: its means are zero
        assert abs(g.zP_c.mean()) < 1e-9 and abs(g.log_odds_c.mean()) < 1e-9
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert abs(C['beta'] - W['beta']) < 0.4, 'centred slope must agree with the stratified fit'
    assert C['r2_weighted'] > 0.8, C['r2_weighted']        # one common line once the cohort level is gone
    assert C['OR_per_zP_ci'][0] > 1

def test_fast_plots_render():
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'plots.py'), 'nomogram', 'within_study', 'forest', 'centred', 'risk_change', 'series_pressure', 'series_flow', 'plane', 'resection', 'load_curve', 'interventions'], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for n in ('nomogram', 'within_study', 'forest', 'centred', 'risk_change', 'series_pressure'):
        assert os.path.exists(os.path.join(ROOT, 'results', f'{n}.png'))

def test_readme_numbers_match_results():
    """The README quotes the main correlations and the fit; they must match results/ so that the two cannot drift apart."""
    import pandas as pd, json
    indices.main()                              # full run: never leave a short-bootstrap JSON behind for build_app
    sens = pd.read_csv(os.path.join(ROOT, 'results', 'sensitivity.tsv'), sep='\t')
    readme = open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read()
    main = sens[sens.analysis.str.startswith('main analysis')]
    for _, r in main.iterrows():
        assert f"{int(r.groups)} groups" in readme, f"README lacks the group count for main {r['index']}"
        assert f"rho = {r.spearman_rho:.2f}" in readme, f"README lacks rho for main {r['index']}"
        assert f"p = {r.p:.3f}" in readme, f"README lacks p for main {r['index']}"
    W = json.load(open(os.path.join(ROOT, 'results', 'within_study_fit.json')))
    assert f"{W['OR_per_zP']:.2f}" in readme and f"{W['OR_per_mmHg']:.2f}" in readme, 'README odds ratios differ from the fit'
    assert f"{W['beta']:.2f}" in readme and f"{W['beta_ci'][0]:.2f}" in readme, 'README slope differs from the fit'
    assert '−6.41' not in readme and '−6.40' not in readme, 'stale pooled intercept in README'

def test_app_builds():
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    html = open(os.path.join(ROOT, 'sfss_calculator.html')).read()
    assert '__BETA__' not in html and 'const BETA=' in html
    assert 'id="gnorm"' not in html            # the normal gradient is fixed by the fitted model, not user-editable
    assert 'cvpRaw===\'\'' in html or "cvpRaw===''" in html   # blank CVP is distinguished from an explicit zero
    assert 'function tagP(' in html and 'function tagF(' in html   # pressure and flow have separate labels
    assert 'must be a positive number' in html         # flow and weight inputs are validated
    assert 'does <b>not</b> estimate' in html          # no absolute-risk claim without a baseline rate
    assert 'const BETA=' in html and 'const B0=' not in html
    assert 'id="pvp2"' in html and 'id="base"' in html and 'id="ref"' in html and 'id="basegrad"' in html
    import json as _json
    C = _json.load(open(os.path.join(ROOT, 'results', 'meta_slope.json')))['centred']
    assert f"const BETA={C['beta']:.4f}" in html, 'the calculator must use the intercept-free centred slope'
    import re as _re
    opts = _re.findall(r'<option value="([^"]+)"', html)
    assert opts[0] == 'cohort' and opts[1] == 'own' and len(opts) >= 5, opts     # cohort average, own rate, fitted strata
    assert '"zmin"' in html and 'extrapolation of its intercept' in html        # reference range is checked
    assert 'const CENTRED=' in html and 'anchored on ' in html                  # one clean curve, anchored or centred
    assert 'STUDIES.forEach((st,i)=>{let pts' not in html                       # no five-curve plot any more

if __name__ == '__main__':
    for t in (test_series_indices_match_stored_results, test_incomplete_edit_is_caught, test_cut_off_groups_excluded_from_main_analysis, test_index_definitions, test_within_study_fit_reproduces_published_effect, test_meta_analysis_and_diagnostics, test_centring_removes_the_cohort_level, test_risk_after_is_an_odds_shift, test_fast_plots_render, test_readme_numbers_match_results, test_app_builds): t(); print('ok', t.__name__)
