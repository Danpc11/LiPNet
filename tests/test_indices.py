"""Minimal checks: index definitions, reproducibility of the pooled fit, and that every plot renders.
    python -m pytest tests -q          (or)          python tests/test_indices.py
"""
import os, sys, json, subprocess
ROOT = os.path.join(os.path.dirname(__file__), '..'); sys.path.insert(0, os.path.join(ROOT, 'src'))
import indices

def test_index_definitions():
    assert indices.zP(15, 5) == 2.0            # consensus threshold: PVP 15 at CVP 5
    assert indices.zP(10, 5) == 1.0            # normal gradient
    assert abs(indices.zF(250, 120) - 2.083) < 1e-3   # 250 mL/min/100 g with the donor reference of Troisi 2005
    assert indices.zF(90, 90) == 1.0

def test_pooled_fit_reproduces_published_coefficients():
    indices.main(boot=50)                       # short bootstrap for speed; point estimate is deterministic
    L = json.load(open(os.path.join(ROOT, 'results', 'logit_zP.json')))
    assert abs(L['b0'] - (-6.41)) < 0.02 and abs(L['b1'] - 1.88) < 0.02
    assert L['groups'] == 5 and L['patients'] == 362 and L['events'] == 46

def test_fast_plots_render():
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'src'), MPLBACKEND='Agg')
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'plots.py'), 'nomogram', 'risk_curve', 'series_pressure', 'series_flow', 'plane', 'resection', 'load_curve', 'interventions'], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    for n in ('nomogram', 'risk_curve', 'series_pressure', 'series_flow'):
        assert os.path.exists(os.path.join(ROOT, 'results', f'{n}.png'))

def test_app_builds():
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'src', 'build_app.py')], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    html = open(os.path.join(ROOT, 'sfss_calculator.html')).read()
    assert '__CURVE__' not in html and 'const B0=' in html

if __name__ == '__main__':
    for t in (test_index_definitions, test_pooled_fit_reproduces_published_coefficients, test_fast_plots_render, test_app_builds): t(); print('ok', t.__name__)
