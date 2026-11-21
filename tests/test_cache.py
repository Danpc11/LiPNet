"""Cached fits are tied to the table and the model code that produced them (src/cache.py, src/manifest.py)."""
import json, os, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
import cache, manifest


def test_the_cached_fits_are_those_of_this_table_and_code():
    names = ('bayes_main', 'bayes_iecv', 'spec_curve', 'bayes_priors', 'bayes_measurement', 'bayes_recovery', 'bayes_extra', 'bayes_overlap')
    assert cache.stale(names) == [], 'run python src/bayes.py: the fits on disk were computed from another table or other code'


def test_a_change_of_the_table_or_a_rebuild_makes_the_cache_stale(tmp_path, monkeypatch):
    out = tmp_path / 'results'; out.mkdir(); src = tmp_path / 'src'; src.mkdir(); data = tmp_path / 'data'; data.mkdir()
    for f in cache.INPUTS: (tmp_path / f).write_text('a')
    monkeypatch.setattr(cache, 'ROOT', str(tmp_path)); monkeypatch.setattr(cache, 'OUT', str(out)); monkeypatch.setattr(cache, 'STAMPS', str(out / 'bayes_inputs.json'))
    monkeypatch.delenv('LIPNET_REBUILD', raising=False)
    (out / 'bayes_main.json').write_text('{}'); assert not cache.fresh('bayes_main'), 'a file with no stamp is not trusted'
    cache.stamp('bayes_main'); assert cache.fresh('bayes_main')
    (data / 'series.tsv').write_text('b'); assert not cache.fresh('bayes_main'), 'the table changed'
    cache.stamp('bayes_main'); assert cache.fresh('bayes_main')
    monkeypatch.setenv('LIPNET_REBUILD', '1'); assert not cache.fresh('bayes_main'), 'a rebuild ignores every cache'


def test_the_lock_file_pins_every_package_the_pipeline_imports():
    lock = dict(l.strip().split('==') for l in open(os.path.join(ROOT, 'requirements-lock.txt')) if '==' in l)
    for p in manifest.PACKAGES: assert p in lock, p
    assert all('>' not in l and '<' not in l for l in open(os.path.join(ROOT, 'requirements-lock.txt')) if not l.startswith('#')), 'exact versions only'


def test_manifest_records_inputs_results_and_versions():
    m = manifest.build()
    assert 'data/series.tsv' in m['inputs'] and 'src/bayes.py' in m['inputs'] and m['packages']['numpyro'] and len(m['inputs']) > 30
    assert all(len(h) == 64 for h in m['inputs'].values())


def _same(a, b, tol=2e-3, path=''):
    """Recursive comparison; numbers within tol (the stored file rounds, and percentiles move in the last digit between
    versions of numpy), everything else exactly. Returns the paths that differ."""
    if isinstance(a, dict) and isinstance(b, dict):
        return [f'{path}/{k}' for k in set(a) ^ set(b)] + [x for k in set(a) & set(b) for x in _same(a[k], b[k], tol, f'{path}/{k}')]
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return [f'{path}: length'] if len(a) != len(b) else [x for i, (u, v) in enumerate(zip(a, b)) for x in _same(u, v, tol, f'{path}[{i}]')]
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str) or isinstance(b, str): return [] if a == b else [path]
    return [] if abs(float(a) - float(b)) <= tol * max(1.0, abs(float(a))) else [path]


def test_the_stored_predictions_are_the_ones_this_code_gives():
    """results/predictions.json is kept in the repository. Everything in it that does not come out of a sampler is
    recomputed here from the table and the code and compared, and the sampled part must be the cached fits themselves.
    A stored file that was not regenerated after a change of the code fails this test."""
    import pandas as pd, indices, predictions as P
    stored = json.load(open(os.path.join(ROOT, 'results', 'predictions.json'))); d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    now = json.loads(json.dumps(dict(P1_allometry=P.p1_allometry(), P2_discordance=P.p2_discordance(d), P4_thresholds=P.p4_thresholds()), default=float))
    m = P.p2_model_against_observed(d); now['P2_discordance']['model_against_observed'] = json.loads(json.dumps(m, default=float))
    u = P.p2_uncertainty(d, m); now['P2_discordance']['uncertainty'] = json.loads(json.dumps(u, default=float))
    now['P2_discordance']['performance'] = json.loads(json.dumps(P.p2_performance(d, m, u), default=float))
    now['P2_discordance']['reference_sweep'] = json.loads(json.dumps(P.p2_reference_sweep(d, model_lowest=m['lowest_predicted']), default=float))
    bad = [x for k in now for x in _same(stored[k], now[k], path=k)]
    assert not bad, 'results/predictions.json is not what the code gives; run python src/predictions.py. First differences: ' + ', '.join(bad[:8])
    p2 = stored['P2_discordance']; assert 'sign_test_all_p' not in p2 and 'reconstructed' not in p2, 'the stored file has the structure of an older release'
    p3 = stored['P3_within_cohort']; R = lambda n: json.load(open(os.path.join(ROOT, 'results', n + '.json')))
    main = R('bayes_main'); assert p3['primary'] == main['primary (SFSS or early dysfunction)'] and p3['secondary'] == main['secondary (mortality or graft loss)'] and p3['exploratory'] == main['exploratory (all outcomes)']
    assert p3['specification_curve'] == R('spec_curve') and p3['design_recovery'] == R('bayes_recovery') and p3['index_contrast'] == R('bayes_extra')['index_contrast'] and p3['measurement_monte_carlo'] == R('bayes_measurement')
    assert json.loads(json.dumps(P.same_cohort_contrast(d), default=float)) == p3['same_cohort_contrast']


def test_the_comparison_is_not_fooled_by_an_old_file():
    assert _same(dict(a=1.0, b=[0.5, 'x']), dict(a=1.0005, b=[0.5, 'x'])) == []
    assert _same(dict(cohorts=5, p=0.03125), dict(groups=6, p=0.015625)) != []
    assert _same(dict(p=0.03125), dict(p=0.015625)) == ['/p']


def test_the_documents_quote_the_cohort_level_test_and_an_exploratory_p4():
    """What a reader of THEORY.md and README.md must find, and what must be gone: the sign test over six groups, the
    comparison of thresholds stated as a prediction, and a second version of the primary slope."""
    theory = open(os.path.join(ROOT, 'THEORY.md'), encoding='utf-8').read(); readme = open(os.path.join(ROOT, 'README.md'), encoding='utf-8').read(); both = theory + readme
    assert '0.016' not in both and '0.015625' not in both, 'the sign test is over five cohorts (0.031 one-sided, 0.06 two-sided), not six groups'
    assert '0.031' in theory and r'p=0.063\ \text{two-sided}' in theory and 'five of five' in theory
    assert '### P4. Exploratory comparison of haemodynamic thresholds across liver diseases and surgical settings' in theory
    p4 = theory[theory.index('### P4.'):theory.index('## What LiPNet is not')]
    assert 'not a prediction' in p4 or 'not something the model predicts' in p4
    assert 'does not show that the three settings share a physiological regime' in p4 and 'common haemodynamic regime' not in both and 'LiPNet predicts' not in p4
    main = json.load(open(os.path.join(ROOT, 'results', 'bayes_main.json')))['primary (SFSS or early dysfunction)']
    assert f"{main['mu_beta']:+.2f} ({main['mu_beta_ci'][0]:+.2f}, {main['mu_beta_ci'][1]:+.2f})" in theory and '+1.56 (+0.35' not in both, 'one primary estimate'


def test_the_release_manifest_is_the_code_of_this_tree():
    """docs/MANIFEST-release.json is the record of the clean rebuild of the release. Its fingerprint must be the one of the
    code and data of this tree, file by file; if it names a commit present here, that commit must hold the same files."""
    import manifest, subprocess
    m = json.load(open(os.path.join(ROOT, 'docs', 'MANIFEST-release.json')))
    now = {f: manifest.sha256(os.path.join(ROOT, f)) for f in manifest.files(manifest.INPUTS)}
    assert m['inputs_fingerprint'] == manifest.fingerprint(m['inputs']), 'the fingerprint is the one of the inputs it lists'
    diff = sorted(f for f in set(now) | set(m['inputs']) if now.get(f) != m['inputs'].get(f))
    assert not diff, 'the release manifest is not the rebuild of this code: ' + ', '.join(diff[:8])
    c = m.get('commit')
    if c and manifest._git('cat-file', '-e', f'{c}^{{commit}}') is not None:
        assert manifest.verify_release() == 0
