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
