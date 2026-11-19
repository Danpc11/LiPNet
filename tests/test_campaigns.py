"""The cached tables of the model are the ones the published code writes, and the manifest check knows about them
and about the figures."""
import json, os, subprocess, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))


def test_the_cached_model_tables_are_reproduced_by_the_code():
    """The two smallest sizes of each campaign and the exact certificate, recomputed and compared (about 20 s).
    python src/model/verify_campaigns.py does every size."""
    from model import verify_campaigns as vc
    f = os.path.join(ROOT, 'results', 'campaign_check.json'); keep = open(f).read() if os.path.exists(f) else None
    try:
        r = vc.main(quick=True)
    finally:
        if keep is not None: open(f, 'w').write(keep)             # the record of the full check is not replaced by the quick one
    assert r['agrees'] and vc.TOLERANCE <= 1e-6
    for name, t in r['tables'].items():
        assert t['rows_recomputed'] >= 3 and not t['disagreements'] and t['largest_relative_difference'] <= vc.TOLERANCE, name
    assert r['tables']['certificate_small']['rows_recomputed'] == r['tables']['certificate_small']['rows_cached'] == 3


def test_a_changed_number_in_a_cached_table_is_caught():
    import pandas as pd
    from model import verify_campaigns as vc
    old = pd.read_csv(os.path.join(ROOT, 'data', 'graph_scaling_all.tsv'), sep='\t'); new = old.copy(); new.loc[0, 'D_in'] *= 1 + 1e-6
    c = vc.compare('graph_scaling_all', new, vc.TABLES['graph_scaling_all'])
    assert not c['agrees'] and c['disagreements'][0]['column'] == 'D_in'
    assert vc.compare('graph_scaling_all', old, vc.TABLES['graph_scaling_all'])['agrees']


def test_every_figure_of_the_article_declares_what_it_reads():
    import paper_figures as pf
    assert set(pf.READS) == set(pf.FIGURES)
    for k, paths in pf.READS.items():
        assert paths and all(os.path.exists(os.path.join(ROOT, p)) or p.startswith('results/') for p in paths), k
    f = os.path.join(ROOT, 'results', 'figures_provenance.json')
    if os.path.exists(f):                                           # written when the figures are drawn
        rec = json.load(open(f)); assert set(rec) == {f'Figure_{k}' for k in pf.FIGURES}
        import manifest
        n, bad = manifest.figures(); assert n == 10 and not bad, bad
