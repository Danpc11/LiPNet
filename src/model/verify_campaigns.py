"""Are the cached tables of the model the ones the published code writes?

    python src/model/verify_campaigns.py            every size of the three tables (about 4 min)
    python src/model/verify_campaigns.py quick      the two smallest sizes of each campaign and the certificate (about 20 s)

data/graph_scaling_all.tsv, data/graph3d_all.tsv and data/certificate_small.tsv are read by the figures and are not
rebuilt by run_all.sh. This script runs the code that produced them (graph_scaling.py, network3d.py, certify.py) in a
scratch directory, merges the output as merge_campaigns.py does, and compares every number with the cached one. The
columns that are results (dissipations, counts, gaps) must agree within TOLERANCE, relative; the run times are not
compared. The verdict goes to results/campaign_check.json with the hash of each table and of the code that was run.
"""
import hashlib, json, os, subprocess, sys, tempfile, time
import numpy as np, pandas as pd

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); ROOT = os.path.dirname(SRC)
TOLERANCE = 1e-9                                   # relative; the optimiser is deterministic at fixed seeds, the rest is floating point
ABSOLUTE = 1e-12                                   # for the columns that are differences from zero (gaps, checks)
CODE = ('model/network.py', 'model/network3d.py', 'model/graph_scaling.py', 'model/certify.py', 'model/merge_campaigns.py')
TABLES = dict(graph_scaling_all=dict(script='model/graph_scaling.py', key=['R', 'b'], near_zero=[]),
              graph3d_all=dict(script='model/network3d.py', key=['R', 'b'], near_zero=['check']),
              certificate_small=dict(script='model/certify.py', key=['R', 'b'], near_zero=['rel_gap']))

sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()


def run(args, cwd):
    env = dict(os.environ, PYTHONPATH=SRC, MPLBACKEND='Agg')
    r = subprocess.run([sys.executable] + args, cwd=cwd, env=env, capture_output=True, text=True)
    if r.returncode: raise RuntimeError(f"{' '.join(args)} failed:\n{r.stderr[-2000:]}")


def compare(name, new, spec):
    old = pd.read_csv(f'{ROOT}/data/{name}.tsv', sep='\t'); key = spec['key']
    old['_k'] = [tuple(np.round(r, 6)) for r in old[key].values]; new = new.copy(); new['_k'] = [tuple(np.round(r, 6)) for r in new[key].values]
    m = old.merge(new, on='_k', suffixes=('', '_new')); worst = 0.0; bad = []
    for c in old.columns:
        if c in key + ['_k', 't']: continue
        a, b = m[c].values.astype(float), m[c + '_new'].values.astype(float)
        both_nan = np.isnan(a) & np.isnan(b)
        if c in spec['near_zero']: err = np.where(both_nan, 0.0, np.abs(a - b)); ok = err <= 1e-9
        else: err = np.where(both_nan, 0.0, np.abs(a - b) / np.maximum(np.abs(a), ABSOLUTE)); ok = err <= TOLERANCE; worst = max(worst, float(np.nanmax(err)))
        for i in np.flatnonzero(~ok): bad.append(dict(column=c, key=[float(x) for x in m['_k'].iloc[i]], cached=float(a[i]), recomputed=float(b[i])))
    return dict(rows_cached=int(len(old)), rows_recomputed=int(len(m)), sizes_recomputed=sorted(float(x) for x in new.R.unique()),
                largest_relative_difference=worst, disagreements=bad, agrees=bool(not bad and len(m) > 0), sha256=sha(f'{ROOT}/data/{name}.tsv'))


def main(quick=False):
    t0 = time.time(); out = {}
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(f'{tmp}/results'); os.makedirs(f'{tmp}/data')
        for name, script in (('graph_scaling_all', 'model/graph_scaling.py'), ('graph3d_all', 'model/network3d.py')):
            sizes = sorted(pd.read_csv(f'{ROOT}/data/{name}.tsv', sep='\t').R.unique())
            for R in (sizes[:2] if quick else sizes): run([f'{SRC}/{script}', str(R)] + (['results'] if 'network3d' in script else []), tmp); print(f'{name} R = {R} done ({time.time() - t0:.0f} s)', flush=True)
        run([f'{SRC}/model/merge_campaigns.py'], tmp); run([f'{SRC}/model/certify.py'], tmp); print(f'certificate done ({time.time() - t0:.0f} s)', flush=True)
        new = dict(graph_scaling_all=pd.read_csv(f'{tmp}/data/graph_scaling_all.tsv', sep='\t'), graph3d_all=pd.read_csv(f'{tmp}/data/graph3d_all.tsv', sep='\t'),
                   certificate_small=pd.read_csv(f'{tmp}/results/certificate_small.tsv', sep='\t'))
    for name, spec in TABLES.items(): out[name] = compare(name, new[name], spec)
    res = dict(agrees=bool(all(v['agrees'] for v in out.values())), every_size=bool(not quick), tolerance_relative=TOLERANCE, tables=out,
               code={c: sha(f'{SRC}/{c}') for c in CODE})
    os.makedirs(f'{ROOT}/results', exist_ok=True); json.dump(res, open(f'{ROOT}/results/campaign_check.json', 'w'), indent=1)
    print(f'{time.time() - t0:.0f} s')
    for name, v in out.items():
        print(f"{name}: {v['rows_recomputed']} of {v['rows_cached']} rows recomputed, largest relative difference {v['largest_relative_difference']:.1e}, "
              f"{'agrees' if v['agrees'] else 'DISAGREES: ' + json.dumps(v['disagreements'][:3])}")
    return res


if __name__ == '__main__':
    sys.exit(0 if main(quick=sys.argv[1:] == ['quick'])['agrees'] else 1)
