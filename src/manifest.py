"""A record of what produced the results on disk, and a check that they still match.

    python src/manifest.py            write results/MANIFEST.json
    python src/manifest.py check      exit 1 if an input or a result is no longer the one recorded
    python src/manifest.py compare docs/MANIFEST-release.json      this tree and its results against the clean rebuild of the release

The manifest holds the SHA-256 of every input (data/, src/), of every numeric result (results/*.json, *.tsv), the
versions of Python and of the packages, and the commit if the tree is a git checkout. run_all.sh writes it last;
`bash run_all.sh --clean` deletes every generated result first and ignores every cache, which is the run to do before
quoting a number.

Figures. src/paper_figures.py records, when it draws a figure, the hash of every file the figure reads, of the two
scripts that draw it and of the two images (results/figures_provenance.json). The check compares those with the files
now on disk: a figure whose input, script or image has changed since it was drawn is named. The same goes for the
tables of the model in data/ that no script of run_all.sh rebuilds: results/campaign_check.json
(src/model/verify_campaigns.py) holds the hash of each table and of the code it was recomputed with.
"""
import glob, importlib.metadata as md, json, os, platform, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cache import ROOT, OUT, sha256

PACKAGES = ('numpy', 'scipy', 'pandas', 'matplotlib', 'jax', 'jaxlib', 'numpyro', 'networkx', 'adjustText')


def files(patterns):
    return sorted(os.path.relpath(f, ROOT).replace(os.sep, '/') for p in patterns for f in glob.glob(f'{ROOT}/{p}', recursive=True) if '__pycache__' not in f and os.path.isfile(f))


def build():
    inputs = files(('data/**/*', 'src/**/*.py', 'src/app_template.html', 'run_all.sh', 'requirements*.txt'))
    results = [f for f in files(('results/*.json', 'results/*.tsv')) if not f.endswith('MANIFEST.json')]
    try: commit = subprocess.run(['git', '-C', ROOT, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip() or None
    except OSError: commit = None
    def version(p):
        try: return md.version(p)
        except md.PackageNotFoundError: return None
    return dict(python=platform.python_version(), platform=platform.platform(), commit=commit, rebuilt_from_scratch=bool(os.environ.get('LIPNET_REBUILD') or os.environ.get('LIPNET_REBUILT')),
                packages={p: version(p) for p in PACKAGES}, inputs={f: sha256(f'{ROOT}/{f}') for f in inputs}, results={f: sha256(f'{ROOT}/{f}') for f in results})


def figures():
    """Each figure of the article against the record made when it was drawn. Returns (number of figures, complaints)."""
    f = f'{OUT}/figures_provenance.json'
    if not os.path.exists(f): return 0, ['no results/figures_provenance.json: run python src/paper_figures.py']
    rec = json.load(open(f)); bad = []; now = lambda p: sha256(f'{ROOT}/{p}') if os.path.isfile(f'{ROOT}/{p}') else None
    for name, r in rec.items():
        for kind, what in (('reads', 'input'), ('drawn_by', 'script')):
            for p, h in r[kind].items():
                if now(p) != h: bad.append(f'{name}: {what} {p} is not the one the figure was drawn from')
        for e, h in r['files'].items():
            if now(f'results/{name}.{e}') != h: bad.append(f'{name}.{e} is not the image recorded when the figure was drawn')
    return len(rec), bad


def campaigns():
    """The cached tables of the model against the record of their last recomputation."""
    f = f'{OUT}/campaign_check.json'
    if not os.path.exists(f): return ['no results/campaign_check.json: run python src/model/verify_campaigns.py']
    c = json.load(open(f)); bad = [] if c['agrees'] else ['the cached tables of the model did not agree with the code when last recomputed']
    for name, t in c['tables'].items():
        if sha256(f'{ROOT}/data/{name}.tsv') != t['sha256']: bad.append(f'data/{name}.tsv has changed since it was checked against the code')
    for p, h in c['code'].items():
        if sha256(f'{ROOT}/src/{p}') != h: bad.append(f'src/{p} has changed since the cached tables were checked against it')
    return bad


def check():
    f = f'{OUT}/MANIFEST.json'
    if not os.path.exists(f): print('no results/MANIFEST.json: run bash run_all.sh'); return 1
    old, new = json.load(open(f)), build(); n_fig, bad = figures(); bad += campaigns()
    for kind in ('inputs', 'results'):
        for name in sorted(set(old[kind]) | set(new[kind])):
            if old[kind].get(name) != new[kind].get(name): bad.append(f"{kind[:-1]} {'changed' if name in old[kind] and name in new[kind] else 'added' if name in new[kind] else 'missing'}: {name}")
    for p, v in new['packages'].items():
        if old['packages'].get(p) != v: bad.append(f"package {p}: {old['packages'].get(p)} when the results were written, {v} now")
    print('\n'.join(bad) if bad else f"results match the manifest ({len(new['results'])} results from {len(new['inputs'])} inputs; {n_fig} figures are the ones their inputs give; "
                                     f"the cached tables of the model are the ones the code writes)"); return 1 if bad else 0


def compare(path, strict=False):
    """This tree and its results against the manifest of another run, for instance the clean rebuild recorded with the
    release (docs/MANIFEST-release.json). If the inputs differ, the two runs are not runs of the same code and data:
    this is said, and with --strict the exit status is 1. Results are then counted: deterministic ones are identical;
    a sampled one (bayes_*, spec_curve, ipd_ecological and what is derived from them) can differ in its last digits on
    another machine or another version of the numerical libraries, which tests/test_article.py bounds."""
    other, new = json.load(open(path)), build(); bad = [n for n in sorted(set(other['inputs']) | set(new['inputs'])) if other['inputs'].get(n) != new['inputs'].get(n)]
    if bad:
        print(f'this tree is not the one of {path}: {len(bad)} inputs differ (' + ', '.join(bad[:12]) + (' ...' if len(bad) > 12 else '') + ')'); return 1 if strict else 0
    names = sorted(set(other['results']) & set(new['results'])); same = [n for n in names if other['results'][n] == new['results'][n]]; diff = [n for n in names if n not in same]
    missing = sorted(set(other['results']) - set(new['results']))
    print(f"inputs identical ({len(new['inputs'])}); results identical: {len(same)} of {len(names)}" + (f"; different: {', '.join(diff)}" if diff else '') + (f"; not computed here: {', '.join(missing)}" if missing else ''))
    pk = {p: (other['packages'].get(p), v) for p, v in new['packages'].items() if other['packages'].get(p) != v}
    if pk or other.get('python') != new['python']: print('environment differs: ' + ', '.join([f"python {other.get('python')} -> {new['python']}"] * (other.get('python') != new['python']) + [f'{p} {a} -> {b}' for p, (a, b) in pk.items()]))
    return 0


if __name__ == '__main__':
    if sys.argv[1:] == ['check']: sys.exit(check())
    if len(sys.argv) >= 3 and sys.argv[1] == 'compare': sys.exit(compare(sys.argv[2], strict='--strict' in sys.argv))
    os.makedirs(OUT, exist_ok=True); m = build(); json.dump(m, open(f'{OUT}/MANIFEST.json', 'w'), indent=1); print(f"wrote results/MANIFEST.json ({len(m['results'])} results, {len(m['inputs'])} inputs)")
