"""Cached fits are reused only while the table and the model code are the ones that produced them.

    results/bayes_inputs.json     for each cached fit, the fingerprint of the inputs it was computed from
    LIPNET_REBUILD=1              ignore every cache

The fingerprint is the SHA-256 of data/series.tsv, src/bayes.py and src/indices.py. A cache whose fingerprint is not
the current one is stale: predictions.py then runs bayes.py again instead of reading it.
"""
import hashlib, json, os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); OUT = f'{ROOT}/results'
INPUTS = ('data/series.tsv', 'src/bayes.py', 'src/indices.py')
STAMPS = f'{OUT}/bayes_inputs.json'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()


def fingerprint():
    h = hashlib.sha256()
    for f in INPUTS: h.update(f.encode()); h.update(sha256(f'{ROOT}/{f}').encode())
    return h.hexdigest()


def _stamps():
    return json.load(open(STAMPS)) if os.path.exists(STAMPS) else {}


def stamp(name):
    """Record that results/<name>.json was computed from the current inputs."""
    s = _stamps(); s[name] = fingerprint(); os.makedirs(OUT, exist_ok=True); json.dump(dict(sorted(s.items())), open(STAMPS, 'w'), indent=1)


def fresh(name):
    """True if results/<name>.json exists, was computed from the current inputs and no rebuild was asked for."""
    return (not os.environ.get('LIPNET_REBUILD')) and os.path.exists(f'{OUT}/{name}.json') and _stamps().get(name) == fingerprint()


def stale(names):
    return [n for n in names if not fresh(n)]
