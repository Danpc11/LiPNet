# Reproduction record

What was run to produce the numbers of the article, what is checked automatically, and what a third party can run to
confirm it. The record of the run is [MANIFEST-release.json](MANIFEST-release.json).

## The run of the release

| | |
|---|---|
| Command | `bash run_all.sh --clean`, on a copy of the repository with an empty `results/` |
| Environment | Python 3.13.18 on Linux x86-64, a new virtual environment installed from `requirements-lock.txt` (26 packages at exact versions) |
| What it does | Empties `results/`, samples every hierarchical fit again, recomputes the model campaigns behind the cached tables in `data/`, rebuilds every result, figure and the calculator, writes the manifest |
| Then | `python -m pytest tests` with `LIPNET_EXACT=1`, and `python src/manifest.py check` |

`MANIFEST-release.json` holds the SHA-256 of every input (`data/`, `src/`, `run_all.sh`, the two requirements files) and
of every numeric result of that run, with the versions of Python and of the packages.

## What is compared with what

| Object | Kept in the repository | How it is tied to the code |
|---|---|---|
| Hierarchical fits (`results/bayes_*.json`, `spec_curve.json`) | yes, they take minutes | Each carries the fingerprint of `data/series.tsv`, `src/bayes.py` and `src/indices.py`; a fit with another fingerprint is not read (`tests/test_cache.py`) |
| `results/predictions.json` | yes | `tests/test_cache.py` recomputes everything in it that is not sampled and compares; the sampled part must equal the cached fits |
| Tables of the model (`data/graph_scaling_all.tsv`, `graph3d_all.tsv`, `certificate_small.tsv`) | yes | `src/model/verify_campaigns.py` reruns the campaigns and compares every number (relative tolerance $10^{-9}$); `run_all.sh --clean` runs every size, the tests the two smallest |
| Every other result | no | Rebuilt by `run_all.sh`; `tests/test_article.py` builds the two it needs if they are missing |
| Figures | no | `results/figures_provenance.json`: hash of what each figure reads, of the scripts and of the images, checked by `python src/manifest.py check` |
| Numbers quoted in the article | in `tests/test_article.py` | Counts and arithmetic exactly; sampled numbers within 0.05 (posterior mean), 0.10 (limit of a 95% interval) or 0.30 (interval of a new centre, log score) |

## Continuous integration

`.github/workflows/tests.yml` has three jobs.

| Job | When | Environment | What it runs |
|---|---|---|---|
| `locked` | every push | Python 3.13, `requirements-lock.txt` | `run_all.sh`, the tests, `manifest.py check`, and the comparison with the manifest of the release |
| `latest` | every push | Python 3.11, `requirements.txt` at whatever versions it resolves to | the tests, on a clone where nothing has been run |
| `clean` | a tag `v*`, or on request (Actions, Run workflow, clean) | Python 3.13, `requirements-lock.txt` | `run_all.sh --clean`, the tests, both checks; the manifest and the results are kept as an artifact |

The `clean` job is the independent run: a machine that is not the authors', starting from nothing. It is meant to be
run on the tag of the release.

## What identical means

On one machine and one environment the pipeline is deterministic: two clean rebuilds give the same bytes for every
numeric result and every figure. Across machines the deterministic results (counts, the network model, the tables of
the campaigns) are still identical. A sampled result can differ in its last digits, because the sampler amplifies the
last bit of a floating-point operation, and a figure can differ in its bytes if the machine has another font.
`python src/manifest.py compare docs/MANIFEST-release.json` reports how many results are identical and names the
others; `tests/test_article.py` bounds how far a sampled number may move.

## For an audit

Audit a commit, not the branch. Files read from the `main` branch through the web can be served from a cache for
some minutes after a push, so two files of one reading can belong to two commits. Clone the tag, or read the files at
the permanent address of the commit:

```bash
git clone --branch v0.7.5 --depth 1 https://github.com/Danpc11/LiPNet && cd LiPNet
python -m venv .venv && . .venv/bin/activate && pip install -r requirements-lock.txt
python -m pytest tests -q                                  # on the clone as it is
bash run_all.sh --clean                                    # everything from nothing, about 30 minutes
python -m pytest tests -q && python src/manifest.py check && python src/manifest.py compare docs/MANIFEST-release.json
```
