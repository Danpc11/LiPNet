# LiPNet - Liver Perfusion Network model

### A network model of hepatic perfusion and haemodynamic load in partial liver grafts

![Python](https://img.shields.io/badge/Python-3.11%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/LiPNet/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/LiPNet/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/LiPNet/)
![Version](https://img.shields.io/badge/version-0.11.5-1f6feb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22861276.svg)](https://doi.org/10.5281/zenodo.22861276)

Code and data for the article *Vascular resistance decouples portal flow from pressure in small-for-size liver grafts*.

LiPNet represents the liver as a dissipation-optimal perfusion network and expresses the burden of a graft as three
loads normalised to the healthy donor, related by

$$z_P = z_F\, z_R \qquad \text{pressure load} = \text{flow load} \times \text{resistance load}.$$

The model, its equations and its assumptions are in **[THEORY.md](THEORY.md)**.

<p align="center">
  <img src="assets/graphical_abstract.png" alt="Portal pressure in a partial graft is portal flow times the resistance of the graft, and that resistance has three terms: architecture, blood and distension" width="760">
</p>

## Installation

Python 3.11 or later.

```bash
python3 -m venv .venv && . .venv/bin/activate    # inside the environment the interpreter is `python`
pip install -r requirements-lock.txt     # exact environment of the release (reproduces every printed digit)
# or: pip install -r requirements.txt   # recent versions; sampled numbers may move in the second decimal
```

## Running

```bash
python -m pytest tests -q                # about 3 minutes on a fresh clone
bash run_all.sh                          # every result, figure and the calculator; reuses valid cached fits (~5 min)
bash run_all.sh --clean                  # empties results/ and recomputes everything, fits included (~30 min)
python src/manifest.py check             # results, figures and cached tables are the ones their inputs give
```

Single steps:

| Command | Output |
|---|---|
| `python src/indices.py` | normalised loads, stratified fits, `results/sensitivity.tsv` |
| `python src/bayes.py` | hierarchical outcome model and its sensitivity analyses, `results/bayes_*.json`, `spec_curve.json` |
| `python src/predictions.py` | predictions P1 to P4, `results/predictions.json` |
| `python src/model/liver.py`, `human.py` | the network in a human liver and in a graft, `results/model_liver.json`, `model_human.json` |
| `python src/model/distension.py` | passive distensibility relations, `results/model_distension.json` |
| `python src/model/viscosity.py` | resistance of the network with the blood of the recipient, `results/model_viscosity.json` |
| `python src/model/hepatectomy.py` | published hepatectomies against the model, `results/model_hepatectomy.json` |
| `python src/fem/run_fem.py all` | lobule finite elements, `results/fem_*` |
| `python src/ipd/run_ipd.py all` | patient-level tables, within-group dispersion model, `results/ipd_*.json` |
| `python src/paper_figures.py [1 2 S3 ...]` | the figures of the article, `results/Figure_*.pdf/.png` |
| `python src/plots.py all` | stand-alone working plots |
| `python src/build_app.py` | `sfss_calculator.html` |

The tables of the network campaigns in `data/` (`graph_scaling_all.tsv`, `graph3d_all.tsv`, `certificate_small.tsv`)
are cached. `python src/model/verify_campaigns.py` recomputes them and compares every number; the commands that write
them are at the end of `run_all.sh`.

## Methodology

- **Network model** (`src/model/`). Lobules between a portal and a hepatic venous tree; the tree that minimises
  dissipation at fixed maintenance cost $\sum_e m_e^b$, with $b\in[0.5,1]$ reported rather than fitted. Closed-form
  canopy, explicit 2D and 3D graphs, exhaustive check on seven lobules. Calibrated in physical units against a human
  portal cast.
- **Partial grafts.** A graft is a sub-network; the hepatic artery joins at the sinusoid. The network predicts $z_R$ of
  a graft from its flow load alone (`predictions.p2_model_against_observed`).
- **Passive distensibility** (`model/distension.py`, `model/compartments.py`). Three published relations of the
  normal liver applied to the whole liver or vessel by vessel, as scenarios; uncertainty drawn jointly.
- **Blood viscosity** (`model/viscosity.py`). Resistance of the same network at the haematocrit of the recipient,
  from a published in vivo viscosity law applied vessel by vessel; checked against perfused livers
  (`data/kamikado2004.tsv`).
- **External check** (`model/hepatectomy.py`, `data/hepatectomy.tsv`). Published hepatectomies, the same liver
  before and after resection, against the rigid network and the three relations.
- **Clinical data** (`data/series.tsv`). Groups of published adult living-donor series, one row per reported group,
  with a provenance label on every value. Groups defined only by a cut-off and second partitions of a counted cohort
  are excluded from the main analysis.
- **Outcome model** (`src/bayes.py`). Hierarchical binomial regression of outcome on $z_P$ centred within study, NUTS;
  prior sensitivity, specification curve, internal-external validation, recovery simulation.
- **Lobule** (`src/fem/`). Finite elements of one hexagonal lobule: anisotropic Darcy flow, transport with
  haemoglobin binding, sinusoid-hepatocyte exchange, nitric oxide wall model.
- **Patient level** (`src/ipd/`). The series that print their recipients one by one; group events integrated over
  the within-group spread of $z_P$. Kaplan-Meier curve reconstruction is implemented and tested on simulated cohorts,
  and not applied to any series.

## Repository layout

```text
src/
  indices.py            loads, provenance checks, stratified fits, recalibration
  bayes.py              hierarchical outcome model
  predictions.py        P1-P4
  audit.py              sampling intervals and checks of the extraction table
  model/                network model (network, network3d, scaling, liver, human, whole_graft,
                        distension, compartments, viscosity, hepatectomy, certify, graph_scaling,
                        merge_campaigns, verify_campaigns)
  fem/                  lobule finite elements
  ipd/                  patient-level layer
  paper_figures.py      figures of the article; figstyle.py, their style
  plots.py              working plots
  cache.py, manifest.py fingerprints of cached fits; manifest of a run
  build_app.py          browser calculator
data/                   series.tsv, interventions.tsv, patients.tsv, patient_level/, hepatectomy.tsv,
                        kamikado2004.tsv, cohort_template.tsv, curves/ (template), cached network tables
results/                generated; the cached fits and predictions.json are kept
tests/                  pytest suite
docs/MANIFEST-release.json   record of the clean rebuild of the release
run_all.sh              whole pipeline
```

## Data

`data/series.tsv` holds published values only; normalised loads are computed by `src/indices.py`. Each value has a
`*_basis` column (`reported`, `derived`, `imputed`, `threshold`, `not applicable`) and the `notes` column gives the
table and page it comes from. Structural columns (`centre_id`, `cohort_id`, `recruitment_start/end`, `overlap_set`,
`configuration`, `outcome_definition`, `outcome_horizon`) separate publications, cohorts and centres.

| Series | Journal | DOI |
|---|---|---|
| Troisi 2003 | Liver Transpl 2003;9:S36–41 | [10.1053/jlts.2003.50200](https://doi.org/10.1053/jlts.2003.50200) |
| Troisi 2005 | Am J Transplant 2005;5:1397–1404 | [10.1111/j.1600-6143.2005.00850.x](https://doi.org/10.1111/j.1600-6143.2005.00850.x) |
| Yagi 2005 | Liver Transpl 2005;11:68–75 | [10.1002/lt.20317](https://doi.org/10.1002/lt.20317) |
| Yagi 2006 | Transplantation 2006;81:373–378 | [10.1097/01.tp.0000198122.15235.a7](https://doi.org/10.1097/01.tp.0000198122.15235.a7) |
| Yamada 2008 | Am J Transplant 2008;8:847–853 | [10.1111/j.1600-6143.2007.02144.x](https://doi.org/10.1111/j.1600-6143.2007.02144.x) |
| Botha 2010 | Liver Transpl 2010;16:649–657 | [10.1002/lt.22043](https://doi.org/10.1002/lt.22043) |
| Ogura 2010 | Liver Transpl 2010;16:718–728 | [10.1002/lt.22059](https://doi.org/10.1002/lt.22059) |
| Ou 2010 | Transplant Proc 2010;42:876–878 | [10.1016/j.transproceed.2010.02.064](https://doi.org/10.1016/j.transproceed.2010.02.064) |
| Chan 2011 | Liver Transpl 2011;17:115–121 | [10.1002/lt.22034](https://doi.org/10.1002/lt.22034) |
| Wu 2011 | Liver Transpl 2011;17:1035–1045 | [10.1002/lt.22326](https://doi.org/10.1002/lt.22326) |
| Ishizaki 2012 | Liver Transpl 2012;18:305–314 | [10.1002/lt.22440](https://doi.org/10.1002/lt.22440) |
| Vasavada 2014 | Int J Surg 2014;12:177–180 | [10.1016/j.ijsu.2013.12.006](https://doi.org/10.1016/j.ijsu.2013.12.006) |
| Wang 2014 | Surg Today 2015;45:979–985 | [10.1007/s00595-014-0999-9](https://doi.org/10.1007/s00595-014-0999-9) |
| Alim 2016 | Liver Transpl 2016;22:1643–1648 | [10.1002/lt.24523](https://doi.org/10.1002/lt.24523) |
| Uemura 2016 | Surgery 2016;159:1623–1630 | [10.1016/j.surg.2016.01.009](https://doi.org/10.1016/j.surg.2016.01.009) |
| Kanetkar 2017 | J Clin Exp Hepatol 2017;7:235–246 | [10.1016/j.jceh.2017.01.114](https://doi.org/10.1016/j.jceh.2017.01.114) |
| Osman 2017 | Hepatol Res 2017;47:293–302 | [10.1111/hepr.12727](https://doi.org/10.1111/hepr.12727) |
| Yao 2018 | Transplantation 2018;102:623–631 | [10.1097/TP.0000000000002047](https://doi.org/10.1097/TP.0000000000002047) |

`data/hepatectomy.tsv` holds the published hepatectomies (one row per state, before and after resection, with the
animals it was measured in, the measurement before it uses and the place in the paper of each value) and `data/kamikado2004.tsv` the resistances of perfused livers by haematocrit,
digitised from the figure of the paper.

To add a series, append rows to `data/series.tsv` with every `*_basis` filled; `tests/test_indices.py` rejects a
value without provenance. `data/cohort_template.tsv` is the format for a centre's own recipients.

## Tests and reproducibility

| Test file | Checks |
|---|---|
| `test_indices.py` | loads, provenance, fits, calculator |
| `test_liver.py`, `test_human.py` | network in a human liver, distensibility relations and solvers |
| `test_external.py` | viscosity law, perfused livers, published hepatectomies |
| `test_campaigns.py`, `test_fem.py`, `test_ipd.py` | network campaigns, lobule model, patient-level layer |
| `test_cache.py` | cached fits and `predictions.json` are those of the current code; release manifest |
| `test_figures.py` | figures and their provenance |
| `test_article.py` | each number quoted in the article against `results/` |

Sampled numbers are reproducible to the last digit on `requirements-lock.txt`; elsewhere `test_article.py` allows a
stated Monte Carlo tolerance (`LIPNET_EXACT=1` removes it).

Only the cached fits and `predictions.json` are kept in `results/`. The other result files (`model_*.json`,
`ipd_*.json`, figures) are written by `run_all.sh`; those of the clean rebuild of a release are attached to its GitHub
release (`LiPNet-<version>-results.zip`) and their hashes are in `docs/MANIFEST-release.json`.

A cached fit is used only if its fingerprint matches `data/series.tsv`, `src/bayes.py` and `src/indices.py`.
`run_all.sh` ends by writing `results/MANIFEST.json` (hashes of inputs and results, package versions).
`docs/MANIFEST-release.json` is the manifest of the clean rebuild of the release; to compare a run with it:

```bash
python src/manifest.py compare docs/MANIFEST-release.json
python src/manifest.py verify-release    # on a full clone: the commit named in the manifest holds the rebuilt code
```

Making a release: rebuild with `bash run_all.sh --clean`, copy `results/MANIFEST.json` to
`docs/MANIFEST-release.json`, commit the code, run `python src/manifest.py stamp-release`, commit the manifest, tag.

## Calculator

[danpc11.github.io/LiPNet](https://danpc11.github.io/LiPNet/) (`sfss_calculator.html`, built by `src/build_app.py`
from `src/app_template.html`) computes $z_F$, $z_P$ and $z_R$; places $z_R$ against the rigid network and the passive
relations at the flow load entered, with normal blood and at the haematocrit or haemoglobin entered; labels each value
against the thresholds in clinical use (gradient 10 mmHg, portal pressure 15 and 20 mmHg, flow 250 mL/min/100 g); and
gives the odds ratio for a change of gradient from the primary fit of `bayes.py`, converted to a risk only with a
baseline rate the user supplies or the level of a published series.
Everything runs in the browser. Research tool, not a medical device.

## Citation

Cite the article and the archived release (`CITATION.cff`). Concept DOI for all versions:
[10.5281/zenodo.22861276](https://doi.org/10.5281/zenodo.22861276).

## Licence

Code: [PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0). Clinical values in
`data/` belong to the authors of the original publications.
