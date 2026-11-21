# Changelog

Changes to the code and data. Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.9.0] – 2026-11-21
- `model/viscosity.py`: resistance of the network at the haematocrit of the recipient (in vivo viscosity law, vessel by vessel); `published_check` against perfused livers (`data/kamikado2004.tsv`); `clinical`.
- `model/hepatectomy.py`, `data/hepatectomy.tsv`: published hepatectomies against the rigid network and the three relations.
- `distension.bracket(viscosity=)`, `distension.joint(haematocrit=)`; without the argument the results are unchanged.
- Figure S8; `tests/test_external.py`.

## [0.8.0] – 2026-11-21
- Documentation reduced to `README.md`, `THEORY.md` and this file; `docs/` keeps only `MANIFEST-release.json`.
- Removed the documentation-content tests.
- Calculator: resistance load card, with the range of the rigid network (two inflows) and of the three passive relations for the flow load and CVP entered, computed in the page with the constants of `model/`; a reference series uses its own level with the pooled slope (`bayes` keys `alphas_pooled_slope`, `pred_band_pooled_slope`), curve and band from the same draws; credible intervals; the odds ratio is labelled as a within-cohort association; footer with the 18 series and the authors; stray character in the flow card removed. `src/model/sfss_calculator.html` (old copy) removed.

## [0.7.6] – 2026-11-20
- `manifest.py`: `inputs_fingerprint`; `stamp-release` and `verify-release`. CI verifies the release manifest on tags.
- `distension.MEASURED_UP_TO`; Figure 2A marks the highest measured flow.

## [0.7.5] – 2026-11-20
- `compartments.solve` returns convergence, iterations and residual and raises `NotConverged`; `solver_check`.
- Root bracketing and residual checks in `distension.distensible`, `_bisect` and `global_variant`.
- `distension.wang_contrast` compares distributions; Figure 2C.
- `bayes.censored_survival`.

## [0.7.4] – 2026-11-20
- Pressure-flow relation on the portocaval gradient (`GREENWAY`, `f` 0.30 to 0.35).
- `bracket` places point estimate and interval with one rule; `consistency`, `heterogeneity`, `wang_contrast`.

## [0.7.3] – 2026-11-20
- `model/compartments.py`: distensibility vessel by vessel.
- `distension.robustness`, `joint`; errors of the predictions in `published_check`.
- `predictions.p2_performance`; midpoint RMSE removed.

## [0.7.2] – 2026-11-20
- Third distensibility relation (linear pressure-flow); compliance doubling pressure 8 to 13 mmHg.

## [0.7.1] – 2026-11-20
- Cube law with $P_d$ the mean of upstream and downstream pressures; `distension.published_check`.

## [0.7.0] – 2026-11-20
- `model/distension.py`; Figure 4 shows model predictions.

## [0.6.2] – 2026-11-19
- `test_article.py` builds the results it needs; `test_cache.py` recomputes `predictions.json`; CI in three jobs; `docs/MANIFEST-release.json`.

## [0.6.1] – 2026-11-19
- `predictions.p2_uncertainty`, sign-test sensitivity; `model/verify_campaigns.py`; figure provenance; tolerances in `test_article.py`.

## [0.6.0] – 2026-11-19
- Column `configuration` in `data/series.tsv`; sign test over cohorts; `p2_model_against_observed`, `same_cohort_contrast`, `p2_reference_sweep`.
- `cache.py`, `manifest.py`, `run_all.sh --clean`, `requirements-lock.txt`, `test_article.py`.

## [0.5.2] – 2026-11-19
- Anatomical colours for the vessels.

## [0.5.1] – 2026-11-19
- `figstyle.py` and `paper_figures.py` replace `main_figures.py` and `human_figures.py`; `test_figures.py`.

## [0.5.0] – 2026-11-19
- Wu 2011 added to `data/series.tsv` (32 groups, 18 series).

## [0.4.5] – 2026-11-19
- Measured generations of the human cast; `human.two_inflows_measured`, `first_week`.

## [0.4.4] – 2026-11-18
- `model/human.py`: the network in physical units, growth and regrowth; branch-length factor; shared lobular drop.

## [0.4.3] – 2026-11-18
- `model/liver.py`, `audit.py`; fits independent of the number of cores.

## [0.4.2] – 2026-10-25
- `predictions.p2_reference_sweep`; `plots.whole_graft`. Fixed: Ogura 2010 dispersions stored as standard deviations; execute bit of `run_all.sh`.

## [0.4.1] – 2026-10-21
- Font-manager notices silenced.

## [0.4.0] – 2026-10-21
- `src/fem/` (lobule finite elements) and `src/ipd/` (patient-level layer).

## [0.3.1] – 2026-09-23
- Calculator: one model for slope and intercepts, validated case export, independent outcome selectors.

## [0.3.0] – 2026-09-22
- `bayes.py` (hierarchical model, NUTS); `predictions.py`; `model/whole_graft.py`.

## [0.2.3] – 2026-09-21
- `indices.baseline_from_rate`, `recalibrate`, `meta_slope`, `overdispersion`, `attenuation`, `centred`.

## [0.2.2] – 2026-09-21
- Second partition of Wang 2014 excluded; intercepts per study and outcome; column `gradient_mmHg`.

## [0.2.1] – 2026-09-21
- 31 groups from 17 series; `indices.within_study_fit`, `risk_after`.

## [0.2.0] – 2026-09-20
- Raw values with provenance columns in `data/series.tsv`; `data/interventions.tsv`; `CITATION.cff`, `.zenodo.json`.

## [0.1.0] – 2026-09-20
- Extraction table, indices, network model, plots, calculator, tests.
