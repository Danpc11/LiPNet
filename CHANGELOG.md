# Changelog

Changes to the code and data. Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.11.5] – 2026-11-22
- Calculator: a value is placed against the thresholds in clinical use and not against multiples of the donor: gradient 10 mmHg, portal pressure 15 and 20 mmHg, flow 250 mL/min/100 g (`viscosity.PVP_THRESHOLDS` added to `FLOW_THRESHOLD` and `GRADIENT_THRESHOLD`). A gradient below the donor's is no longer labelled as steal unless the flow is also below the donor's.
- Calculator: haemoglobin accepted in place of the haematocrit (`viscosity.HB_TO_HCT`); the effect of each manoeuvre is read from `data/interventions.tsv` and shown as a reported effect, not a target; the page names the release it was built from; fields no longer overflow their card.
- `tests/test_indices.py`: the labels, the thresholds against the cut-offs of the series, the conversion and the release shown.
- No change to the model or to the data; results unchanged.

## [0.11.4] – 2026-11-22
- `.zenodo.json`: description and keywords of the current content (viscosity, distension, hepatectomies); third grant in the funding note. `CITATION.cff`: concept DOI and release date.
- README: the overview image is the one of the article as it stands (resistance load as architecture, blood and distension); the environment is created with `python3`.
- `requirements.txt` lists pytest, which the README runs; `requirements-lock.txt` pointed to a document that no longer exists.
- `.gitignore` names the result files that are kept.
- No change to the code that computes or to the data; results unchanged.

## [0.11.3] – 2026-11-22
- `data/hepatectomy.tsv`: columns `animal_group` and `baseline_group`. `hepatectomy.group_of` reads the first and `baseline_of` the second; groups are no longer inferred from study, flow before and n. The three states of Golriz 2022 are three groups of animals with one pooled state before. Results unchanged.
- Figure 3E prints the exact P of the trend with one state per group of animals.
- README: where the result files of a release are found.

## [0.11.2] – 2026-11-22
- Fixed: `ipd.ecological.reported_sigma` gave every group of a series the mean dispersion of that series; each group now takes its own (`gradient_sd` of the table). `results/ipd_ecological.json: reported_sigma_fit` changes; no other result does.
- One definition instead of several: the networks of the rigid range (`liver.B_CANOPY`, `human.B_CALIBRATED`; `viscosity.B_HUMAN` removed), the regions of a group (`distension.regions`, used by `bracket` and by `viscosity`), the rigid network of a remnant (`hepatectomy.rigid_each`), the donor references (`indices.NORMAL_GRADIENT`, `DEFAULT_DONOR_REF`, `DEFAULT_CVP`, read by `human`, `ipd`, `audit`, `predictions`, `distension` and `build_app`), `audit.sd_from_range` (the copy in `ipd.patient_level` removed).
- Removed, called by nothing: `network.Loaded.support_ok`, `src_component`, `cycle_rank`, `d_norm`, `cost` and `LN3`; `fem.darcy.calibrate`, `fem.mesh.elem_to_node`, `integrate`, `elem_mean`, `fem.network_link.share_vs_cost_exponent`; the `bed` branch of `paper_figures.plane`; unused imports, variables and arguments (`pooled_sigma(d)`, `reconstruct(tol=)`, `figstyle.bare(keep=)`).
- `build_app.py` no longer reads `results/within_study_fit.json`.
- Tests: `WITHIN_GROUP_SD` against the table and the patient-level tables; the calculator tests are skipped, not passed, without node.
- README, module docstrings and comments brought in line with the code.

## [0.11.1] – 2026-11-21
- `test_article.py`: builders for `results/ipd_patient_level.json` and `results/model_liver.json`, which the claims added in 0.11.0 read and a clone without results does not have; a test that every result the claims read is kept in the repository or has a builder. No change to code, data or results.

## [0.11.0] – 2026-11-21
- `hepatectomy.trend`: no P value for the seven states (five groups of animals); exact permutation P of Spearman's rho for every choice of one state per group (`exact_p`), and rho over the draws.
- `hepatectomy.state(corr=, base=)`, `run`, `covariance_sensitivity`: correlation between the group means before and after (0, 0.5, 0.8); the state before is shared by the states of the same animals (`group_of`, `baseline`).
- `hepatectomy.errors`, `configurations`: error of every admissible configuration and distance to the range. Removed `rms_log_error_of_the_midrange`, `observed_over_midrange` and the numeric `trend.p_value`; states gain `observed_over_model`.
- `viscosity.vessel_haematocrit`, `factor(micro=)`, `microvascular`: haematocrit of the vessels below 70 µm at 0.75, 0.5 and 0.25 of the systemic one. `clinical.rms_distance_to_the_rigid_range` replaces `rmse_rigid_midrange`.
- Figure 3E, F: rho for the seven states and for one state per group; observed over the range of each model.
- Tests for the four analyses. `test_article.py`: every claim labelled by where the article prints it (main text or supplement); claims added for the numbers that had none (two-sided sign test, network alone, whole graft, reported spreads, individual recipients, span of the intercepts, counts); `up` rounds a P value half up.

## [0.10.0] – 2026-11-21
- `viscosity.thresholds`: the clinical thresholds of flow and gradient on the scale of the loads, and the gradient at the flow threshold by haematocrit; unrounded regions for the six groups (`viscosity._regions`).
- `hepatectomy.band`, `animals_resected`.
- Figures: 2 (terms of the resistance load, thresholds, haematocrit), 3 (perfused livers and hepatectomies), 4 (outcome, risk inside a cohort, whole graft); Figure S8 of 0.9.0 is now Figure 3.
- Calculator: optional haematocrit (the ranges with the blood of the recipient) and the resistance load at which the two clinical thresholds agree; tested against `model/viscosity.py` with node.
- `test_article.py`: the numbers of the reframed article.

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
