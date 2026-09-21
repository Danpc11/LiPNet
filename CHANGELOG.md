# Changelog

All notable changes to this repository are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow the GitHub releases, each archived at Zenodo under the concept DOI [10.5281/zenodo.22861276](https://doi.org/10.5281/zenodo.22861276).

## [Unreleased]

### Fixed
- **The calculator mixed two models.** The intercepts came from the stratified fit (slope 1.96) while the predictions used the centred slope (1.66), so the curve did not reproduce the series it was anchored on. Slope, interval and intercepts now come from the same stratified model, and a user-supplied reference enters as logit p(z) = logit p_ref + beta·(z − z_ref), so the curve and its band pass through the anchor: at the reference gradient the band is exactly the rate given. Each fitted stratum now reproduces its own observed events (Wang 39, Ogura 20, Osman 7, Uemura 30, Yagi 8).
- **The centred analysis is no longer presented as equivalent to the stratified fit.** It works on transformed proportions and its interval ignores the dependence induced by centring within a study, so it is documented and used as a descriptive view on a common scale; the effect is taken from the stratified fit, with the meta-analysis as sensitivity.
- **The per-patient extrapolation is gone.** `attenuation()` no longer reports an "implied individual slope" from an assumed within-group SD; recovering an individual-level relation from group means is not possible here, and the claim that aggregation understates the effect has been removed from the README.
- Provenance labels now match the stored values: rows that carry only a gradient have `not applicable` in `PVP_basis` and `CVP_basis`, and the Uemura and Botha notes no longer mention a derived CVP that is not stored. A test asserts that no basis claims a value that is absent.
- `.zenodo.json` updated to 31 groups from 17 series; the calculator footer lists all 17; `0.2..2` in this file; `risk_curve` removed from the `plots.py` usage line; `recalibrate()` documented as taking individual outcomes only.


## [0.3.0] – 2026-09-21

### Added
- `indices.baseline_from_rate()` and `indices.recalibrate()`: calibration in the large, the standard first step of prediction-model updating (Steyerberg; Vergouwe et al., Stat Med 2017; Janssen et al., Can J Anesth 2009). The first gives the intercept from a centre's overall outcome rate and average gradient; the second fits only the intercept on a cohort with the slope held as an offset. README gains a section on how to obtain a baseline.
- The graphical abstract, lost from the README in an earlier rewrite, is back.
- `indices.meta_slope()`: random-effects (DerSimonian-Laird) meta-analysis of the per-stratum slopes, with tau2, Q and I2. Result: OR 6.69 per unit of zP (95% CI 3.12-14.4), tau2 = 0, Q = 1.47 on 4 df (p = 0.83), I2 = 0%, i.e. the five series estimate one common effect. It agrees with the stratified fit and does not spend degrees of freedom on the intercepts.
- `indices.overdispersion()`: Pearson chi2/df of the stratified fit (0.31 on 5 df), with the quasi-binomial interval; there is no extra-binomial spread.
- `indices.attenuation()`: regression dilution. The sampling error of the group means costs 3% of the slope; the within-group spread of the gradient implies a per-patient slope of about 3.3 (lambda 0.59), so aggregation understates rather than inflates the effect.
- `indices.centred()`: the within (fixed-effects) transformation. Subtracting each cohort's own mean log-odds and mean zP cancels the stratum intercept exactly and puts all 11 groups on one scale; they fall on a single line (weighted R2 = 0.91, slope 1.66 +- 0.32, OR 5.25 per unit of zP, 1.39 per mmHg), with plot `centred`.
- The calculator now uses the intercept-free centred slope and offers a **cohort-average anchor**: enter your cohort's overall outcome rate and its average portocaval gradient, and it places the patient relative to that average, which is the quantity the centring defines and the pair of numbers a centre can actually state.
- The calculator's plot is one curve now, not five: with a reference chosen it shows the absolute risk implied by that anchor with its 95% band and two markers (current gradient and gradient after the change), with the y axis scaled to the range in use; without an anchor it shows the published groups centred on their own cohort, which is the relation the slope comes from. The five-curve version is gone.
- `within_study` simplified: curves labelled at their right end, no legend box and no statistics box (those live in `forest` and `centred`); `risk_change` now uses the same centred slope as the calculator, so the two no longer disagree.
- Plot `forest`: per-stratum slope with its interval and the pooled estimate.
- `results/meta_slope.json`; three more rows in `results/sensitivity.tsv`; tests for all three.

## [0.2.2] – 2026-09-21

Code and methodology audit: three corrections that change the estimated effect.

### Fixed
- **A cohort was counted twice.** The Wang 2014 pressure groups (PVP at closure ≥20 / <20, 292 patients) are a second partition of the same cohort as its splenectomy groups (276 patients); both entered the model, inflating that study to 568 patients. The pressure partition is now excluded from the main analysis, and a test asserts that the Wang cohort contributes 276 patients once.
- **Intercepts are now per study AND outcome definition.** A series reporting two outcomes (primary graft dysfunction and graft loss) no longer shares a single intercept.
- **The denominator of Ishizaki 2012** is the outcome denominator (42, no SFSS), not the 31 recipients in whom pressure was measured.
- Bootstrap replicates that fail to converge are retried with Nelder-Mead and discarded if they still fail; the count is reported (`bootstrap_failures`, currently 0 of 2000). A profile-likelihood interval is computed as an independent check on the bootstrap.

### Added
- The calculator lets you choose the **reference level**: your own rate at the starting gradient, or any of the fitted strata (study x outcome), with its number of recipients, events and the zP range its groups actually span. Choosing a series shows the risk it would go from and to, warns when the entered gradient is outside that range, and highlights its curve.

### Changed
- The long disclaimer moved from the top of the calculator to a short note under the risk curve, and the field hints were trimmed.
- README cut roughly in half: results first, caveats collected at the end instead of scattered through the text.
- Twelve tests folded into five that check the things worth checking (indices and baseline, data and provenance, the fitted effect, plots and calculator, README against results).
- With these corrections the within-stratum effect is larger and less precise: **OR 7.09 per unit of zP (95% CI 3.37–20.5), 1.48 per mmHg** (slope 1.96, bootstrap 1.22–3.02, profile likelihood 1.18–2.83), on 11 groups, 104 events in 734 recipients. It holds by outcome (1.87 for SFSS, 2.03 for mortality) and dropping any study (1.67–2.16).
- `data/series.tsv` gains `gradient_mmHg` and `gradient_basis`: the series that report the portocaval gradient itself (Uemura 2016, Ishizaki 2012, Botha 2010, Ogura 2010, Chan 2011, Yamada 2008) no longer carry invented PVP/CVP pairs, and `indices.zP()` accepts a gradient directly.
- The redundant `CVP_measured` column is removed; `CVP_basis` is the single provenance field.


## [0.2.1] – 2026-09-21

Five series with measured portocaval gradients were added; they refute the pooled common-intercept model and replace it with a within-study effect.

### Added
- `data/series.tsv` grows to 31 groups from 17 series: Uemura 2016 (Surgery 159:1623, three GRWR groups with the final PVP−CVP gradient and SFSS counts), Yao 2018 (Transplantation 102:623), Kanetkar 2017 (J Clin Exp Hepatol 7:235, PVP by direct portal cannulation), Ishizaki 2012 (Liver Transpl 18:305, left lobes without modulation, gradient 12.4 mmHg and no SFSS in 42) and Botha 2010 (Liver Transpl 16:649, hemiportocaval shunt, gradient 18 → 5 mmHg).
- `indices.within_study_fit()`: one intercept per study and a common slope, with a bootstrap interval for the slope and a leave-one-study-out check. Result: logit(p) = alpha_study + 1.34·zP, OR 3.81 per unit of zP (95% CI 2.44–6.37), OR 1.31 per mmHg of gradient; the slope stays between 1.13 and 2.03 when any study is dropped.
- `indices.risk_after(baseline_rate, delta_zP, beta)`: converts a change of gradient into an absolute risk given the user's own baseline rate.
- Plots `within_study` (observed groups and one fitted curve per study) and `risk_change` (odds ratio and absolute risk against the change of gradient, for several baseline rates).
- Tests for the fitted effect, the leave-one-study-out stability, the odds shift, and the absence of an absolute-risk claim in the calculator.

### Changed
- The calculator no longer reports a pooled absolute risk. It asks for the pressure after the planned manoeuvre and returns the odds ratio with its interval, and an absolute risk only when the user supplies their own outcome rate at the starting gradient. Its plot now shows one curve per study.
- Main analysis excludes the groups defined only by a cut-off (now Vasavada 2014, Yao 2018 and Kanetkar 2017).
- README rewritten around what the data support: the absolute level is centre-specific (0% to 17% at the same gradient), the effect of changing the gradient is transferable.

### Fixed
- The `threshold` provenance is now counted as imputed in the sensitivity filters.
- `within_study_fit(boot=0)` no longer fails (used by the leave-one-study-out loop).


## [0.2.0] – 2026-09-20

Methodological and calculator fixes after two rounds of code review.

### Changed
- `data/series.tsv` now holds raw values only, with provenance columns `PVF_basis`, `PVP_basis`, `CVP_basis`, `events_basis` (`reported`, `derived`, `imputed`, `threshold`, `group mean`, `not applicable`) and a column `in_main_analysis`; the stored `zP` and `zF` columns were removed and are recomputed by `src/indices.py`, the single source used by the fit, the plots, the tests and the calculator.
- Groups defined only by a flow cut-off (Vasavada 2014, PVF above and below 190 mL/min/100 g) are excluded from the main analysis, drawn hollow in `series_flow` and reported in sensitivity. Main flow result: 11 groups, Spearman ρ = 0.35, p = 0.29 (13 groups, ρ = 0.23, p = 0.46 when included).
- Sensitivity analysis is now per index: a group is excluded for `zP` only if PVP or CVP is imputed or filled, and for `zF` only if PVF or the donor reference is. `results/sensitivity.tsv` lists all analyses. The pooled logistic fit cannot be repeated without imputed CVP (only Yamada 2008 has fully reported pressures, with no events); this is stated in the README as the main limitation.
- Notes column translated to English and made explicit about every derivation (PVF divided by graft weight, PVP from gradient plus CVP, deaths from survival percentages with the unrounded value).
- Calculator: the normal portocaval gradient is fixed at 5 mmHg (the value the coefficients were fitted on) instead of user-editable; the curve grid extends to `zP` 0.5–8 and the range of the fitted groups (1.43–2.68) is stored so that values outside it are labelled as extrapolation with a band computed for the same `zP`.
- Calculator text ("4 of the 5 groups had no reported CVP") and coefficients are injected by `build_app.py` from the analysis outputs.
- README results updated to the current analysis (intercept −6.40; flow result with and without cut-off groups; per-index sensitivity).

### Added
- `indices.compute()` records every value it fills in (`CVP_filled`, `donor_ref_filled`), counts it as imputed regardless of the label, and raises an error when a label contradicts the data (e.g. a blank CVP marked `reported`).
- Calculator input validation: a blank CVP is taken as 5 mmHg and flagged, an explicit 0 is a measured zero; PVP and CVP must be non-negative; graft weight, GRWR, recipient weight and donor reference must be positive numbers (a blank donor reference falls back to 90 mL/min/100 g and says so); a portal flow of 0 gives `zF = 0`, is labelled "no portal inflow", and the ratio `zP/zF` and its interpretation are omitted.
- Separate labels for pressure (`tagP`: zero or negative gradient, steal risk, window, above window, high) and flow (`tagF`: no inflow, below donor, 1–2×, 2–3×, >3× donor).
- Tests: incomplete edits are caught, cut-off groups are excluded, README figures are checked against `results/`, calculator has no editable gradient and validates inputs. Tests run the full bootstrap so they never leave a short-bootstrap JSON that would alter the built calculator.
- `data/interventions.tsv` (pressure and flow change per inflow-modulation manoeuvre, with source) replaces values hard-coded in `plots.py`.
- `src/model/merge_campaigns.py` to merge per-size campaign outputs; `CITATION.cff` with version, date and concept DOI; `.zenodo.json`; `index.html` redirect for GitHub Pages; `CHANGELOG.md`.

### Fixed
- A blank CVP was read as 0 by the calculator (risk jumped from 13.1% to 49.6% with the default inputs).
- Changing the reference gradient in the calculator invalidated the fitted coefficients.
- Recomputed indices were not the ones used by the correlations and the fit.
- Uncertainty band outside the curve range was taken from the curve end-point instead of the actual `zP`.
- Negative or zero graft weight produced a negative `zF`; donor reference 0 was silently replaced by 90.
- A shared `tag()` labelled a zero pressure gradient as "no portal inflow".
- `plane` and `risk_curve` plots included cut-off groups.
- The token `n/a` in provenance columns was read by pandas as missing; replaced by `not applicable`.
- Unused functions removed from `src/model/network.py` (430 → 244 lines); leftover paper-figure scripts and the derived `data/graph3d_exponents.tsv` removed.

### Removed
- `data/graph3d_exponents.tsv` (exponents are computed on the fly), stored `zP`/`zF` columns, editable gradient field in the calculator.

## [0.1.0] – 2026-09-20

First public release. Archived as [10.5281/zenodo.22861277](https://doi.org/10.5281/zenodo.22861277).

### Added
- Extraction table of 12 published LDLT series (22 groups, 1,026 recipients) with the source table or page for every value.
- Indices `zP = (PVP − CVP)/5` and `zF = graft PVF per 100 g / donor PVF per 100 g`; Spearman correlations with outcome; pooled binomial logistic fit of SFSS on `zP` with 2000-resample bootstrap band.
- 14 stand-alone plots (`src/plots.py`): nomogram, risk curve, series by pressure and by flow, flow–pressure plane, resection, load curve, interventions, 2D and 3D optimal networks, scaling, allometry, shear profile, certificate.
- Perfusion-network model of the liver (`src/model/`): 2D and hemispheric 3D graphs, dissipation-optimal tree at fixed maintenance cost, closed-form scaling and allometric closure, campaign scripts with cached outputs.
- Browser calculator `sfss_calculator.html`, built by `src/build_app.py`.
- `run_all.sh`, tests with GitHub Actions, PolyForm Noncommercial 1.0.0 licence, patient-level table and cohort template.
