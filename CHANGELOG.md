# Changelog

All notable changes to this repository are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow the GitHub releases, each archived at Zenodo under the concept DOI [10.5281/zenodo.22861276](https://doi.org/10.5281/zenodo.22861276).

## [0.4.2] – 2026-10-25

### Added
- `predictions.p2_reference_sweep`: the two references behind $z_R$ are swept rather than assumed. The default donor flow (80 to 100 mL/min/100 g, in steps of 5; Troisi and Chan keep the donor value they report) is crossed with the normal gradient (4, 5, 6 mmHg), and the assumed CVP (3 to 9 mmHg) is swept separately at the default references. All five paired groups stay below one in 14 of the 15 flow-by-gradient combinations; the exception is the extreme corner, a reference of 100 with a normal gradient of 4 mmHg, where Ishizaki 2012 reaches 1.03. Crossing all three at their extremes (100, 4 mmHg, CVP 3) lifts Wang 2014 without splenectomy to 1.14, which is reported as the worst corner. The result sits in `results/predictions.json` under `P2_discordance.reference_sweep` and is asserted by `tests/test_indices.py`.
- `plots.whole_graft`: the whole-graft extension as one figure (resistance chain, loads along a narrowing of each anastomosis, and the two trajectories in the flow-pressure plane), which the manuscript cites as Fig. S2.

### Fixed
- `data/series.tsv` stored Ogura 2010's printed dispersions (0.3 and 0.4 mmHg) as standard deviations, so the measurement Monte Carlo in `bayes.py`, which divides the stored spread by sqrt(n), understated that series' uncertainty about ninefold. They are standard errors (see 0.4.0) and are now stored as the converted standard deviations, 2.78 and 2.62 mmHg, with the reason in the row notes. The same rows now note that Table 3 of the paper lists 4 and 12 deaths over the whole follow-up, and why the 2-year counts derived from the survival curve are the ones used.
- `run_all.sh` was committed without its execute bit, so `./run_all.sh` failed with "Permission denied".
- `bayes.py` asked XLA for four host devices whatever the machine had; on a two-core machine the thirty refits of the recovery simulation exhausted memory during compilation. It now asks for `min(4, cpu_count)`.
- `docs/IPD.md` quoted a pooled within-group spread of 0.76 $z_P$ units and a reported-spread fit of +1.70 (+0.34 to +3.40); the code gives 0.69 (3.45 mmHg) and +1.71 (+0.34 to +3.28).
- `results/bayes_recovery.json` is regenerated with the ten replicates the script asks for (the cached file had eight).

## [0.4.1] – 2026-10-21

### Fixed
- The figure scripts silence the matplotlib font-manager notices. `Liberation Sans` heads the font list and is absent on macOS, so every lookup logged a notice while the Arial fallback drew the figure correctly; `warnings.filterwarnings` never reached them because they come from the logger rather than the warnings system. The font list is unchanged, so figures built on Linux keep the typeface of the published ones.

## [0.4.0] – 2026-10-21

### Added
- `src/fem/`: a finite-element model of one hepatic lobule, which opens up the single graph edge the network model uses for a lobule. Anisotropic Darcy perfusion of the hexagonal lobule reproduces a perfusion of 1.05 mL/min/g, a sinusoidal transit time of 8.6 s and a wall shear stress of 0.19 Pa; transport carries oxygen with its haemoglobin-bound pool (periportal 56 mmHg to pericentral 33, 33% extraction), ethanol through ADH and CYP2E1, acetaldehyde, and a pericentrally cleared drug, in one continuum or with the sinusoid and the hepatocyte separated by a membrane of finite permeability. Nitric oxide is solved as a wall sub-model rather than a lobule field, because its diffusion length is 0.18 µm in the lumen. `coupling.py` returns what the network needs: a geometric shape factor (kappa = 0.996, so the single edge costs the graph nothing) and a flow-dependent resistance load from the shear-NO loop. `network_link.py` then shows how little of that survives the rest of the series: the dissipation-optimal network puts 21% of its resistance in the sinusoids, so a graft at four times the donor's flow per gram reaches a measured zR of 0.91, against the 0.31 to 0.74 of the published groups.
- `src/ipd/`: recovery of patient-level information from the published literature. `patient_level.py` carries the sixteen deaths of Ogura 2010 Table 3 with their final portal pressures, the sixteen recipients of Botha 2010 Table 5 with their portocaval gradient on both sides of a hemiportocaval shunt, and the outcome contrast of Chan 2011 Table 3 (12 events in 46 recipients, a null); every transcription is checked against the statistics its own paper reports. `ecological.py` replaces the group-mean likelihood with one that integrates over the within-group spread, and `km_reconstruct.py` implements the Guyot reconstruction, validated against simulated cohorts.

### Fixed
- Ogura 2010's printed dispersions of the portocaval gradient are standard errors, not standard deviations: the same paper gives that cohort's MELD as 23.2 ± 1.1 (n = 86) and 22.9 ± 1.7 (n = 43) while its Table 1 gives 22.9 ± 10.6. Converted, they are 2.78 and 2.62 mmHg. With that correction and two dispersions computed from a published patient table and a published range, all eight within-group standard deviations fall between 2.6 and 4.4 mmHg, which is what licenses a single pooled value of 3.4 mmHg; before it, they looked bimodal and no pooling would have been defensible.

### Notes
- The ecological correction runs one way: the patient-level slope is 1 to 11% larger than the group-mean slope over the range the series report, so the published estimate is a lower bound rather than an optimistic one.

## [0.3.1] – 2026-09-23

### Fixed
- The case export still rebuilt its own values, so an input the page rejected (a graft weight of zero or a negative one) could reach the file. `calc()` now publishes exactly what it validated and `caseRow()` reads only that; a rejected field is exported blank.
- One outcome selector wrote both `graft_loss_90d` and `death_90d` when "death or graft loss" was chosen. The three outcomes are now independent tick boxes, so each column carries what actually happened; all three stay blank while the outcome is unknown, and a CVP that was assumed rather than measured is exported blank.
- Removed two unused imports and a stale path to a file from another repository.
- The calculator mixed two models: it used the Bayesian slope (1.547) with the intercepts and bands of the stratified fit (1.959), so a reference series could show a curve and a band that disagreed. `bayes.fit()` now returns the intercept of each stratum on the absolute zP scale and the prediction band implied by the same draws, and `build_app.py` takes slope, interval, intercepts, bands and counts from that one posterior. Wang 2014 at zP 2.40 now reads 11.4% with a band of 9.5 to 17.4%, against 14.1% observed.
- `caseRow()` rebuilt the row from raw fields, so a rejected input could still be exported (an empty portal pressure gave a gradient of -5 and zP -1). It now reuses only validated values, exports blanks for anything the page rejected, and writes the 21 columns of `data/cohort_template.tsv` in order rather than 13 of its own.
- The cohort-average anchor is labelled as an approximation, and the exact option (a rate measured at one gradient) is offered first.
- The calculator was applying the exploratory slope (all outcomes combined, OR 6.02) while the paper reports the primary analysis; it now uses the primary estimate (OR 4.70 per unit of zP, 1.36 per mmHg, 7 groups from 3 series) and offers only the three SFSS series as reference levels, so page and manuscript state the same effect for the same outcome.
- `src/build_app.py` had duplicate placeholder keys, so the counts shown in the page came from a different fit than the slope.
- `THEORY.md` had been overwritten with the changelog; the theory document is restored, and a test now checks that each document is the document it claims to be.
- `index.html`, the GitHub Pages redirect, still pointed at the old repository name.

## [0.3.0] – 2026-09-22

The release that turns the repository into the model it is named after.

### The model
- **LiPNet**, the Liver Perfusion Network model, with its three normalised loads and their identity:

      zP = zF x zR        pressure load = flow load x resistance load

  `indices.zR()`, a `zR` column and a test that asserts the identity holds exactly on every group measured on both scales (0.31, 0.54, 0.55, 0.72, 0.74, all below one).
- `src/predictions.py` states the model's predictions and tests each against data it was not fitted to: P1 the shear set-point invariant across mammals, P2 the loads separating when the outflow is enlarged, P3 the within-cohort ordering, P4 three clinical thresholds meeting on one scale.
- `src/model/whole_graft.py` extends the model to a whole graft: inlet and outlet anastomoses in series, with the prediction that they fail in opposite directions (outlet stenosis raises zR above 1, inlet stenosis leaves it at 1 while both loads fall).

### The analysis
- `src/bayes.py`: the primary analysis, a hierarchical binomial model with zP centred within study, fitted by NUTS with four chains, R-hat, effective sample size, divergences and posterior predictive checks; sensitivity to the prior on tau; the Kyoto overlap sets; a conditional internal-external validation that scores the slope without fitting anything on the held-out centre; a Monte Carlo that propagates the uncertainty of the indices; a specification curve over 15 analysis choices; and a recovery simulation that shows what this design can and cannot identify.
- The superseded implementations are gone: the hand-written Metropolis sampler and the observed-to-expected validation that recalibrated the held-out intercept, which was not an independent evaluation.
- The calculator uses the estimate the paper reports, so page and manuscript cannot disagree, and carries a card that copies one case as a row for the cohort that is still missing.

### Data
- 31 groups from 17 series with provenance per cell, structural columns (centre, cohort, recruitment period, overlap set, outflow, gradient timing, outcome definition and horizon) and a `gradient_mmHg` column for the series that report the gradient itself.

### Documentation
- README rewritten around the model, with the equations as rendered mathematics, the full accounting of what enters each analysis, and the limitations confined to the prediction they belong to.
- The accompanying article is *Vascular resistance modulates portal flow-pressure decoupling in partial liver grafts*, recorded in `CITATION.cff` as the preferred citation.

## [0.2.3] – 2026-09-21

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
