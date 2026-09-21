# Liver pressure index
### Portal pressure integrates haemodynamic risk after living donor liver transplantation

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live%20on%20GitHub%20Pages-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/liver_pressure_index/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22861276.svg)](https://doi.org/10.5281/zenodo.22861276)

<!-- When the preprint is out, replace the grey badge with, e.g.:
[![medRxiv](https://img.shields.io/badge/medRxiv-10.1101%2FXXXX-b31b1b)](https://doi.org/10.1101/XXXX)
and, on acceptance:
[![DOI](https://img.shields.io/badge/J%20Hepatol-10.1016%2Fj.jhep.XXXX-0f4c81)](https://doi.org/10.1016/j.jhep.XXXX)
-->

<p align="center">
  <img src="assets/graphical_abstract.png" alt="Graphical abstract: a partial liver graft as a fraction of the donor perfusion network; the portocaval gradient normalised to its physiological value orders outcome within each centre, while the absolute level of risk is centre-specific" width="820">
</p>

Two dimensionless indices of haemodynamic load on a partial liver graft, the published series used to test them, the code that reproduces every plot as a separate file, and a calculator.

```
zF = graft PVF per 100 g / donor PVF per 100 g     (donor reference from the same series if measured, otherwise 90 mL/min/100 g)
zP = (PVP − CVP) / 5 mmHg                          (CVP = 5 when not reported; flagged in the data)
```
Both equal 1 in a healthy donor, and `zP/zF` estimates the graft's effective outflow resistance relative to the donor.

## What the published data support, and what they do not

**The absolute level of risk is not transferable between centres.** At a gradient of about 5 mmHg (zP ~ 1) the reported small-for-size syndrome (SFSS) rate is 13-17% in Uemura 2016 and about 6% in Botha 2010, while Ishizaki 2012 reports 0% at zP ~ 2.5. A single logistic model with a common intercept therefore fails: fitted on the SFSS groups it gives a slope of 0.11 (95% CI -0.25 to 0.46), and the rank correlation between zP and outcome across the main-analysis groups is weak (16 groups, rho = 0.23, p = 0.399). Flow per gram does no better (12 groups, rho = 0.27, p = 0.388; rho = 0.15, p = 0.614 with the groups defined by a flow cut-off).

**The effect of changing the gradient is transferable.** Fitted with one intercept per study and a common slope on the 11 groups from the 5 series that contribute a within-stratum contrast (Ogura 2010, Osman 2017, Uemura 2016, Wang 2014, Yagi 2005; 104 events in 734 recipients). A stratum is a study **and** an outcome definition, so a series reporting two outcomes does not share one intercept, and a cohort partitioned twice (the Wang 2014 pressure groups are the same patients as its splenectomy groups) enters only once:

```
logit(p) = alpha_stratum + 1.96 * zP     (95% CI for the slope 1.22 to 3.02; profile likelihood 1.18 to 2.83)
odds ratio 7.09 per unit of zP (3.37-20.46), i.e. 1.48 per mmHg of gradient (1.28-1.83)
```

**Removing the cohort level.** Centring within stratum (subtracting from each group the mean log-odds and the mean zP of its own cohort) cancels the intercept exactly and puts every group on one scale: the outcome odds and the gradient, each relative to that cohort's own average. All 11 groups then fall on a single line, weighted R2 = 0.91, slope 1.66 +- 0.32, **OR 5.25 per unit of zP (2.79-9.89), 1.39 per mmHg** (plot `centred`). This is the same effect obtained without estimating any intercept, and it is the cleanest statement of the result: a cohort at twice its own usual gradient has about five times its own usual odds.

Three checks support that estimate. A random-effects meta-analysis of the five per-stratum slopes gives the same answer, **OR 6.69 per unit of zP (3.12-14.4)**, with **tau2 = 0, Q = 1.47 on 4 df (p = 0.83), I2 = 0%**: the five series estimate one common effect. There is no extra-binomial spread (Pearson chi2/df = 0.31, so the standard error needs no inflation). And the aggregation into group means costs only 3% of the slope, while the within-group spread of the gradient (SD about 3.5 mmHg) implies that a per-patient analysis would find a *larger* slope, about 3.3 (attenuation factor 0.59). The slope stays between 1.67 and 2.16 when any one study is dropped, and is 1.87 for the SFSS outcomes and 2.03 for the mortality outcomes fitted separately. The intercepts differ by more than two logits between centres; that difference is recipient severity, outcome definition and technique, not haemodynamics. The model therefore predicts **how much a given reduction of the gradient changes the odds**, not the baseline risk of a patient.

## Calculator

**Live:** https://danpc11.github.io/liver_pressure_index/ (`sfss_calculator.html`, runs entirely in the browser).

Enter the current PVP and CVP, the pressure expected after the planned manoeuvre, and optionally PVF and graft weight. The calculator implements the centred result directly: it uses the intercept-free slope (OR 5.25 per unit of zP, 1.39 per mmHg) and asks for a reference level rather than assuming one. It returns zP, zF, the ratio zP/zF, and the odds ratio implied by the change of gradient with its bootstrap interval. To turn that into an absolute risk you choose a reference level: **your cohort's average** (its overall outcome rate and its average gradient, the two numbers a centre actually knows, which is exactly the anchor the centring defines), your own rate at the starting gradient, or one of the fitted strata (Ogura 2010, Osman 2017, Uemura 2016, Wang 2014, Yagi 2005), in which case the app shows what that series would go from and to, warns when the entered gradient falls outside the range that series actually spans, and the plot shows a single curve anchored on that reference, with the current and post-change gradients marked; with no anchor it shows the published groups centred on their own cohort. Example: lowering the gradient from 12 to 7 mmHg gives OR 0.14 (0.05-0.31); on a 10% baseline that is about 1.5%. A research tool, not a validated individual predictor and not a medical device.

## How to obtain your baseline

The model supplies the slope; the level has to come from the centre that uses it. This is the standard problem of **recalibration in the large**, the first and least demanding step of prediction-model updating: keep the coefficients and re-estimate only the intercept for the new setting (Steyerberg, *Clinical Prediction Models*; Vergouwe et al., Stat Med 2017;36:4529-39; Janssen et al., *A simple method to adjust clinical prediction models to local circumstances*, Can J Anesth 2009;56:194-201). It is what has been done for decades with Framingham, EuroSCORE and TAVI mortality models before local use, and geographic validation studies show exactly our pattern: heterogeneity of I2 = 0% for discrimination and the slope, but 74% for the calibration intercept (Jacobs et al., J Clin Epidemiol 2023).

Three ways to get it, in increasing order of effort:

1. **From two numbers you already have** (what the calculator does by default). With your cohort's overall outcome rate p and its average gradient, the intercept is `alpha = logit(p) - beta * mean(zP)`, and the risk of a patient is `logistic(alpha + beta * zP)`. This is calibration in the large computed from the event rate alone, and it needs no individual data.
2. **From a retrospective audit.** With 50-100 consecutive recipients and their final PVP, CVP and outcome, fit only the intercept with `beta` held fixed as an offset (`indices.recalibrate()`); this also gives a standard error for the level and lets you check the slope in your own data.
3. **From a prospective cohort.** Fill `data/cohort_template.tsv` and refit intercept and slope; the closed testing procedure of Vergouwe et al. tells you when your sample justifies re-estimating more than the intercept.

Two cautions from the same literature. Recalibrating the intercept corrects the average level but does not repair discrimination if the case mix differs (Debray et al., J Clin Epidemiol 2015;68:279-89), and recalibration should not be reflexive: if the level differs because of something the model omits, updating the intercept hides that rather than fixing it (Van Calster et al., 2025).

## Contents

```
sfss_calculator.html     the calculator (also served at the Pages URL above)
index.html               redirects the Pages root to the calculator
assets/                  graphical abstract and other images for this README
data/series.tsv          31 groups from 17 series, raw values with provenance columns
data/patients.tsv        individual patients where the papers tabulate them (Ou 2010, Yamada 2008, Alim 2016)
data/interventions.tsv   pressure and flow changes produced by each inflow-modulation manoeuvre, with source
data/cohort_template.tsv column layout for a future individual-patient cohort
data/graph*.tsv, certificate_small.tsv   cached outputs of the long model campaigns
src/indices.py           indices, correlations, within-study model, meta-analysis of the stratum slopes, centring that removes the cohort level, overdispersion and regression-dilution checks, risk_after()
src/plots.py             17 plots, each saved on its own as results/<name>.pdf and .png
src/build_app.py         regenerates sfss_calculator.html from the fitted effect
src/app_template.html    the calculator with placeholders
src/model/               the perfusion-network model (2D, 3D, closed-form scaling) and its campaigns
tests/                   index definitions, reproducibility of the fit, plots and app build (run by CI)
run_all.sh               indices -> plots -> app
```

## Reproduce

```bash
pip install -r requirements.txt
./run_all.sh                       # everything, ~2 min
python src/plots.py risk_curve     # or any single plot
python src/plots.py                # prints the list of plots
python -m pytest tests -q          # what the CI runs
```

Plots: `nomogram`, `within_study`, `forest`, `centred`, `risk_change`, `series_pressure`, `series_flow`, `plane`, `resection`, `load_curve`, `interventions`, `network_2d`, `network_3d`, `scaling`, `allometry`, `shear_profile`, `certificate`.

## Data

`data/series.tsv` holds **raw values and their provenance only**, including a `gradient_mmHg` column for the series that report the portocaval gradient itself (Uemura 2016, Ishizaki 2012, Botha 2010, Ogura 2010, Chan 2011, Yamada 2008); `zP` uses that gradient when present and PVP − CVP otherwise; the indices are recomputed by `src/indices.py`, which is the single source used by the plots, the tests and the calculator. Columns: `study`, `group`, `n`, `source` (table or section and page in the original PDF), `GRWR`, `PVF_per_100g`, `PVF_basis`, `donor_PVF_per_100g_ref`, `PVP`, `PVP_basis`, `CVP`, `CVP_basis`, `CVP_measured`, `outcome`, `outcome_type`, `events`, `events_basis`, `pct`, `notes`.

Each `*_basis` column states how the value was obtained:

| basis | meaning | where it occurs |
|---|---|---|
| `reported` | copied from a table or the text | most values |
| `derived` | computed from reported values (e.g. PVF divided by graft weight; PVP = gradient + CVP; deaths from a survival percentage) | Yagi 2006 PVF; Chan 2011 CVP; Ogura 2010 PVP and deaths; Yagi 2005 deaths |
| `imputed` | not reported; a representative value was assumed | CVP = 5 mmHg in Wang 2014, Osman 2017, Yagi 2005; PVP 21 and 16 for the Wang 2014 groups defined by a 20 mmHg cut |
| `threshold` / `group mean` | group defined by a cut-off; the cut-off or the mean of the complementary group stands in for the group value. **Excluded from the main analysis** (`in_main_analysis = False`), shown hollow in the plots and reported in sensitivity | Vasavada 2014 |
| `not applicable` | the series did not measure that variable | flow-only or pressure-only series |

`indices.compute()` records every value it fills in itself (`CVP_filled`, `donor_ref_filled`) and treats it as imputed regardless of the label, and it raises an error if a label contradicts the data (e.g. a blank CVP marked `reported`).

**Sensitivity** (`results/sensitivity.tsv`): main analysis, zP 16 groups rho = 0.23 (p = 0.399), zF 12 groups rho = 0.27 (p = 0.388); including the cut-off groups, zP 21 groups rho = 0.34 (p = 0.137) and zF 14 groups rho = 0.15 (p = 0.614); restricted to groups with no imputed or filled value in the index itself, zP 10 groups rho = -0.28 (p = 0.434). The within-study model uses only measured or derived gradients and is the result the analysis supports.

Series: Troisi 2003, Troisi 2005, Ou 2010, Vasavada 2014, Alim 2016, Chan 2011, Yagi 2005, Yagi 2006, Wang 2014, Osman 2017, Ogura 2010, Yamada 2008, Uemura 2016 (Surgery 159:1623), Yao 2018 (Transplantation 102:623), Kanetkar 2017 (J Clin Exp Hepatol 7:235), Ishizaki 2012 (Liver Transpl 18:305), Botha 2010 (Liver Transpl 16:649). Several come from the same centre and overlapping periods (Kyoto: Yagi 2005, Yagi 2006, Ogura 2010, Uemura 2016, Yao 2018), so the groups are not independent observations and the confidence intervals do not account for that.

## Adding a cohort

Fill `data/cohort_template.tsv` (one row per patient). `src/indices.py` exposes `zF()`, `zP()`, `within_study_fit()` (which takes a `study` column, so a new cohort enters as one more stratum) and `risk_after(baseline_rate, delta_zP, beta)`. `python src/build_app.py` then updates the calculator with the new effect. With individual-patient data the same model can be fitted with a random intercept per centre and adjusted for severity, which is the analysis this repository is meant to enable.

## Model

`src/model/network.py` builds the liver as a graph (hilum → portal tree → hexagonal lobules → hepatic venous tree), solves flows and dissipation for prescribed lobule demand and finds the dissipation-optimal tree at fixed maintenance cost `Σ m_e^b`. `network3d.py` does the same on a hemisphere. `scaling.py` gives the closed form `D* ∝ F² C^(−2/b) N^((2−b)/b) L³`, the set-point `τ0 = √(D*/C)` and the allometric closure. `graph_scaling.py`, `network3d.py <R>` and `certify.py` are the campaigns whose outputs are cached in `data/`.

## Versions

| version | date | changes |
|---|---|---|
| 0.3.0 | 2026-09-21 | Random-effects meta-analysis of the stratum slopes (I² = 0%), overdispersion and regression-dilution checks, forest plot. See [CHANGELOG.md](CHANGELOG.md) |
| 0.2.2 | 2026-09-21 | Audit fixes: a cohort counted twice removed, intercepts per study × outcome, gradient column, convergence and profile-likelihood checks; effect OR 7.09 per unit z_P |
| 0.2.1 | 2026-09-21 | Five series with measured gradients; within-study model replaces the pooled fit; calculator rebuilt around the change of gradient |
| 0.2.0 | 2026-09-20 | Provenance columns, single-source indices, calculator input validation |
| 0.1 | 2026-09-20 | First public release: extraction table of 12 series, indices, pooled fit, 14 plots, calculator, GitHub Pages, tests, PolyForm Noncommercial licence. Archived at Zenodo, version DOI [10.5281/zenodo.22861277](https://doi.org/10.5281/zenodo.22861277) |

Every GitHub release is archived at Zenodo. The concept DOI [10.5281/zenodo.22861276](https://doi.org/10.5281/zenodo.22861276) always resolves to the latest version; each version has its own DOI (table above).

## Citation

See `CITATION.cff` (GitHub shows a "Cite this repository" button). Until the preprint appears, cite the archived software: Vázquez-Victorio G, Pérez-Calixto M, Escutia-Guadarrama L, Cervera A, Tovar H, Pérez-Calixto D. liver_pressure_index: portal pressure and flow indices of haemodynamic load on partial liver grafts (v0.1). Zenodo, 2026. https://doi.org/10.5281/zenodo.22861276


## Licence

[PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0): free to use, modify and share for research, teaching, personal and non-profit purposes (universities, public research organisations and hospitals included, whatever their funding source); any commercial use requires a separate licence from the authors. Values in `data/` are transcribed from the cited publications and remain the property of their authors. The calculator is a research tool, not a medical device.
