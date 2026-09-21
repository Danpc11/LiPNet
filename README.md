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

Two dimensionless indices of haemodynamic load on a partial liver graft, the published series used to test them, the code that reproduces every plot as a separate file, and a calculator.

```
zF = graft PVF per 100 g / donor PVF per 100 g     (donor reference from the same series if measured, otherwise 90 mL/min/100 g)
zP = (PVP − CVP) / 5 mmHg                          (CVP = 5 when not reported; flagged in the data)
```
Both equal 1 in a healthy donor, and `zP/zF` estimates the graft's effective outflow resistance relative to the donor.

## What the published data support, and what they do not

**The absolute level of risk is not transferable between centres.** At a gradient of about 5 mmHg (zP ~ 1) the reported small-for-size syndrome (SFSS) rate is 13-17% in Uemura 2016 and about 6% in Botha 2010, while Ishizaki 2012 reports 0% at zP ~ 2.5. A single logistic model with a common intercept therefore fails: fitted on the SFSS groups it gives a slope of 0.13 (95% CI -0.22 to 0.49), and the rank correlation between zP and outcome across the main-analysis groups is weak (18 groups, rho = 0.31, p = 0.214). Flow per gram does no better (11 groups, rho = 0.35, p = 0.290; rho = 0.23, p = 0.459 with the two groups defined by a flow cut-off).

**The effect of changing the gradient is transferable.** Fitted with one intercept per study and a common slope on the 13 groups from the 5 series that contribute a within-study contrast (Ogura 2010, Osman 2017, Uemura 2016, Wang 2014, Yagi 2005; 134 events in 1,026 recipients):

```
logit(p) = alpha_study + 1.34 * zP        (95% CI for the slope 0.89 to 1.85)
odds ratio 3.81 per unit of zP (2.44-6.37), i.e. 1.31 per mmHg of gradient (1.20-1.45)
```

The slope stays between 1.13 and 2.03 when any one study is dropped. The intercepts differ by more than two logits between centres; that difference is recipient severity, outcome definition and technique, not haemodynamics. The model therefore predicts **how much a given reduction of the gradient changes the odds**, not the baseline risk of a patient.

## Calculator

**Live:** https://danpc11.github.io/liver_pressure_index/ (`sfss_calculator.html`, runs entirely in the browser).

Enter the current PVP and CVP, the pressure expected after the planned manoeuvre, and optionally PVF and graft weight. It returns zP, zF, the ratio zP/zF, the odds ratio implied by the change of gradient with its bootstrap interval, and, only if you supply **your own** outcome rate at the starting gradient, the corresponding absolute risk. Example: lowering the gradient from 12 to 7 mmHg gives OR 0.26 (0.16-0.41); on a 10% baseline that is about 2.8%. A research tool, not a validated individual predictor and not a medical device.

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
src/indices.py           indices, correlations, within-study model, pooled model (for the record), risk_after()
src/plots.py             15 plots, each saved on its own as results/<name>.pdf and .png
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

Plots: `nomogram`, `within_study`, `risk_change`, `series_pressure`, `series_flow`, `plane`, `resection`, `load_curve`, `interventions`, `network_2d`, `network_3d`, `scaling`, `allometry`, `shear_profile`, `certificate`.

## Data

`data/series.tsv` holds **raw values and their provenance only**; the indices are recomputed by `src/indices.py`, which is the single source used by the plots, the tests and the calculator. Columns: `study`, `group`, `n`, `source` (table or section and page in the original PDF), `GRWR`, `PVF_per_100g`, `PVF_basis`, `donor_PVF_per_100g_ref`, `PVP`, `PVP_basis`, `CVP`, `CVP_basis`, `CVP_measured`, `outcome`, `outcome_type`, `events`, `events_basis`, `pct`, `notes`.

Each `*_basis` column states how the value was obtained:

| basis | meaning | where it occurs |
|---|---|---|
| `reported` | copied from a table or the text | most values |
| `derived` | computed from reported values (e.g. PVF divided by graft weight; PVP = gradient + CVP; deaths from a survival percentage) | Yagi 2006 PVF; Chan 2011 CVP; Ogura 2010 PVP and deaths; Yagi 2005 deaths |
| `imputed` | not reported; a representative value was assumed | CVP = 5 mmHg in Wang 2014, Osman 2017, Yagi 2005; PVP 21 and 16 for the Wang 2014 groups defined by a 20 mmHg cut |
| `threshold` / `group mean` | group defined by a cut-off; the cut-off or the mean of the complementary group stands in for the group value. **Excluded from the main analysis** (`in_main_analysis = False`), shown hollow in the plots and reported in sensitivity | Vasavada 2014 |
| `not applicable` | the series did not measure that variable | flow-only or pressure-only series |

`indices.compute()` records every value it fills in itself (`CVP_filled`, `donor_ref_filled`) and treats it as imputed regardless of the label, and it raises an error if a label contradicts the data (e.g. a blank CVP marked `reported`).

**Sensitivity** (`results/sensitivity.tsv`): main analysis, zP 18 groups rho = 0.31 (p = 0.214), zF 11 groups rho = 0.35 (p = 0.290); including the cut-off groups, zP 21 groups rho = 0.34 (p = 0.137) and zF 13 groups rho = 0.23 (p = 0.459); restricted to groups with no imputed or filled value in the index itself, zP 10 groups rho = -0.28 (p = 0.434). The within-study model uses only measured or derived gradients and is the result the analysis supports.

Series: Troisi 2003, Troisi 2005, Ou 2010, Vasavada 2014, Alim 2016, Chan 2011, Yagi 2005, Yagi 2006, Wang 2014, Osman 2017, Ogura 2010, Yamada 2008, Uemura 2016 (Surgery 159:1623), Yao 2018 (Transplantation 102:623), Kanetkar 2017 (J Clin Exp Hepatol 7:235), Ishizaki 2012 (Liver Transpl 18:305), Botha 2010 (Liver Transpl 16:649). Several come from the same centre and overlapping periods (Kyoto: Yagi 2005, Yagi 2006, Ogura 2010, Uemura 2016, Yao 2018), so the groups are not independent observations and the confidence intervals do not account for that.

## Adding a cohort

Fill `data/cohort_template.tsv` (one row per patient). `src/indices.py` exposes `zF()`, `zP()`, `within_study_fit()` (which takes a `study` column, so a new cohort enters as one more stratum) and `risk_after(baseline_rate, delta_zP, beta)`. `python src/build_app.py` then updates the calculator with the new effect. With individual-patient data the same model can be fitted with a random intercept per centre and adjusted for severity, which is the analysis this repository is meant to enable.

## Model

`src/model/network.py` builds the liver as a graph (hilum → portal tree → hexagonal lobules → hepatic venous tree), solves flows and dissipation for prescribed lobule demand and finds the dissipation-optimal tree at fixed maintenance cost `Σ m_e^b`. `network3d.py` does the same on a hemisphere. `scaling.py` gives the closed form `D* ∝ F² C^(−2/b) N^((2−b)/b) L³`, the set-point `τ0 = √(D*/C)` and the allometric closure. `graph_scaling.py`, `network3d.py <R>` and `certify.py` are the campaigns whose outputs are cached in `data/`.

## Versions

| version | date | changes |
|---|---|---|
| 0.1 | 2026-09-20 | First public release: extraction table of 12 series, indices, pooled fit, 14 plots, calculator, GitHub Pages, tests, PolyForm Noncommercial licence. Archived at Zenodo, version DOI [10.5281/zenodo.22861277](https://doi.org/10.5281/zenodo.22861277) |

Every GitHub release is archived at Zenodo. The concept DOI [10.5281/zenodo.22861276](https://doi.org/10.5281/zenodo.22861276) always resolves to the latest version; each version has its own DOI (table above).

## Citation

See `CITATION.cff` (GitHub shows a "Cite this repository" button). Until the preprint appears, cite the archived software: Vázquez-Victorio G, Pérez-Calixto M, Escutia-Guadarrama L, Cervera A, Tovar H, Pérez-Calixto D. liver_pressure_index: portal pressure and flow indices of haemodynamic load on partial liver grafts (v0.1). Zenodo, 2026. https://doi.org/10.5281/zenodo.22861276


## Licence

[PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0): free to use, modify and share for research, teaching, personal and non-profit purposes (universities, public research organisations and hospitals included, whatever their funding source); any commercial use requires a separate licence from the authors. Values in `data/` are transcribed from the cited publications and remain the property of their authors. The calculator is a research tool, not a medical device.
