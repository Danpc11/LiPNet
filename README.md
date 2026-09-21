# Liver pressure index
### Portal pressure and haemodynamic load on a partial liver graft

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/liver_pressure_index/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22861276.svg)](https://doi.org/10.5281/zenodo.22861276)

<!-- preprint badge goes here when it is out -->
Two donor-normalised indices of the haemodynamic load on a partial liver graft, the published series used to test them, and a calculator.

```
zF = graft PVF per 100 g / donor PVF per 100 g     (donor reference from the same series, else 90 mL/min/100 g)
zP = (PVP − CVP) / 5 mmHg                          (CVP = 5 when not reported)
```

Both equal 1 in a healthy donor. Their ratio, `zP/zF`, estimates the graft's outflow resistance relative to the donor, which is why flow and pressure can disagree.


<p align="center">
  <img src="assets/graphical_abstract.png" alt="Graphical abstract: a partial liver graft as a fraction of the donor perfusion network; within each centre the portocaval gradient tracks outcome, while the level of risk is centre-specific" width="520">
</p>

## What the data show

**The level of risk does not travel.** At a gradient near 5 mmHg the reported small-for-size rate is 13–17% in Uemura 2016 and about 6% in Botha 2010, while Ishizaki 2012 reports none at 12 mmHg. A single logistic curve fitted across centres is therefore flat (slope 0.11, 95% CI −0.25 to 0.46), and the rank correlation between `zP` and outcome across series is weak (16 groups, ρ = 0.23, p = 0.399). Flow per gram does no better (12 groups, ρ = 0.27, p = 0.388).

**The effect of changing the gradient does travel.** Fitted with one intercept per stratum (study and outcome definition) and a common slope, on eleven groups from five series, 104 events in 734 recipients:

```
OR 7.09 per unit of zP (95% CI 3.37–20.5)   =   1.48 per mmHg of gradient
```

A random-effects meta-analysis of the five per-series slopes agrees (OR 6.69, 95% CI 3.12–14.4) and finds no heterogeneity (τ² = 0, I² = 0%, Q = 1.47 on 4 df, p = 0.83). The slope holds for both outcome definitions (1.87 and 2.03) and when any study is dropped (1.67–2.16). Centring each cohort on its own average puts all eleven groups on one scale (plot `centred`, weighted R² = 0.91); that view is descriptive, the effect comes from the stratified fit. This is the familiar behaviour of prediction models across sites: the coefficients transport, the intercept does not.

## Calculator

**https://danpc11.github.io/liver_pressure_index/** — a single HTML file, runs in the browser.

Enter the current PVP and CVP, the pressure expected after the planned manoeuvre, and an anchor for the level: your cohort's overall outcome rate and average gradient, your rate at the current gradient, or one of the fitted series. It returns `zP`, `zF`, their ratio, the odds ratio for the change, and the absolute risks that anchor implies. Slope, interval and intercepts all come from the same stratified model, so the curve passes through the anchor you give it.

## Getting your own baseline

The slope comes from the model; the level has to come from you. Re-estimating only the intercept for a new setting is *recalibration in the large*, the first step of prediction-model updating (Steyerberg, *Clinical Prediction Models*; Vergouwe et al., Stat Med 2017;36:4529; Janssen et al., Can J Anesth 2009;56:194). In a geographic validation of TAVI mortality models across 16 hospitals the heterogeneity was 0% for the slope and 74% for the intercept — the same split we see here (Jacobs et al., J Clin Epidemiol 2023).

- From two numbers: `baseline_from_rate(rate, mean_zP, beta)` gives `alpha = logit(rate) − beta·mean(zP)`. No individual data needed.
- From an audit of 50–100 recipients: `recalibrate(outcomes, zP, beta)` fits the intercept with the slope held fixed, and returns its standard error.
- From a prospective cohort: fill `data/cohort_template.tsv` and refit both.

## Contents

```
sfss_calculator.html   the calculator          data/series.tsv   31 groups from 17 series, raw values with provenance
index.html             redirect for Pages      data/patients.tsv individual patients where the papers tabulate them
src/indices.py         indices and models      data/interventions.tsv  effect of each inflow-modulation manoeuvre
src/plots.py           17 stand-alone plots    data/cohort_template.tsv  columns for a new cohort
src/build_app.py       rebuilds the calculator src/model/        the perfusion-network model behind the indices
```

## Reproduce

```bash
pip install -r requirements.txt
./run_all.sh                     # indices, plots and calculator, about two minutes
python src/plots.py centred      # or any single plot; with no argument it lists them
```

## Data

`data/series.tsv` holds raw values only; the indices are recomputed by `src/indices.py`. Each `*_basis` column says where the value came from: `reported`, `derived` (e.g. PVP from a reported gradient), `imputed` (CVP assumed at 5 mmHg), `threshold` (the group is defined only by a cut-off) or `not applicable`. Groups defined only by a cut-off, and the second partition of a cohort already counted, are excluded from the main analysis and reported in `results/sensitivity.tsv`.

Series: Troisi 2003, Troisi 2005, Ou 2010, Vasavada 2014, Alim 2016, Chan 2011, Yagi 2005, Yagi 2006, Wang 2014, Osman 2017, Ogura 2010, Yamada 2008, Uemura 2016, Yao 2018, Kanetkar 2017, Ishizaki 2012, Botha 2010.

## Limitations

Five series from four centres, 104 events. Several groups come from Kyoto over overlapping periods, so they are not independent. CVP was assumed in three series. Everything here is group-level: it shows how outcome moves with the gradient inside a cohort, not what a given patient's risk is. The calculator is a research tool, not a medical device.

## Citation and licence

`CITATION.cff`; until the preprint appears, cite the Zenodo record [10.5281/zenodo.22861276](https://doi.org/10.5281/zenodo.22861276). Code under [PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0); the values in `data/` belong to the authors of the cited papers. Funded by DGAPA-PAPIIT IN234029 and SECIHTI CBF-2025-G-789.
