# Liver pressure index
### A perfusion-network model of the liver, and the predictions it makes about partial grafts

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/liver_pressure_index/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22861276.svg)](https://doi.org/10.5281/zenodo.22861276)

## The model

The liver is treated as a transport network: a portal tree feeding a very large number of near-identical lobules, and a hepatic venous tree collecting from them. Minimising dissipated power at a fixed vascular maintenance cost `Σ m_e^b` gives, with no free parameters,

```
D* ∝ F² C^(−2/b) N^((2−b)/b) L³        the optimum, a space-filling tree with one sinusoid per lobule
τ0 = √(D*/C)                            a wall shear set-point fixed by the whole network
ΔP = f · R                              the pressure drop across a lobule
```

A partial graft is the fraction `g` of that network that remains, perfused with `h` times the donor inflow, so each lobule receives `h/g` times its donor flow. Written in measurable quantities:

```
zF = graft PVF per 100 g / donor PVF per 100 g      flow per lobule, relative to the donor
zP = (PVP − CVP) / 5 mmHg                           sinusoidal pressure, relative to normal
zP = zF · (R_graft / R_donor)                       they are the same variable only at donor outflow resistance
```

The model contains haemodynamics and nothing else: no recipient severity, no donor age, no steatosis. That is what makes it falsifiable.

<p align="center">
  <img src="assets/graphical_abstract.png" alt="Graphical abstract: the liver as a dissipation-optimal perfusion network; a partial graft is a fraction of that network, and the portocaval gradient normalised to its physiological value is the load per lobule" width="820">
</p>

## Four predictions, and what the data say

Run `python src/predictions.py`; every number below comes from `results/predictions.json`.

**P1. The shear set-point must be invariant across mammals.** Portal pressure is 6–11 mmHg from mouse to human, so if the network is at its optimum, τ0 cannot scale with body mass. With the measured scalings (portal flow ~ M^0.77, portocentral distance ~ M^0.10), the predicted exponent of τ0 is within ±0.08 of zero for every maintenance exponent, and the predicted vascular mass exponent is 0.89–0.90 at b = 0.85–1.0 against the 0.86 observed for hepatic blood volume. This uses no transplant data at all.

**P2. Flow and pressure must separate when the outflow is enlarged.** Five published groups report both indices. The four whose outflow was reconstructed (middle hepatic vein included, venoplasty, left lobe with the middle and left hepatic vein trunk, splenectomy) have `zP/zF` of 0.31, 0.55, 0.54 and 0.74, all below one as predicted, and outcomes of 0–10% despite flows of three to four times the donor's. No flow threshold and no pressure threshold predicts this; the ratio does.

**P3. Within a cohort, outcome must rise with zP, with a common slope and a centre-specific level.** The primary analysis is a hierarchical binomial model with the pressure index centred within each study, so a study's level and its association are not traded off against each other:

```
E_gs ~ Binomial(n_gs, p_gs),  logit(p_gs) = alpha_s + beta_s (zP_gs − mean zP_s),  beta_s ~ N(mu_beta, tau_beta²)
```

NUTS, four chains (R-hat 1.001, minimum effective size 5114, 0 divergences, posterior predictive check 7/7). In the prespecified order: primary outcome, SFSS or early dysfunction, 7 groups from 3 studies, **mu_beta +1.55 (95% CrI +0.33 to +2.79), OR 4.70 per unit of zP, P(mu > 0) = 0.990**; secondary outcome, mortality or graft loss, +1.74 (+0.60, +2.86); both combined, +1.79 (+0.97, +2.64). Every within-study pair moves in the predicted direction, and the intercepts span almost four logits, which is what the model implies since it says nothing about what sets a centre's baseline.

Across 15 specifications (outcome definition, imputation, one Kyoto publication at a time, assumed CVP from 3 to 9 mmHg, three priors on tau, one stratum per centre instead of per publication) mu_beta stays between +0.92 and +1.93, always positive. The smallest value is the one that pools each centre's publications into a single stratum, +0.92 (+0.20, +1.65), which is the most conservative way to handle the overlapping Kyoto series. The prior on the between-study variation matters, as it must with five studies: HalfNormal(0.25) gives +1.60 with P(mu>0) = 0.998, HalfNormal(1.0) gives +1.46 (-0.07, +2.84) with 0.970. Propagating the uncertainty of the indices themselves, with the assumed CVP drawn per study over 3-9 mmHg and the posteriors pooled rather than their means, gives +1.74 (+0.89, +2.63).

Validation leaves out one centre, estimates the slope on the others and scores the held-out centre's within-cohort contrast with the conditional likelihood, which cancels its intercept, so nothing is fitted on its outcomes: the log-score gain over no association is +2.26 (Cairo), +2.74 (Fukuoka) and +8.40 (Kyoto), with the direction right in all three. A recovery simulation on the real sizes and gradients shows the limits of the design: with no true effect it never excludes zero, with mu = 1 it does about 60% of the time, with mu = 2 always, and tau_beta is not identifiable with five studies.

**P4. Thresholds from three different fields must coincide on this scale.** Clinically significant portal hypertension (HVPG ≥ 10 mmHg), the graft threshold (PVP 15 mmHg at CVP 5) and the limit above which a cirrhotic remnant is pushed after resection are all zP = 2; variceal bleeding (HVPG ≥ 12) is 2.4. Three thresholds derived independently, in cirrhosis, transplantation and hepatic surgery, land on the same normalised value.

## What this is not

An individual risk model. The level of risk is centre-specific: at a gradient near 5 mmHg the reported small-for-size rate is 13–17% in one series and 6% in another, and a third reports none at 12 mmHg. A single logistic curve fitted across centres is flat (slope 0.11, 95% CI −0.25 to 0.46), and the rank correlation of zP with outcome across series is weak (16 groups, ρ = 0.23, p = 0.399), as is that of flow per gram (12 groups, ρ = 0.27, p = 0.388). The transferable quantity is the effect of changing the gradient, not the baseline.

## Calculator

**https://danpc11.github.io/liver_pressure_index/** — one HTML file, runs in the browser. Enter the current PVP and CVP, the pressure expected after the planned manoeuvre, and an anchor for the level (your cohort's overall rate and average gradient, your rate at the current gradient, or one of the fitted series). It returns zP, zF, their ratio, the odds ratio for the change and the absolute risks that anchor implies. Slope and intercepts come from the same model, so the curve passes through the anchor.

To obtain your own baseline, `indices.baseline_from_rate(rate, mean_zP, beta)` gives `alpha = logit(rate) − beta·mean(zP)`, and `indices.recalibrate(outcomes, zP, beta)` fits the intercept on an audit with the slope held fixed. This is recalibration in the large, the first step of prediction-model updating (Steyerberg; Vergouwe et al., Stat Med 2017;36:4529; Janssen et al., Can J Anesth 2009;56:194).

## Contents

```
src/model/          the network model: 2D and 3D graphs, the optimum, closed-form scaling, allometric closure
src/predictions.py  tests P1-P4 against the data -> results/predictions.json
src/indices.py      zF, zP, stratified and hierarchical fits, internal-external validation, recalibration
src/plots.py        18 stand-alone plots, one file each
src/build_app.py    rebuilds the calculator from the fitted effect
data/series.tsv     31 groups from 17 series: raw values, provenance, centre, cohort, period, outflow, timing
sfss_calculator.html, index.html, tests/, run_all.sh
```

## Reproduce

```bash
pip install -r requirements.txt
./run_all.sh                       # indices, predictions, plots and calculator, about three minutes
python src/predictions.py          # the four predictions with their numbers
python src/plots.py predictions     # the four-panel summary figure
```

Model campaigns behind the cached tables (the optimum on 37–1478 lobules in 2D and 3D, the exhaustive certificate on seven lobules, the allometric closure) are in `src/model/` and listed, commented, in `run_all.sh`.

## Data

`data/series.tsv` holds raw values only; indices are recomputed by `src/indices.py`. Each `*_basis` column says where a value came from (`reported`, `derived`, `imputed`, `threshold`, `not applicable`), and the structural columns (`centre_id`, `cohort_id`, `recruitment_start/end`, `overlap_set`, `outflow`, `gradient_timing`, `gradient_sd`, `outcome_definition`, `outcome_horizon`, `modulation_strategy`) separate study, centre and cohort. Groups defined only by a cut-off, and the second partition of a cohort already counted, are excluded from the main analysis and reported in `results/sensitivity.tsv`.

Series: Troisi 2003, Troisi 2005, Ou 2010, Vasavada 2014, Alim 2016, Chan 2011, Yagi 2005, Yagi 2006, Wang 2014, Osman 2017, Ogura 2010, Yamada 2008, Uemura 2016, Yao 2018, Kanetkar 2017, Ishizaki 2012, Botha 2010.

## Limitations

P3 rests on aggregated group data: five series from four centres, 104 events, three of which report the primary outcome, with Kyoto series overlapping in time. Assumed CVP does not affect the slope (a constant shift within a study is absorbed by its intercept), but it does move where groups sit on the axis. Aggregated data allow the likelihood, the deviance and observed-to-expected ratios; they cannot give a c-statistic, an individual calibration curve, a Brier score or a decision curve. P1, P2 and P4 do not depend on that analysis. The calculator is a research tool, not a medical device.

## Citation and licence

`CITATION.cff`; until the preprint appears, cite the Zenodo record [10.5281/zenodo.22861276](https://doi.org/10.5281/zenodo.22861276). Code under [PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0); values in `data/` belong to the authors of the cited papers. Funded by DGAPA-PAPIIT IN234029 and SECIHTI CBF-2025-G-789.
