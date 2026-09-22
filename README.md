# LiPNet - Liver Perfusion Network model

### A physics-based model of hepatic perfusion and haemodynamic load in partial liver grafts

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/LiPNet/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/LiPNet/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/LiPNet/)
![Version](https://img.shields.io/badge/version-0.3.0-1f6feb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22861276.svg)](https://doi.org/10.5281/zenodo.22861276)

---

## Overview

Partial liver grafts receive portal blood through a smaller amount of tissue than a whole liver.

As a result, portal flow per unit of liver mass can become several times higher than in the donor.

However, high portal flow does not always produce equally high portal pressure.

Some grafts can carry high flow while maintaining relatively low portal pressure, particularly when venous outflow is enlarged or reconstructed.

This raises a simple question:

> **What determines how portal flow is converted into portal pressure?**

LiPNet approaches this problem by treating the liver as a vascular transport network.

The central result is:

$$\boxed{z_P = z_F \, z_R}$$ or, in words, $$\boxed{\text{pressure load}=\text{flow load}\times\text{resistance load}}$$

The model therefore separates three components of hepatic haemodynamics:

- portal flow,
- vascular resistance,
- portal pressure.

This provides a simple framework for understanding why flow and pressure may become uncoupled in partial liver grafts.

---

## The model

The liver is represented as a transport network with three main components:

1. a portal venous tree that delivers blood,
2. a large number of approximately similar liver lobules,
3. a hepatic venous tree that collects blood.

The portal and hepatic venous trees form a **canopy-to-canopy transport network**, with the liver lobules acting as the functional units between both vascular trees.

The network is assumed to operate near an energetic optimum.

We minimise dissipated power while keeping vascular maintenance cost fixed:

$$\sum_e m_e^b$$

where:

- $m_e$ is the vascular material associated with vessel segment $e$,
- $b$ is the vascular maintenance exponent.

This gives:

$$D^* \propto F^{2}\,C^{-2/b}\,N^{(2-b)/b}\,L^{3},$$

together with a characteristic wall-shear set-point:

$$\tau_0=\sqrt{D^*/C},$$

and the pressure drop across a lobule:

$$\Delta P=fR.$$

Here:

- $D^*$ is the minimum dissipation of the network,
- $F$ is total flow,
- $C$ represents vascular maintenance cost,
- $N$ is the number of lobules,
- $L$ is a characteristic transport length,
- $\tau_0$ is the wall-shear set-point,
- $f$ is lobular flow,
- $R$ is vascular resistance.

The optimum corresponds to a space-filling vascular network in which portal blood reaches the lobular bed and then drains through the hepatic venous network.

---

## Partial grafts

Consider a partial graft that contains a fraction $g$ of the donor liver network.

If the graft receives $h$ times the donor inflow, each remaining lobule receives approximately:

$$\frac{h}{g}$$

times its original donor flow.

This is the basic reason why small grafts can experience very high haemodynamic load.

For example, a graft containing half of the original liver mass but receiving the full donor inflow has approximately twice the flow per remaining lobule.

LiPNet expresses this behaviour using measurable clinical quantities.

---

## Three normalized haemodynamic loads

### Flow load

The normalized portal flow load is:

$$z_F=\frac{\text{graft PVF per }100\,\text{g}}{\text{donor PVF per }100\,\text{g}}.$$

Here, PVF is portal venous flow.

A value of:

$$z_F=1$$

means that the graft receives the same portal flow per unit mass as the healthy donor liver.

A value of:

$$z_F=3$$

means that each unit of graft tissue receives approximately three times the donor flow.

---

### Pressure load

The normalized pressure load is:

$$z_P=\frac{\text{PVP}-\text{CVP}}{5\ \text{mmHg}}.$$

where:

- PVP is portal venous pressure,
- CVP is central venous pressure.

The denominator of 5 mmHg represents a normal portal-to-central venous pressure gradient.

Therefore:

$$z_P=1$$

represents the normal donor state.

---

### Resistance load

The normalized resistance load is:

$$z_R=\frac{z_P}{z_F}=\frac{R_{\text{graft}}}{R_{\text{donor}}}.$$

This gives the central identity of LiPNet:

$$\boxed{z_P = z_F \, z_R}$$

All three quantities equal 1 in a healthy donor:

$$z_F=z_P=z_R=1.$$

This identity is the main point of the model.

Flow and pressure are equivalent normalized loads only when:

$$z_R=1.$$

When graft resistance changes, flow and pressure separate.

---

### Why resistance matters

Suppose a partial graft receives three times the donor flow:

$$z_F=3.$$

If resistance remains equal to the donor value:

$$z_R=1,$$

then:

$$z_P=3.$$

But if venous reconstruction reduces effective resistance to:

$$z_R=0.5,$$

then:

$$z_P=3\times0.5=1.5.$$

The graft therefore carries very high flow without developing the pressure expected from flow alone.

This is the mechanism behind **portal flow-pressure decoupling** in LiPNet.

Any intervention that enlarges effective outflow can lower $z_R$.

Examples include:

- middle hepatic vein reconstruction,
- venoplasty,
- enlarged hepatic venous drainage,
- haemodynamic interventions that increase effective compliance.

The model therefore interprets pressure as the combined result of flow and resistance.

---

## What the model includes

LiPNet contains haemodynamics only.

It does not include:

- recipient severity,
- donor age,
- steatosis,
- inflammation,
- graft quality,
- immunological injury,
- metabolic dysfunction.

This is intentional.

The model is designed to isolate the physical contribution of portal flow, vascular resistance and portal pressure.

This also makes its predictions falsifiable.

<p align="center">
  <img src="assets/graphical_abstract.png"
       alt="Graphical abstract of LiPNet showing the liver as a perfusion network and the relation between flow load, resistance load and pressure load"
       width="520">
</p>

---

## Four predictions

All numerical results below are produced by:

```bash
python src/predictions.py
```

and stored in:

```text
results/predictions.json
```

The four predictions test different parts of the model.

---

### P1. The hepatic shear set-point should be nearly invariant across mammals

The first prediction does not use transplant data.

Portal pressure remains within a relatively narrow physiological range across mammals despite large differences in body mass.

If the liver vascular network operates near the predicted optimum, the characteristic wall-shear set-point $\tau_0$ should therefore change very little with body size.

Using published allometric relationships:

$$\text{portal flow}\sim M^{0.77}$$

and:

$$\text{portocentral distance}\sim M^{0.10},$$

LiPNet predicts that the body-mass exponent of $\tau_0$ remains between:

$$-0.077\quad\text{and}\quad+0.071$$

over the full prespecified range of vascular maintenance exponents.

This range is close to zero.

The result therefore supports an approximately invariant hepatic shear set-point across mammals.

The vascular mass prediction provides a second, weaker test.

Across the same range of maintenance exponents, the predicted vascular mass exponent is:

$$0.89-1.08.$$

The reported range for hepatic blood volume is approximately:

$$0.80-0.92.$$

About 18% of the prespecified model range overlaps with the observed range.

No value of $b$ was selected to improve agreement with the observation.

#### Interpretation

The allometric analysis tests the physical structure of LiPNet independently of liver transplantation.

It asks whether the same network principle can remain physiologically reasonable from small mammals to humans.

---

### P2. Flow and pressure should separate when outflow resistance decreases

The central mechanistic prediction follows directly from:

$$z_P=z_Fz_R.$$

If venous outflow is enlarged:

$$z_R<1.$$

Portal pressure should then be lower than expected from portal flow alone.

Five published clinical groups reported both flow load and pressure load.

Four of these groups had reconstructed or functionally enlarged outflow:

- right lobe with middle hepatic vein reconstruction,
- outflow venoplasty,
- left lobe with middle and left hepatic vein trunk,
- splenectomy associated with increased haemodynamic compliance.

Their resistance loads were:

$$z_R=0.31,\;0.54,\;0.55,\;0.74.$$

All four values were below 1, as predicted.

The exact sign test gives:

$$p=0.06.$$

When the reported uncertainty in flow and pressure is propagated, each reconstructed group remains below:

$$z_R=1$$

with probability approximately 1.00.

These grafts had portal flow loads approximately three to four times the donor value, while reported adverse outcome rates remained between 0% and 10%.

### Interpretation

The result shows why portal flow and portal pressure should not be treated as interchangeable measures.

The same graft can have:

- very high flow,
- relatively moderate pressure,
- low effective resistance.

The most direct evidence for flow-pressure separation therefore comes from grafts in which both quantities were measured.

---

### P3. Within a cohort, outcome should increase with pressure load

If $z_P$ reflects the haemodynamic burden experienced by the graft, higher pressure load should be associated with worse outcomes within the same clinical cohort.

The complete dataset contains:

- 17 clinical series,
- 10 centres,
- 31 groups,
- 21 groups with pressure load,
- 14 groups with flow load.

A within-study slope requires at least two groups from the same study with the same outcome measured at different haemodynamic levels.

This information was available for:

- 5 series for $z_P$,
- 4 series for $z_F$.

Six additional series contained only one group and therefore inform the baseline level but not the within-study slope.

To compare pressure and flow using the same statistical structure, we fitted the same model to both loads and restricted the comparison to the primary outcome.

| Load | Groups / series | Events / patients | $\mu_\beta$ (95% CrI) | $P(\mu_\beta>0)$ |
|---|---:|---:|---:|---:|
| $z_P$, pressure | 7 / 3 | 76 / 573 | +1.56 (+0.35, +2.81) | 0.993 |
| $z_F$, flow | 8 / 4 | 48 / 321 | +0.86 (-0.58, +2.36) | 0.885 |

The interval for pressure load excludes zero.

The interval for flow load does not.

However, this does **not** prove that pressure is statistically superior to flow.

The probability that the pressure slope is larger than the flow slope is:

$$P(\mu_{z_P}>\mu_{z_F})=0.77.$$

The flow analysis also contains less information:

- 48 events,
- a median group size of 10 patients,
- two groups with zero events.

Wang 2014 contributes data to both analyses.

When the series defined by a pressure-gradient cut-off is added, the analysis contains:

- 13 groups,
- 6 strata,

with:

$$\mu_\beta=+1.96\quad(95\%\ \text{CrI }+1.13\text{ to }+2.81).$$

The main evidence that flow and pressure are not interchangeable remains Prediction 2, where both quantities are measured in the same graft groups.

---

### Hierarchical model

The main pressure analysis uses a hierarchical binomial model.

Within each study, pressure load is centred around the study mean:

```text
E_gs ~ Binomial(n_gs, p_gs)

logit(p_gs) =
    alpha_s
    + beta_s (zP_gs - mean(zP_s))

beta_s ~ Normal(mu_beta, tau_beta²)
```

where:

- $E_{gs}$ is the number of events in group $g$ of study $s$,
- $n_{gs}$ is the number of patients,
- $\alpha_s$ is the study-specific baseline,
- $\beta_s$ is the study-specific association with pressure load,
- $\mu_\beta$ is the average association across studies,
- $\tau_\beta$ describes between-study variation.

Centering $z_P$ within each study separates two questions:

1. What is the baseline risk in this centre?
2. How does risk change as pressure load changes within that centre?

This is important because baseline outcome rates vary strongly across centres.

---

### Primary clinical result

For the primary outcome of small-for-size syndrome or early dysfunction:

- 7 groups,
- 3 studies,
- 76 events,
- 573 patients.

The estimated mean slope was:

$$\mu_\beta=+1.55$$

with:

$$95\%\ \text{CrI}=+0.33\text{ to }+2.79.$$

This corresponds to an odds ratio of:

$$OR=4.70$$

per unit increase in $z_P$.

The posterior probability of a positive association was:

$$P(\mu_\beta>0)=0.990.$$

For the secondary outcome of mortality or graft loss:

$$\mu_\beta=+1.74\quad(95\%\ \text{CrI }+0.60\text{ to }+2.86).$$

When primary and secondary outcomes were combined:

$$\mu_\beta=+1.79\quad(95\%\ \text{CrI }+0.97\text{ to }+2.64).$$

Every within-study comparison moved in the predicted direction.

---

### Model diagnostics

The hierarchical model was sampled using NUTS with four chains.

Diagnostics were:

- maximum $\hat{R}$: 1.001,
- minimum effective sample size: 5114,
- divergences: 0,
- posterior predictive checks inside the expected interval: 7/7.

These results indicate stable sampling for the primary model.

---

### Sensitivity analyses

The pressure-load association was tested across 15 alternative specifications.

These included changes in:

- outcome definition,
- CVP assumptions,
- handling of overlapping Kyoto publications,
- prior on between-study variation,
- study versus centre stratification,
- imputation strategy.

Across these specifications:

$$\mu_\beta$$

remained between:

$$+0.92\quad\text{and}\quad+1.93.$$

The association remained positive in every specification.

The most conservative analysis combined publications from the same centre into a single stratum:

$$\mu_\beta=+0.92\quad(95\%\ \text{CrI }+0.20\text{ to }+1.65).$$

The prior on between-study variation has an effect, as expected with a small number of studies.

Using:

$$\tau\sim\text{HalfNormal}(0.25)$$

gives:

$$\mu_\beta=+1.60,\qquad P(\mu_\beta>0)=0.998.$$

Using:

$$\tau\sim\text{HalfNormal}(1.0)$$

gives:

$$\mu_\beta=+1.46\quad(-0.07,\,+2.84),$$

with:

$$P(\mu_\beta>0)=0.970.$$

When uncertainty in the calculated indices is propagated by drawing CVP values between 3 and 9 mmHg:

$$\mu_\beta=+1.74\quad(95\%\ \text{CrI }+0.89\text{ to }+2.63).$$

---

### Internal-external validation

Validation was performed by leaving out one centre at a time.

The slope was estimated using the remaining centres.

The held-out centre was then evaluated using a conditional likelihood that removes its intercept.

This means that no outcome information from the held-out centre is used to estimate its baseline risk.

The log-score improvement over a model with no association was:

- Cairo: +2.26,
- Fukuoka: +2.74,
- Kyoto: +8.40.

The predicted direction was correct in all three centres.

A recovery simulation using the real sample sizes and gradients showed the limits of the available evidence.

When the true effect was zero, the model did not incorrectly exclude zero.

For:

$$\mu=1,$$

the interval excluded zero in about 63% of simulations.

For:

$$\mu=2,$$

it excluded zero in all simulations.

The between-study variance $\tau_\beta$ cannot be estimated precisely with only five studies.

---

### P4. Clinical thresholds from different liver fields converge on the same scale

LiPNet predicts that haemodynamic thresholds developed independently in different clinical settings should become comparable when expressed as normalized pressure load.

This is what is observed.

#### Portal hypertension

Clinically significant portal hypertension is defined by:

$$HVPG\geq10\ \text{mmHg}.$$

On the LiPNet scale:

$$z_P=\frac{10}{5}=2.$$

---

#### Partial liver grafts

A commonly used graft threshold is:

$$PVP=15\ \text{mmHg}$$

with:

$$CVP=5\ \text{mmHg}.$$

Therefore:

$$z_P=\frac{15-5}{5}=2.$$

---

#### Hepatic resection

For a cirrhotic liver with a baseline pressure gradient of approximately 10 mmHg, reducing functional liver mass by resection increases the load carried by the remaining tissue.

The same normalized threshold appears around:

$$z_P=2.$$

---

#### Variceal bleeding

A pressure gradient associated with variceal bleeding is approximately:

$$HVPG\geq12\ \text{mmHg}.$$

This corresponds to:

$$z_P=\frac{12}{5}=2.4.$$

---

#### Interpretation

Thresholds developed independently in:

- cirrhosis,
- portal hypertension,
- partial liver transplantation,
- hepatic surgery,

map onto a similar normalized pressure-load scale.

The repeated appearance of:

$$z_P\approx2$$

suggests a common haemodynamic regime.

---

## What LiPNet is not

LiPNet is **not an individual risk-prediction model**.

Absolute outcome rates differ strongly between centres.

For example, at similar pressure gradients, different published series report very different rates of small-for-size syndrome or early dysfunction.

A single logistic curve fitted across all centres gives little useful information:

$$\text{slope}=0.11$$

with:

$$95\%\ \text{CI}=-0.25\text{ to }0.46.$$

The across-series rank correlation between $z_P$ and outcome is also weak (16 groups):

$$\rho=0.23,\qquad p=0.399.$$

For portal flow per gram (12 groups):

$$\rho=0.27,\qquad p=0.388.$$

This is expected from the model.

LiPNet does not explain why one centre has a higher baseline event rate than another.

The transferable quantity is the **change in risk associated with a change in haemodynamic load within a clinical setting**, not the absolute baseline risk.

---

## Whole-graft extension

The same framework can also be applied to whole grafts.

A whole graft keeps the donor lobular bed, but vascular resistance may change at the inlet or outlet anastomoses.

A simplified resistance model is:

```text
splanchnic inflow
        |
        v

     inlet
       |
   lobular bed
       |
     outlet
       |
       v

      CVP
```

with possible collateral pathways.

For the vascular anastomoses:

$$R_{\text{in}}\propto d_{\text{in}}^{-4}$$

and:

$$R_{\text{out}}\propto d_{\text{out}}^{-4}.$$

This produces asymmetric haemodynamic effects.

#### Inlet narrowing

If portal inflow is restricted, both:

$$z_F$$

and:

$$z_P$$

decrease together.

This behaviour is consistent with low-flow states and portal steal.

#### Outlet narrowing

If hepatic venous outflow is restricted, effective resistance increases.

Portal pressure therefore rises more strongly than portal flow.

This produces:

$$z_R>1$$

and is consistent with venous congestion or outflow obstruction.

Thus, inlet and outlet lesions move the graft in different directions in the:

$$(z_F,z_P)$$

plane.

---

## Calculator

A browser-based research calculator is available at:

**https://danpc11.github.io/LiPNet/**

The calculator runs locally in the browser.

Users can enter:

- current PVP,
- current CVP,
- expected pressure after a planned manoeuvre,
- a reference event rate or cohort-specific baseline.

The calculator returns:

- $z_F$,
- $z_P$,
- $z_R$,
- the estimated odds ratio associated with a pressure change,
- an absolute risk estimate anchored to the user-provided baseline.

The slope and uncertainty come from the hierarchical model described above.

If a user provides a reference probability:

$$p_{\mathrm{ref}}$$

at pressure load:

$$z_{\mathrm{ref}},$$

the model uses:

$$\operatorname{logit}p(z)=\operatorname{logit}p_{\mathrm{ref}}+\mu_\beta(z-z_{\mathrm{ref}}).$$

The resulting curve therefore passes through the supplied clinical anchor.

---

### Recalibration

A centre can estimate its own baseline using:

```python
indices.baseline_from_rate(rate, mean_zP, beta)
```

which computes:

```text
alpha = logit(rate) - beta * mean(zP)
```

For local outcome data, the intercept can be recalibrated using:

```python
indices.recalibrate(outcomes, zP, beta)
```

while keeping the slope fixed.

This corresponds to recalibration-in-the-large, a standard first step in prediction-model updating.

---

## Repository structure

```text
src/model/
    Network model:
    - 2D and 3D vascular graphs
    - optimization
    - closed-form scaling
    - allometric closure

src/predictions.py
    Tests predictions P1-P4
    Output:
    results/predictions.json

src/indices.py
    Calculates:
    - zF
    - zP
    - zR
    - stratified models
    - hierarchical models
    - internal-external validation
    - recalibration

src/plots.py
    Generates 18 stand-alone figures

src/build_app.py
    Rebuilds the browser calculator

data/series.tsv
    31 groups from 17 published clinical series

sfss_calculator.html
index.html
tests/
run_all.sh
```

---

## Reproduce the analysis

Install the dependencies:

```bash
pip install -r requirements.txt
```

Run the complete workflow:

```bash
./run_all.sh
```

This reproduces:

- haemodynamic indices,
- predictions,
- statistical analyses,
- figures,
- calculator.

Typical runtime is approximately three minutes.

To run only the four main predictions:

```bash
python src/predictions.py
```

To generate the summary prediction figure:

```bash
python src/plots.py predictions
```

The model campaigns used to generate the cached tables are stored in:

```text
src/model/
```

They include:

- optimized networks from 37 to 1478 lobules,
- 2D and 3D networks,
- exhaustive optimization on seven lobules,
- allometric closure.

The corresponding commands are documented in:

```text
run_all.sh
```

---

## Data

The main clinical dataset is:

```text
data/series.tsv
```

It contains raw published values.

The normalized indices are not stored as fixed results.

They are recalculated by:

```text
src/indices.py
```

Each variable includes provenance information.

The `*_basis` columns identify whether each value was:

- reported,
- derived,
- imputed,
- defined by a threshold,
- not applicable.

Structural variables include:

- `centre_id`,
- `cohort_id`,
- `recruitment_start`,
- `recruitment_end`,
- `overlap_set`,
- `outflow`,
- `gradient_timing`,
- `gradient_sd`,
- `outcome_definition`,
- `outcome_horizon`,
- `modulation_strategy`.

These fields are used to distinguish publications, cohorts, centres and overlapping patient populations.

Groups defined only by a clinical cut-off are excluded from the primary analysis.

Second partitions of already-counted cohorts are also excluded from the primary analysis.

They are reported separately in:

```text
results/sensitivity.tsv
```

---

### Included clinical series

The dataset includes the following published series:

- Troisi 2003
- Troisi 2005
- Ou 2010
- Vasavada 2014
- Alim 2016
- Chan 2011
- Yagi 2005
- Yagi 2006
- Wang 2014
- Osman 2017
- Ogura 2010
- Yamada 2008
- Uemura 2016
- Yao 2018
- Kanetkar 2017
- Ishizaki 2012
- Botha 2010

---

## Limitations

The main clinical analysis uses aggregated group-level data.

Prediction 3 is based on five series from four centres, with 104 clinical events.

Only three series report the primary outcome.

Several Kyoto publications overlap in recruitment period.

These limitations reduce the amount of independent information available for estimating between-study heterogeneity.

CVP is not available in all studies.

When CVP is assumed as a constant within a study, it shifts all groups from that study by the same amount on the $z_P$ axis.

This does not change the within-study slope because the study-specific intercept absorbs the shift.

It can, however, change the absolute position of the groups on the normalized pressure scale.

Because the available data are aggregated, the analysis can estimate:

- likelihood,
- deviance,
- observed-to-expected ratios,
- study-level slopes.

It cannot provide reliable individual-level measures such as:

- c-statistic,
- individual calibration curves,
- Brier score,
- decision-curve analysis.

These limitations apply mainly to Prediction 3.

Predictions 1, 2 and 4 do not depend on the hierarchical outcome model.

The calculator is therefore a **research tool**.

It is not a medical device and should not be used as a stand-alone tool for clinical decision-making.

---

## Main interpretation

LiPNet proposes a simple physical interpretation of hepatic haemodynamics.

Portal flow alone does not define the haemodynamic load of a liver graft.

The same portal flow can produce different portal pressures depending on vascular resistance.

The central relationship is:

$$\boxed{z_P = z_F \, z_R}$$

where:

- $z_F$ describes flow load,
- $z_R$ describes resistance load,
- $z_P$ describes the resulting pressure load.

This framework explains how a partial graft can carry three to four times the donor portal flow while maintaining a much smaller increase in portal pressure when effective vascular resistance is reduced.

It also provides a common normalized language for:

- portal hyperperfusion,
- venous outflow reconstruction,
- small-for-size graft physiology,
- portal hypertension,
- hepatic resection,
- vascular obstruction.

LiPNet is therefore intended as a mechanistic framework rather than an individual risk score.

---

## Publication

The article associated with this repository is currently in preparation:

**Vascular resistance modulates portal flow-pressure decoupling in partial liver grafts**

Until a preprint is available, please cite the archived software release.

---

## Citation

The software is archived on Zenodo:

**DOI: 10.5281/zenodo.22861276**

A machine-readable citation is also available in:

```text
CITATION.cff
```

---

## Licence

Source code is distributed under:

**PolyForm Noncommercial 1.0.0**

https://polyformproject.org/licenses/noncommercial/1.0.0

Clinical values stored in `data/` remain the intellectual property of the authors of the original publications.
