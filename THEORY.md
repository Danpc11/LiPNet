# LiPNet: the model and its predictions

This document holds the physics and the four predictions. The [README](README.md) covers what the repository contains and how to use it.

---

## Contents

- [The question](#the-question)
- [The model](#the-model)
- [Partial grafts](#partial-grafts)
- [Three normalized haemodynamic loads](#three-normalized-haemodynamic-loads)
  - [Flow load](#flow-load)
  - [Pressure load](#pressure-load)
  - [Resistance load](#resistance-load)
  - [Why resistance matters](#why-resistance-matters)
- [What the model includes](#what-the-model-includes)
- [Four predictions](#four-predictions)
  - [P1. The scale of the design should be nearly invariant across mammals](#p1-the-scale-of-the-design-should-be-nearly-invariant-across-mammals)
  - [P2. Flow load and pressure load are not interchangeable](#p2-flow-load-and-pressure-load-are-not-interchangeable)
  - [P3. Within a cohort, outcome should increase with pressure load](#p3-within-a-cohort-outcome-should-increase-with-pressure-load)
  - [Hierarchical model](#hierarchical-model)
  - [Primary clinical result](#primary-clinical-result)
  - [Model diagnostics](#model-diagnostics)
  - [Sensitivity analyses](#sensitivity-analyses)
  - [Internal-external validation](#internal-external-validation)
  - [P4. Exploratory comparison of haemodynamic thresholds across liver diseases and surgical settings](#p4-exploratory-comparison-of-haemodynamic-thresholds-across-liver-diseases-and-surgical-settings)
- [What LiPNet is not](#what-lipnet-is-not)
- [Whole-graft extension](#whole-graft-extension)

---

## The question

Partial liver grafts receive portal blood through less tissue than a whole liver, so portal flow per unit of mass can be several times the donor's. High flow, however, does not always produce equally high pressure: some grafts carry it with a relatively low portal pressure, above all when the venous outflow is wide or reconstructed.

> **What determines how portal flow is converted into portal pressure?**

The three loads and the identity $z_P = z_F \, z_R$ are stated on the [front page](README.md#what-this-is). What follows is where they come from and what they predict.

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

in the limit of many lobules (for $b<1$; at the size of a human liver the local exponent of $D^*$ with $N$ is 1.2 to 3.0 across $b$, [docs/LIVER.md](docs/LIVER.md)), together with the scale of the design:

$$\tau_0=\sqrt{D^*/C},$$

the dissipation per unit of maintenance cost, which is the same on every vessel of the optimum. The wall shear stress of a vessel of vascular mass $m_e$ is

$$\tau_e=\tau_0\,m_e^{(b-1)/2},$$

so $\tau_0$ is a wall shear stress only at $b=1$ (Murray's law); for $b<1$ the shear rises from the trunks to the sinusoids.

and the pressure drop across a lobule:

$$\Delta P=fR.$$

Here:

- $D^*$ is the minimum dissipation of the network,
- $F$ is total flow,
- $C$ represents vascular maintenance cost,
- $N$ is the number of lobules,
- $L$ is a characteristic transport length,
- $\tau_0$ is the scale of the design,
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

Vascular resistance per unit mass is:

$$R=\frac{\Delta P}{Q}\quad\text{(per 100 g)}.$$

The normalized resistance load is defined as the graft resistance relative to the donor:

$$z_R=\frac{R_{\text{graft}}}{R_{\text{donor}}}.$$

$R_{\text{donor}}$ uses the reference gradient of 5 mmHg and the donor reference flow. Dividing the two other loads:

$$\frac{z_P}{z_F}=\frac{\Delta P_{\text{graft}}/\Delta P_{\text{donor}}}{Q_{\text{graft}}/Q_{\text{donor}}}=\frac{R_{\text{graft}}}{R_{\text{donor}}}=z_R.$$

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

But if the effective resistance of the graft is half the donor's:

$$z_R=0.5,$$

then:

$$z_P=3\times0.5=1.5.$$

The graft therefore carries very high flow without developing the pressure expected from flow alone.

This is **portal flow-pressure decoupling** in LiPNet: the two loads separate by exactly the resistance load.

The network itself keeps a partial graft close to the donor's resistance. A graft is the sub-network below a first-order branch, and loses only the trunks it left behind: $z_R$ is 0.93 to 1.00 on the closed-form network and 0.90 to 0.96 on explicit three-dimensional graphs. The hepatic artery joins the portal blood at the sinusoid, so a load read from portal flow alone overstates what the sinusoids and hepatic veins carry: with the arterial flow measured in donor and graft, $z_R$ is 0.82 to 0.90. The derivations are in [docs/LIVER.md](docs/LIVER.md) and [docs/HUMAN.md](docs/HUMAN.md). The model therefore expects $z_R$ between 0.8 and 1, and does not name the element that would lower it further.

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

The four predictions test different parts of the model, and they are not evidence of the same kind. Three levels are kept apart throughout:

| Level | What it is | Where |
|---|---|---|
| Identity | $z_P=z_Fz_R$, exact because $z_R$ is defined as $z_P/z_F$. It cannot support a mechanism and is not tested. | [Three normalized haemodynamic loads](#three-normalized-haemodynamic-loads) |
| Prediction of the network | A number the model gives without the pressure of any recipient: the radii of a human portal tree, and $z_R$ of a partial graft from its geometry and its two inflows. Each can be wrong, and the second is. | P2, [docs/HUMAN.md](docs/HUMAN.md), [docs/LIVER.md](docs/LIVER.md) |
| Association in the published series | What the cohorts report: $z_R$ by surgical configuration, and outcome against pressure load between groups of one cohort. Observational and aggregated. | P2, P3 |

P1 is a check of consistency and P4 an exploratory comparison of scales whose agreement is arithmetic; neither is a validation.

---

### P1. The scale of the design should be nearly invariant across mammals

The first prediction does not use transplant data.

Portal pressure remains within a relatively narrow physiological range across mammals despite large differences in body mass.

If the liver vascular network operates near the predicted optimum, the scale of the design $\tau_0$ should therefore change very little with body size.

Using published allometric relationships:

$$\text{portal flow}\sim M^{0.77}$$

and:

$$\text{portocentral distance}\sim M^{0.10},$$

LiPNet predicts that the body-mass exponent of $\tau_0$ remains between:

$$-0.077\quad\text{and}\quad+0.071$$

over the full prespecified range of vascular maintenance exponents.

This range is close to zero.

The scale of the design is therefore approximately invariant across mammals. The wall shear stress it implies depends on the vessel: its body-mass exponent is $-0.077$ to $-0.003$ in the sinusoid and $-0.30$ to $-0.03$ in the trunk ([docs/LIVER.md](docs/LIVER.md)).

The vascular mass prediction provides a second, weaker test.

Across the same range of maintenance exponents, the predicted vascular mass exponent is:

$$0.89-1.08.$$

The reported range for hepatic blood volume is approximately:

$$0.80-0.92.$$

About 18% of the prespecified model range overlaps with the observed range.

No value of $b$ was selected to improve agreement with the observation.

#### What P1 shows

P1 is a check of consistency, not an independent validation. The exponent of $\tau_0$ is not compared with any measurement: it follows from the inputs of the closure (an invariant portal pressure and the measured scalings of flow and lobule size) and from the exponent of the optimum, and the check is that the closure does not force the design scale to change strongly with body mass. The one comparison with data that the closure did not use is the exponent of the vascular mass against the scaling of hepatic blood volume, and it is met in 5 of the 28 combinations of inputs.

The independent test of the network is anatomical and is in [docs/HUMAN.md](docs/HUMAN.md): with no fitted scale the optimum reproduces the eleven measured radii of a human portal cast within 10% at $b=0.75$, and neither end of the prespecified range of $b$ does.

---

### P2. Flow load and pressure load are not interchangeable

The central prediction follows directly from:

$$z_P=z_Fz_R.$$

Flow and pressure carry the same information only if $z_R=1$. The network keeps a partial graft between 0.8 and 1 (see [Why resistance matters](#why-resistance-matters)), so the model expects the pressure load to fall short of the flow load, by at most a fifth.

Six published clinical groups from five cohorts reported both flow load and pressure load (two of the groups are the two arms of Wang 2014). Their resistance loads were:

$$z_R=0.31,\;0.54,\;0.55,\;0.65,\;0.72,\;0.74.$$

**Description.** All are below 1, at flow loads of 3.3 to 4.1. The unit of the test is the cohort: five of five are below 1, and the exact sign test in the direction the network gives is

$$p=0.031\ \text{one-sided},\qquad p=0.063\ \text{two-sided}.$$

The direction is the one the model gives, since a graft sheds trunks of the donor tree and cannot gain resistance from that. The test was not registered before the extraction, so both values are reported and the two-sided one is the one to read without that premise.

The six values are six descriptions and five independent units. Dropping either arm of Wang 2014 leaves five of five, and counting by centre gives the same five. With five units the exact one-sided sign test cannot give less than 0.031, and leaving any one cohort out gives 0.0625. The test says that five of five agree in sign and nothing about the size of the departure.

**Uncertainty of each value, by source.** Two sources are kept apart. The first is the sampling error of the two group means. It is computed from the dispersions each paper reports and, where none is reported, from the largest reported one, at zero correlation between the two means; where a dispersion is not reported it is then swept (coefficient of variation of the flow 0.25 to 0.70, standard deviation of the gradient 2.5 to 6 mmHg), together with a correlation of $-0.5$, 0 and $+0.5$ between the sampling errors of the two means, which no paper reports. The second is the references, which are not data of the group: the donor flow (80 to 100 mL/min per 100 g where 90 is assumed, and the same relative range around a value measured in the series), the normal gradient (4 to 6 mmHg) and the central venous pressure where it is assumed (3 to 9 mmHg), each drawn uniformly over its range. Intervals are 95%.

| Group | Dispersion reported | $z_R$ | Sampling | Sampling, widest of the sweep | References | Both | $P(z_R<1)$ | Variance from references |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Chan 2011 | both | 0.31 | 0.25 to 0.37 | 0.24 to 0.38 | 0.24 to 0.40 | 0.23 to 0.42 | 1.00 | 67% |
| Wang 2014, splenectomy | none | 0.54 | 0.48 to 0.60 | 0.45 to 0.64 | 0.32 to 0.73 | 0.31 to 0.74 | 1.00 | 94% |
| Yagi 2006 | flow only | 0.55 | 0.44 to 0.69 | 0.40 to 0.76 | 0.43 to 0.71 | 0.40 to 0.78 | 1.00 | 57% |
| Wu 2011 | both | 0.65 | 0.57 to 0.75 | 0.55 to 0.77 | 0.52 to 0.84 | 0.50 to 0.88 | 1.00 | 78% |
| Wang 2014, no splenectomy | none | 0.72 | 0.64 to 0.80 | 0.60 to 0.85 | 0.45 to 0.95 | 0.45 to 0.97 | 0.98 | 92% |
| Ishizaki 2012 | both | 0.74 | 0.60 to 0.92 | 0.57 to 0.96 | 0.59 to 0.96 | 0.54 to 1.04 | 0.95 | 59% |

Sampling error alone does not bring any group to 1: the upper limit stays below 1 in every one of the 141 combinations of the sweep, and reaches 0.96 at most (Ishizaki 2012). The references are the larger source, 57 to 94% of the variance of $\log z_R$, and most where the central venous pressure is assumed (Wang 2014). With both sources the interval of Ishizaki 2012 reaches 1.04 and that of the arm of Wang 2014 without splenectomy 0.97: five of the six groups, and four of the five cohorts, stay below 1. The uniform ranges are a choice, not a measurement, and the references are shared between groups, so these intervals describe one group at a time; what the references do to all groups at once is in the paragraph on references below.

**Prediction of the network against observation.** The network predicts $z_R$ of each group from the geometry of a lobar graft, the arterial share of the donor and the flow load of the group. No pressure of a recipient enters.

| Group | Configuration | $z_F$ | $z_R$ predicted | $z_R$ observed (interval) | Share of the departure from 1 predicted |
|---|---|---|---|---|---|
| Chan 2011 | outflow reconstruction | 3.93 | 0.81 to 0.89 | 0.31 (0.25 to 0.37) | 0.16 to 0.27 |
| Wang 2014, splenectomy | inflow modulation | 4.12 | 0.81 to 0.89 | 0.54 (0.48 to 0.60) | 0.24 to 0.41 |
| Yagi 2006 | outflow and pressure control | 3.28 | 0.81 to 0.90 | 0.55 (0.44 to 0.69) | 0.22 to 0.42 |
| Wu 2011 | outflow reconstruction | 3.25 | 0.81 to 0.90 | 0.65 (0.57 to 0.75) | 0.29 to 0.54 |
| Wang 2014, no splenectomy | none | 3.74 | 0.81 to 0.89 | 0.72 (0.64 to 0.80) | 0.37 to 0.67 |
| Ishizaki 2012 | outflow reconstruction | 3.34 | 0.81 to 0.90 | 0.74 (0.60 to 0.92) | 0.39 to 0.73 |

The mean prediction is 0.85 against a mean observation of 0.58 (mean error $-0.27$, root-mean-square error 0.31). The prediction is too high in the six groups. The whole sampling interval of the observation lies below it in five (four with the widest dispersions of the sweep), and the interval that also carries the references in three: Chan 2011, Yagi 2006 and the splenectomy arm of Wang 2014. The model has the direction of the departure from 1 and under-predicts its size: it is falsified as a complete account of the resistance of a graft. For Wu 2011, which measured the arterial flow of donor and graft, the prediction narrows to 0.82 to 0.86 against 0.65.

**What the rigid network leaves out.** The network is rigid, and the venous bed of a liver is not. Sinusoids carry 3 to 15% of the pressure drop in the optimum and the outlet anastomosis 5% of the resistance, so neither can supply the further 32% by which the resistance would have to fall. Two published relations for the passive distensibility of the normal liver, resistance as the inverse cube of the distending pressure and a blood volume that doubles for 8 mmHg, give nearly the same answer with nothing fitted: a fully distensible bed at these flow loads has $z_R$ of 0.39 to 0.58. The same relations give pig livers, which they were not derived from, a resistance of 0.48 to 0.57 at 2.7 times the flow, against 0.54 and 0.59 measured. Of the six groups, three lie between the bed and the rigid network, two at the bed and Chan 2011 below both. The missing resistance is therefore of the size passive distension of the venous trees can supply. That is a bound, not an identification: both relations come from cats over the normal range of pressure and are extrapolated, and after cold ischaemia the same pig livers had twice the resistance ([docs/LIVER.md](docs/LIVER.md), section 7).

**Configuration.** The table classifies each group by mechanism (column `configuration`): outflow reconstruction, inflow modulation, portosystemic shunt, combined, none. A splenectomy is an inflow manoeuvre and is not counted as an outflow reconstruction. Among the six groups, three had an outflow reconstruction and no inflow manoeuvre (0.31, 0.65, 0.74), one an outflow venoplasty with pressure control (0.55), one a splenectomy (0.54) and one no manoeuvre (0.72). The group with no manoeuvre is inside the range of the three with outflow reconstruction, so these data do not show that reconstruction lowers $z_R$. Within the model, widening the outlet anastomosis lowers $z_R$ to 0.95 at most.

The only contrast inside one cohort is Wang 2014: $z_R$ was 0.54 with splenectomy and 0.72 without, a ratio of 0.75 (0.64 to 0.88), at a higher flow load (4.12 against 3.74). The arms were not randomised and the manoeuvre was chosen for the recipients who needed it.

Three questions are therefore answered differently. Is $z_R$ below 1 in the published cohorts? Yes, in all five. Is it lower with outflow reconstruction than without? Not shown. Is a difference caused by the manoeuvre? Not addressable with these data.

**References.** Two references enter $z_R$: the donor's gradient, which no series measures and is taken as 5 mmHg, and the donor's portal flow, measured in four of the series (81, 85, 91 and 118 mL/min per 100 g) and taken as 90 elsewhere. With a donor flow of 80 to 100, a normal gradient of 4 to 6 mmHg and a central venous pressure of 3 to 9 mmHg where one is assumed, varied together (60 combinations), every group stays below 1 in 54; the six exceptions all have a normal gradient of 4 mmHg. Being below the network is a weaker statement: every group is below 0.81 in 34 of the 60, in none of the 20 with a normal gradient of 4 mmHg, in 14 of the 20 at 5 mmHg and in all 20 at 6 mmHg. Chan 2011 and Yagi 2006 are below the network in every combination; Ishizaki 2012 and Wu 2011 in 60% and 67%.

#### What P2 shows

Portal flow and portal pressure should not be treated as interchangeable measures. The same graft can have:

- very high flow,
- relatively moderate pressure,
- low effective resistance.

The resistance load is what separates them. At the reference values it is lower in every published group than the network accounts for. That it is below 1 survives the sampling error in every group and the uncertainty of the references in five of six. That it is below the network survives both in three of six; for the three groups closest to the network it depends on a donor gradient of at least 5 mmHg.

---

### P3. Within a cohort, outcome should increase with pressure load

If $z_P$ reflects the haemodynamic burden experienced by the graft, higher pressure load should be associated with worse outcomes within the same clinical cohort.

The complete dataset contains:

- 18 clinical series,
- 11 centres,
- 32 groups,
- 22 groups with pressure load,
- 15 groups with flow load.

A within-study slope requires at least two groups from the same study with the same outcome measured at different haemodynamic levels.

This information was available for:

- 5 series for $z_P$,
- 4 series for $z_F$.

Six additional series contained only one group and therefore inform the baseline level but not the within-study slope.

As an exploratory analysis, the same model was fitted to each load on the series that provide a within-study contrast for it, both restricted to the primary outcome. The row for pressure is the primary fit itself, not a second run of the sampler.

| Load | Groups / series | Events / patients | $\mu_\beta$ (95% CrI) | $P(\mu_\beta>0)$ |
|---|---:|---:|---:|---:|
| $z_P$, pressure | 7 / 3 | 76 / 573 | +1.55 (+0.33, +2.79) | 0.990 |
| $z_F$, flow | 8 / 4 | 48 / 321 | +0.84 (-0.59, +2.36) | 0.872 |

The interval for pressure load excludes zero and the interval for flow load does not. This does **not** show that pressure predicts better than flow, and the two slopes are not compared: no probability that one exceeds the other is reported. The two rows are different cohorts except for Wang 2014, which enters both, so the two posteriors are neither the same recipients nor independent. A comparison would need a joint model of cohorts that report both loads and an outcome, and there is one such cohort, described below.

The flow analysis also contains less information:

- 48 events,
- a median group size of 10 patients,
- 2 groups with zero events.

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

This is an association between groups of one cohort, an ecological one. It is not the effect of pressure on a patient: the groups are means, they were often defined by the pressure itself or by a manoeuvre chosen because of it, and the manoeuvre changes both the pressure and the risk. [docs/IPD.md](docs/IPD.md) integrates over the dispersion inside a group and does not remove confounding by indication.

The mean slope is the average over centres. The slope expected in a centre not in the data is wider: $-0.02$ to $+3.09$ for the primary outcome (probability of a positive slope 0.97) and $+0.60$ to $+2.94$ for all outcomes (0.99).

One cohort reports both loads and an outcome in the same recipients (Wang 2014). The arm with the worse outcome (19.7% against 9.7%, Fisher $p=0.023$) had the higher pressure load (2.68 against 2.22) and the lower flow load (3.74 against 4.12). The hierarchical model fitted to the flow load uses mostly other cohorts and is not a comparison of how well the two loads predict.

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

The result depends on the prior on the between-study standard deviation $\tau_\beta$, which three series cannot estimate. For the primary outcome the mean slope is $+1.60$ ($+0.48$ to $+2.76$) under HalfNormal(0.25), $+1.55$ ($+0.33$ to $+2.79$) under HalfNormal(0.5), the prespecified one, and $+1.46$ ($-0.07$ to $+2.84$) under HalfNormal(1.0): with the widest prior the interval reaches zero, with a probability of a positive slope of 0.97. Convergence of the sampler says nothing about this.

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

$$+0.92\quad\text{and}\quad+1.95.$$

The association remained positive in every specification.

The most conservative analysis combined publications from the same centre into a single stratum:

$$\mu_\beta=+0.92\quad(95\%\ \text{CrI }+0.20\text{ to }+1.67).$$

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

A recovery simulation using the real sample sizes, gradients and strata showed the limits of the available evidence. Outcomes were simulated with a known mean slope and refitted, 20 times per value.

| True $\mu_\beta$ | Mean estimate | Bias | Coverage of the 95% CrI | Interval excludes zero |
|---:|---:|---:|---:|---:|
| 0 | -0.07 | -0.07 | 1.00 (0.83 to 1.00) | 0.00 (0.00 to 0.17) |
| 1 | +0.98 | -0.02 | 1.00 (0.83 to 1.00) | 0.75 (0.51 to 0.91) |
| 2 | +1.68 | -0.32 | 0.95 (0.75 to 1.00) | 1.00 (0.83 to 1.00) |

The proportions are over 20 refits and carry their exact (Clopper-Pearson) 95% interval. When the true effect was zero the interval excluded zero in 0 of 20; for a slope of 1 in 15 of 20; for a slope of 2 in 20 of 20.

The simulation is simpler than the data in four ways, and its answer holds within them:

- each stratum keeps its own level, the logit of its pooled observed rate, and the predictor is taken as exact;
- the heterogeneity of the slope is set beforehand, $\tau_\beta=0.3$, and not varied;
- outcomes are binomial with no overdispersion within a group, which is the assumption of the model itself, so the simulation tests the sampler and the priors on this design and not the model;
- 20 refits per value: a coverage of 1.00 is compatible with a true coverage of 0.83.

The between-study variance $\tau_\beta$ cannot be estimated precisely with only five studies.

---

### P4. Exploratory comparison of haemodynamic thresholds across liver diseases and surgical settings

Thresholds used in different clinical settings are placed on the normalized pressure-load scale. This is a comparison of scales. It is exploratory, it is not a statistical result, and it is not a prediction that was tested.

#### The agreement at 2 is arithmetic

The thresholds already agree before any normalization: each is a gradient of about 10 mmHg.

| Setting | Threshold | Gradient | $z_P$ at 5 mmHg |
|---|---|---:|---:|
| Cirrhosis, clinically significant portal hypertension | $HVPG\geq10$ mmHg | 10 mmHg | 2.0 |
| Partial graft | $PVP=15$ mmHg at $CVP=5$ mmHg | 10 mmHg | 2.0 |
| Resection in cirrhosis, limit | $HVPG=10$ mmHg, the threshold of the first row used as a contraindication | 10 mmHg | 2.0 |
| Variceal bleeding | $HVPG\geq12$ mmHg | 12 mmHg | 2.4 |

Dividing 10 mmHg by the normal gradient of 5 mmHg gives 2 in every row:

$$z_P=\frac{10}{5}=2.$$

So the value 2 is a consequence of the reference chosen for the scale, not something the model predicts and not an independent agreement between fields. With a reference of 4 mmHg the same thresholds sit at 2.5, and with 6 mmHg at 1.67. The third row is not a third threshold: it is the first one applied to surgery.

#### The measurements differ

The thresholds are not the same measurement. The hepatic venous pressure gradient (HVPG) is wedged minus free hepatic venous pressure: it reads the sinusoidal gradient and equals the portocaval gradient only when the resistance is sinusoidal and the hepatic veins and the cava add nothing; it underestimates portal pressure in presinusoidal portal hypertension. The gradient of a graft is portal minus central venous pressure and includes the outflow. No equivalence between these thresholds has been validated.

#### What the scale adds, and what remains untested

The scale gives the number a reading: $z_P=2$ is twice the pressure drop that one lobule of a donor liver carries. The model adds one statement of its own, which is not tested here: a resection that removes a fraction $x$ of a liver with baseline gradient $G$ leaves the remnant at

$$z_P=\frac{G}{5\,(1-x)},$$

if the portal flow is unchanged and the lobules are in parallel (Figure 4B of the article).

In the article this comparison is Supplementary Figure S7 and a sentence of the discussion; it is not a result.

#### What P4 shows

Thresholds developed in cirrhosis, in partial liver transplantation and in hepatic surgery can be written on one scale, and on that scale they are close because they were close in mmHg. That is all. P4 does not show that the three settings share a physiological regime or a common threshold, and it does not validate $z_P$. A test would need the same measurement, with outcomes, in more than one setting.

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

---

Back to the [README](README.md) for the repository, the data and the calculator.
