# Patients inside groups, and curves inside summaries

The clinical analysis in this repository is built on group-level data: one row per reported
group, with a mean pressure load, a denominator and an event count. The two weaknesses a
reviewer will name are both about that, and both have a method attached rather than an
apology.

1. The regression uses the **mean** load of a group where the quantity of interest is the
   load of a **patient**. That is the ecological fallacy, and in a logistic model it has a
   known direction.
2. The series report **summaries** of survival where they published **curves**. A
   Kaplan-Meier step function plus its numbers at risk is close to a complete record, and
   the individual event times can be solved back out of it.

Neither module in `src/ipd/` creates information. The first integrates over a dispersion
that six of the groups already report; the second inverts a function the original authors
published. Synthetic patients are not generated anywhere, and nothing here is presented as
measured data.

---

## 1. The slope for a patient, not for a group mean

### What the series actually report

| series | group | gradient (mmHg) | dispersion as printed | within-group SD |
|---|---|---|---|---|
| Ogura 2010 | final PVP < 15 (n = 86) | 6.2 | ± 0.3 | **2.78** |
| Ogura 2010 | final PVP >= 15 (n = 43) | 9.9 | ± 0.4 | **2.62** |
| Uemura 2016 | GRWR < 0.7 | 5.1 | SD 4.1 | 4.1 |
| Uemura 2016 | GRWR 0.7-0.8 | 4.75 | SD 3.8 | 3.8 |
| Uemura 2016 | GRWR >= 0.8 | 4.73 | SD 3.1 | 3.1 |
| Ishizaki 2012 | left lobe, no modulation | 12.4 | SD 4.4 | 4.4 |
| Botha 2010 | left lobe with HPCS | 6.2 | not printed | **3.91** |
| Chan 2011 | right lobe with MHV (n = 46) | 6 | range −1 to 11 | **2.71** |

Two of those corrections matter and one of them was nearly missed.

**Ogura's dispersions are standard errors.** Its Results section reports the same cohort's
MELD as 23.2 ± 1.1 (n = 86) and 22.9 ± 1.7 (n = 43), while its own Table 1 gives the MELD
of that cohort as 22.9 ± 10.6. Multiplying the first two by the square root of their sample
sizes recovers 10.2 and 11.1. The smaller figures are therefore standard errors of the
mean, whatever the methods section says, and the same conversion applied to the gradient
turns 0.3 and 0.4 mmHg into standard deviations of 2.78 and 2.62. `ec.sem_to_sd` does the
conversion and the arithmetic is checked in `tests/test_ipd.py`.

**Two dispersions are computed rather than read off.** Botha publishes one row per
recipient, so the standard deviation comes straight from the ten measured gradients: 3.91
mmHg. Chan publishes a range for all 46 recipients, which the estimator of Wan et al.
(*BMC Med Res Methodol* 2014) turns into 2.71 mmHg.

With all four corrections, the eight dispersions run from 2.6 to 4.4 mmHg and agree within
a factor of 1.7. That is the result: **the within-group spread of the portocaval gradient is
about 3 to 4 mmHg in every series that measured it, which is as wide as the entire
between-group range of the meta-analysis.** A single pooled value of 3.4 mmHg
(0.69 in z_P units) is therefore defensible, where before the correction the dispersions
looked bimodal and no pooling would have been.

### The model

Patient *i* of group *g* carries an unobserved load `x_i ~ Normal(mu_g, sigma_g^2)`, where
`mu_g` is the reported mean and `sigma_g` the within-group SD, both in units of
`z_P = (PVP - CVP) / 5 mmHg`. The group's event probability is the average of the
individual logistic curve over that distribution,

    pbar_g = Integral logistic(alpha_s + beta_s (x - xbar_s)) Normal(x; mu_g, sigma_g^2) dx

with `alpha_s` the intercept of stratum *s*, `beta_s` its slope and `xbar_s` the stratum
mean; `events_g ~ Binomial(n_g, pbar_g)`. The integral is evaluated by 21-point
Gauss-Hermite quadrature. At `sigma = 0` it collapses to the group-mean model of
`src/bayes.py`, and `beta` then estimates the slope for a patient rather than for a group
mean.

### The sampler is checked before it is used

`src/ipd/ecological.py` samples with an adaptive random-walk Metropolis in NumPy, so the
module needs no new dependency. Run at `sigma = 0` on the primary outcome it must return
the published NUTS fit, and it does:

| | groups | patients | events | mu_beta | 95% CrI | tau |
|---|---|---|---|---|---|---|
| NUTS, `src/bayes.py` | 7 | 573 | 76 | +1.547 | +0.33 to +2.79 | 0.365 |
| this sampler, sigma = 0 | 7 | 573 | 76 | +1.56 | +0.31 to +2.80 | 0.39 |

R-hat 1.005. That agreement is a test in `tests/test_ipd.py`, not a remark.

### Result

| assumed within-group SD | 0 | 1 mmHg | 2 mmHg | 3 mmHg | 4 mmHg | 5 mmHg | 6 mmHg |
|---|---|---|---|---|---|---|---|
| mu_beta | +1.56 | +1.58 | +1.66 | +1.74 | +1.68 | +1.59 | +1.39 |
| odds ratio per unit z_P | 4.8 | 4.8 | 5.3 | 5.7 | 5.3 | 4.9 | 4.0 |

Using each group's own reported SD, and the pooled median of 0.69 z_P units (3.45 mmHg)
where it is missing, gives mu_beta = +1.71 (+0.34 to +3.28).

Two things follow, and the second matters more than the first.

**The association does not depend on the ecological assumption.** The slope stays positive
with probability above 0.99 at every spread tested, including spreads larger than any the
series report.

**Over the plausible range the bias runs one way: against the published estimate.** From 1
to 4 mmHg of within-group spread, the patient-level slope is between 1% and 11% *larger*
than the group-mean slope. Dispersion inside a group flattens a logistic relation, so
fitting group means understates what a patient experiences. The published +1.55 is a
conservative reading, not an optimistic one. The curve turns over beyond 5 mmHg, where the
groups become statistically indistinguishable and the prior takes over; that is the
sampler reporting that the data have run out, which is also worth showing.

```bash
python src/ipd/run_ipd.py ecological
python src/ipd/run_ipd.py figures     # results/ipd_ecological.pdf
```

---

## 2. Individual survival data from a published curve

A Kaplan-Meier step function is not a summary. Every step is an event; the drop in the
printed numbers at risk between steps fixes how many patients were censored in between; and
the reported total event count closes the system. Guyot et al. (*BMC Med Res Methodol*
2012;12:9) showed that the event and censoring times can be solved back out of those three
ingredients, with an error dominated by the digitisation rather than by the algorithm.

`src/ipd/km_reconstruct.py` takes what a reader measures off the figure:

```python
from ipd import km_reconstruct as km
r = km.reconstruct(points=[(t1, S1), (t2, S2), ...],   # read off the step function
                   risk=[(t, n), ...],                  # the numbers-at-risk row
                   total_events=37)                     # when the paper states it
r['time'], r['status']                                  # one row per patient
```

`data/curves/TEMPLATE_study_group.csv` is the input format and `reconstruct_csv` reads it.

A caveat that decides how much this is worth here. **The primary outcome of the analysis is
binary.** Small-for-size syndrome within 14 to 28 days is never plotted as a survival
curve, so reconstruction can only feed the secondary outcome, mortality or graft loss. Of
the series examined so far, Ogura 2010 is the only one whose survival curves are stratified
by the portocaval gradient (Fig. 6, cut at 9 mmHg), and those panels print no
numbers-at-risk row, which is the ingredient that fixes the cohort size. That figure also
reports a contrast that did not reach significance (88.4% against 74.6% at two years,
p = 0.163) while the split on absolute portal pressure did (p = 0.0003), and the authors
conclude that the gradient does not predict outcome. That has to be met head on rather than
reconstructed around.

### Validation
### Validation

The method is checked against cohorts whose truth is known: simulate, compute the curve,
read 40 points and a six-entry risk table off it, throw the patients away, reconstruct.

| true n | true events | recovered n | recovered events | max error on S(t) | RMSE |
|---|---|---|---|---|---|
| 200 | 101 | 200 | 99 | 0.011 | 0.003 |
| 120 | 96 | 120 | 93 | 0.009 | 0.003 |
| 400 | 137 | 400 | 135 | 0.005 | 0.001 |

Cohort size comes back exactly, because the risk table fixes it. Event counts come back to
within three. The reconstructed curve sits on top of the digitised one to about one
percentage point, which is finer than anyone can read a published figure.

### What is still needed

**No curve has been digitised yet, and the better route turned out not to need one.** Three
of the series publish patient-level haemodynamics directly: Botha 2010 as a table (used in
section 1), Chan 2011 as a scatter of graft portal flow against portal pressure for all 46
recipients (Figs. 1C and 3C), and Ishizaki 2012 as a before-and-after line plot of the
gradient for 25 recipients (Fig. 2). Digitising a scatter of points recovers exactly the
pair of loads the model is about, where digitising a survival curve recovers only the
outcome. Chan is the one worth doing first, because it is the only source that would give
z_R for individual recipients rather than for group means.

---

## What this does not fix

Reconstruction and integration both work on information that exists. Neither produces a
within-group contrast where the paper reports only a mean, and neither substitutes for
patient-level data held by the centres. The honest statement of what is missing is
unchanged: discrimination, individual calibration and clinical utility need real individual
records, and two or three centres sharing them would settle in one analysis what this layer
can only bound.


---

## 3. One series that published its patients

Botha et al. (*Liver Transpl* 2010;16:649-657) Table 5 lists all sixteen recipients of a
left lobe graft with a hemiportocaval shunt: the portocaval gradient with the shunt clamped
and ten minutes after unclamping, the graft-to-recipient weight ratio, encephalopathy,
small-for-size syndrome, follow-up and graft loss. The meta-analysis collapsed that into
one row with a median. `data/patient_level/botha2010.tsv` is the transcription and
`src/ipd/patient_level.py` the analysis.

Eleven of the thirteen statistics the paper reports elsewhere in its own tables and text are
recovered exactly from the transcribed rows. The one substantive disagreement is in the
paper: Table 4 gives the minimum pre-shunt gradient as 12 mmHg, while Table 5 lists a
recipient at 10. It is reported rather than corrected.

The shunt is a deliberate, within-patient manipulation of the resistance load, and it is
the cleanest such manipulation in the whole data set. The gradient fell in every one of the
ten recipients in whom it was measured on both sides, from a median of 17.5 to 5 mmHg,
which is a pressure load falling from 3.5 to 1.0 while the graft and its inflow stay where
they are.

| achieved gradient (mmHg) | 1 | 4 | 4 | 5 | 5 | 5 | 5 | 8 | 10 | **15** |
|---|---|---|---|---|---|---|---|---|---|---|
| z_P | 0.2 | 0.8 | 0.8 | 1.0 | 1.0 | 1.0 | 1.0 | 1.6 | 2.0 | **3.0** |
| SFSS | | | | | | | | | | **yes** |
| graft loss | | | | | | | | | | **yes** |

The series' only case of small-for-size syndrome and its only graft loss are the same
recipient, and that recipient has the highest achieved pressure load of the ten, at
z_P = 3.0 against a median of 1.0. The authors say as much in their discussion: a gradient
that is too high can lead to small-for-size syndrome, and in their series it did, at
15 mmHg.

The statistic has to match the sample. With one event among ten recipients no model can be
fitted, so the question is asked exactly: under independence between outcome and load, the
probability that the single event lands on the highest of ten gradients is 1/10. So p = 0.10
and nothing stronger. What the row is worth is not its p value but its structure: it is a
within-patient contrast, on the model's own variable, in a cohort where the flow load was
held roughly constant by design.


---

## 4. The series that reports its outcomes against its haemodynamics, and finds nothing

Chan et al. (*Liver Transpl* 2011;17:115-121) Table 3 splits all 46 recipients by Clavien
grade and reports graft portal flow and portal pressure for each side of the split. That is
**12 events in 46 recipients**, four times the two in-hospital deaths the main extraction
scores for this series, and it requires no digitising.

| | uncomplicated (n = 34) | Clavien grade 2+ (n = 12) |
|---|---|---|
| graft portal flow per 100 g | 325 mL/min | 311 mL/min |
| portal pressure after implantation | 14 mmHg (10-26) | 14 mmHg (11-18) |
| z_F | 4.01 | 3.84 |
| z_P | 1.2 | 1.2 |
| z_R | 0.30 | 0.31 |

Three things follow, and the first two are uncomfortable.

**The contrast is a null.** The medians are the same to the nearest mmHg, so the two
outcome groups differ in neither load. Formally the series cannot enter a slope model at
all: with no variation in the predictor between the groups there is nothing for a slope to
estimate, which is why it does not appear in the hierarchical fit.

**The extremes of pressure sat among the uncomplicated.** The recipients without a
significant complication spanned 10 to 26 mmHg; those with one spanned 11 to 18. The
highest pressures in the cohort belonged to patients who did well. The authors say the same
thing in their own words: complications were not related to portal flow, flow per 100 g, or
pressure before or after implantation.

**It is also the strongest single case for the resistance reading.** This is the series with
the widest venous outflow in the whole set, a right lobe including the middle hepatic vein
with direct caval anastomosis, and it carries a flow load of 4.0 at a pressure load of 1.2:
a resistance load of 0.30, the lowest of the five paired groups. The authors attribute it to
outflow: in their reading, residual portal hypertension after implantation is a
manifestation of insufficient venous outflow, and a graft that drains freely raises the
threshold at which sinusoidal injury appears.

Both readings have to be carried at once. A series where outflow was maximised shows a very
low resistance load and good outcomes across a wide pressure range; the same series shows no
within-cohort gradient-to-outcome association. The first supports the identity, the second
limits how much the pressure load can be asked to carry on its own.

### What is still worth digitising

Figs. 1C and 3C plot graft portal flow, total and per 100 g, against portal pressure after
implantation for each of the 46 recipients. Those are the only published points that would
give **z_R for individual recipients** rather than for group means, which no other source in
the set can do. The limitation is that the figures carry pressure rather than gradient, and
the individual central venous pressures are not published, so a per-patient z_P needs the
same assumed CVP already flagged in the extraction. Thirty to fifty points read off each
panel with WebPlotDigitizer would be enough.
