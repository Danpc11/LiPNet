# LiPNet - Liver Perfusion Network model

### A physics-based model of hepatic perfusion and haemodynamic load in partial liver grafts

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/LiPNet/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/LiPNet/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/LiPNet/)
![Version](https://img.shields.io/badge/version-0.3.0-1f6feb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22861276.svg)](https://doi.org/10.5281/zenodo.22861276)

---

## What this is

LiPNet treats the liver as a vascular transport network and asks a single question: what turns portal flow into portal pressure?

The answer is one identity, and it is the result the repository exists to support:

$$\boxed{z_P = z_F \, z_R}\qquad\text{pressure load}=\text{flow load}\times\text{resistance load}$$

$$z_F=\frac{\text{graft PVF per }100\,\text{g}}{\text{donor PVF per }100\,\text{g}},\qquad
z_P=\frac{\text{PVP}-\text{CVP}}{5\ \text{mmHg}},\qquad
z_R=\frac{z_P}{z_F}=\frac{R_{\text{graft}}}{R_{\text{donor}}}$$

All three loads equal 1 in a healthy donor. Flow and pressure therefore carry the same information only while $z_R=1$: whenever the venous outflow is enlarged, $z_R$ falls and a graft can take three or four times the donor's flow without the pressure that would normally come with it. In the five published groups that report both loads, the four with reconstructed or enlarged outflow had $z_R = 0.31,\;0.54,\;0.55,\;0.74$, all below one, with adverse outcomes of 0 to 10%.

The derivation, the four predictions that follow from it and all the numbers are in **[THEORY.md](THEORY.md)**.

<p align="center">
  <img src="assets/graphical_abstract.png"
       alt="Graphical abstract of LiPNet showing the liver as a perfusion network and the relation between flow load, resistance load and pressure load"
       width="520">
</p>

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

$$\mathrm{logit}\,p(z)=\mathrm{logit}\,p_{\mathrm{ref}}+\mu_\beta(z-z_{\mathrm{ref}}).$$

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

## Included clinical series

Thirty-one groups from seventeen published series of adult living-donor liver transplantation. Every value in `data/series.tsv` carries its source table and page, and a provenance label.

| Series | Journal | Link |
|---|---|---|
| Troisi 2003 | Liver Transpl 2003;9(9):S36–41 | [10.1053/jlts.2003.50200](https://doi.org/10.1053/jlts.2003.50200) |
| Troisi 2005 | Am J Transplant 2005;5:1397–1404 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Troisi+hemi-portocaval+shunts+inflow+modulation+small-for-size+2005) |
| Yagi 2005 | Liver Transpl 2005;11(1):68–75 | [10.1002/lt.20317](https://doi.org/10.1002/lt.20317) |
| Yagi 2006 | Transplantation 2006;81(3):373–378 | [10.1097/01.tp.0000198122.15235.a7](https://doi.org/10.1097/01.tp.0000198122.15235.a7) |
| Yamada 2008 | Am J Transplant 2008;8(4):847–853 | [10.1111/j.1600-6143.2007.02144.x](https://doi.org/10.1111/j.1600-6143.2007.02144.x) |
| Botha 2010 | Liver Transpl 2010;16(5):649–657 | [10.1002/lt.22043](https://doi.org/10.1002/lt.22043) |
| Ogura 2010 | Liver Transpl 2010;16(6):718–728 | [10.1002/lt.22059](https://doi.org/10.1002/lt.22059) |
| Ou 2010 | Transplant Proc 2010;42(3):876–878 | [10.1016/j.transproceed.2010.02.064](https://doi.org/10.1016/j.transproceed.2010.02.064) |
| Chan 2011 | Liver Transpl 2011 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Chan+Lo+portal+inflow+pressure+right+liver+living+donor+middle+hepatic+vein) |
| Ishizaki 2012 | Liver Transpl 2012;18(3):305–314 | [10.1002/lt.22440](https://doi.org/10.1002/lt.22440) |
| Vasavada 2014 | Exp Clin Transplant 2014;12(5):437–442 | [journal](https://www.ectrx.org/detail/archive/2014/12/5/0/437/society.php) |
| Wang 2014 | Surg Today 2015;45(8):979–985 | [10.1007/s00595-014-0999-9](https://doi.org/10.1007/s00595-014-0999-9) |
| Alim 2016 | Liver Transpl 2016 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Alim+graft-to-recipient+weight+ratio+MELD+living+donor+liver+transplantation) |
| Uemura 2016 | Surgery 2016;159(6):1623–1630 | [10.1016/j.surg.2016.01.009](https://doi.org/10.1016/j.surg.2016.01.009) |
| Kanetkar 2017 | J Clin Exp Hepatol 2017;7(3):235–246 | [10.1016/j.jceh.2017.01.114](https://doi.org/10.1016/j.jceh.2017.01.114) |
| Osman 2017 | Hepatol Res 2017;47(4):293–302 | [10.1111/hepr.12727](https://doi.org/10.1111/hepr.12727) |
| Yao 2018 | Transplantation 2018;102(4):623–631 | [10.1097/TP.0000000000002047](https://doi.org/10.1097/TP.0000000000002047) |

Three of these links are PubMed searches rather than DOIs, because the bibliographic record has not been re-verified; the extraction itself is sourced in `data/series.tsv`.


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

The seventeen series and their links are listed in [Included clinical series](#included-clinical-series) above.

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

## Scope

LiPNet is a mechanistic framework, not an individual risk score. It gives a common normalised language for portal hyperperfusion, venous outflow reconstruction, small-for-size physiology, portal hypertension, hepatic resection and vascular obstruction, and it says which part of that language transfers between centres: the change in risk with a change in load does, the baseline risk does not.

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
