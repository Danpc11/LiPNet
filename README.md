# Liver pressure index
### Portal pressure integrates haemodynamic risk after living donor liver transplantation

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)
[![Tests](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/liver_pressure_index/actions/workflows/tests.yml)
[![Calculator](https://img.shields.io/badge/Calculator-live%20on%20GitHub%20Pages-4285F4?logo=googlechrome&logoColor=white)](https://danpc11.github.io/liver_pressure_index/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22861276.svg)](https://doi.org/10.5281/zenodo.22861276)
[![License](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-orange)](LICENSE)

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
Both equal 1 in a healthy donor. In 12 published living donor liver transplantation (LDLT) series, `zP` orders small-for-size syndrome (SFSS) and mortality across groups (13 groups, Spearman ρ = 0.79, p = 0.001) and `zF` does not (11 groups, ρ = 0.35, p = 0.29; ρ = 0.23, p = 0.46 if the two Vasavada 2014 groups defined by a flow cut-off are included). The two coincide only when the graft keeps the donor's outflow resistance; `zP/zF` estimates that ratio. Pooled logistic fit on five SFSS groups (362 recipients, 46 events): logit(SFSS) = −6.40 + 1.88·zP, 5% risk at zP = 1.8 and 10% at 2.2.

<p align="left">
  <img src="assets/graphical_abstract.png" alt="Graphical abstract: a partial liver graft as a fraction of the donor perfusion network; the pressure index zP orders small-for-size risk across published series, the flow index zF does not; the two coincide only at donor outflow resistance" width="410">
</p>

## Calculator

**Live:** https://danpc11.github.io/liver_pressure_index/ (served from this repository by GitHub Pages; the file is `sfss_calculator.html`, it runs entirely in the browser).

Enter PVP, CVP, PVF, graft weight (or GRWR and recipient weight) and the donor flow reference. It returns both indices, the pooled SFSS risk with its bootstrap band, the position relative to the window 1–2, the outflow-resistance ratio `zP/zF`, and what each inflow-modulation manoeuvre achieved in the published series. The normal gradient (5 mmHg) is fixed because the coefficients were fitted on it; a blank CVP is taken as 5 mmHg and flagged, an explicit 0 is a measured zero; graft weight and donor reference must be positive (a blank donor reference falls back to 90 mL/min/100 g and says so); a portal flow of 0 gives z_F = 0 and no `zP/zF` ratio; outside the z_P range of the fitted groups (1.43–2.68) the risk is marked as extrapolation. It is a research tool calibrated on group means from published data, not a validated individual predictor and not a medical device.

## Contents

```
sfss_calculator.html     the calculator (also served at the Pages URL above)
index.html               redirects the Pages root to the calculator
assets/                  graphical abstract and other images for this README
data/series.tsv          22 groups from 12 series, one row per group, source table/page for every value
data/patients.tsv        individual patients where the papers tabulate them (Ou 2010, Yamada 2008, Alim 2016)
data/interventions.tsv   pressure and flow changes produced by each inflow-modulation manoeuvre, with source
data/cohort_template.tsv column layout for a future individual-patient cohort
data/graph*.tsv, certificate_small.tsv   cached outputs of the long model campaigns
src/indices.py           computes zF, zP, Spearman correlations and the pooled logistic fit
src/plots.py             14 plots, each saved on its own as results/<name>.pdf and .png
src/build_app.py         regenerates sfss_calculator.html from the fitted coefficients and the series
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

Plots: `nomogram`, `risk_curve`, `series_pressure`, `series_flow`, `plane`, `resection`, `load_curve`, `interventions`, `network_2d`, `network_3d`, `scaling`, `allometry`, `shear_profile`, `certificate`. Each is a stand-alone figure with its own axes and legend, so they can be combined freely.

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

**Sensitivity** (`results/sensitivity.tsv`, produced by `indices.py`): main analysis, `zP` 13 groups ρ = 0.79 (p = 0.001), `zF` 11 groups ρ = 0.35 (p = 0.29); including the cut-off groups, `zF` 13 groups ρ = 0.23 (p = 0.46); restricted to groups with no imputed or filled value, `zP` 5 groups ρ = 0.70 (p = 0.19) and `zF` 9 groups ρ = 0.46 (p = 0.21), same directions but not significant. The pooled logistic fit cannot be repeated without imputed values: of its five SFSS groups only Yamada 2008 has fully reported pressures (and no events), so **the fitted coefficients rest on groups whose CVP was assumed at 5 mmHg**. Each mmHg of CVP shifts z_P by 0.2. This is the main limitation of the analysis and the reason the calculator is a research tool.

Series: Troisi 2003 (Liver Transpl 9:S36), Troisi 2005 (Am J Transplant 5:1397), Ou 2010 (Transplant Proc 42:876), Vasavada 2014 (Int J Surg 12:177), Alim 2016 (Liver Transpl 22:1643), Chan 2011 (Liver Transpl 17:115), Yagi 2005 (Liver Transpl 11:68), Yagi 2006 (Transplantation 81:373), Wang 2014 (Surg Today 45:979), Osman 2017 (Hepatol Res 47:293), Ogura 2010 (Liver Transpl 16:718), Yamada 2008 (Am J Transplant 8:847).

## Adding a cohort

Fill `data/cohort_template.tsv` (one row per patient). `src/indices.py` exposes `zF()` and `zP()`; a per-patient logistic fit replaces the pooled group fit in `fit_logistic()`, and `python src/build_app.py` then updates the calculator (and, after a push, the live page) with the new coefficients and band.

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
