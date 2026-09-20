# Liver pressure index
### Portal pressure integrates haemodynamic risk after living donor liver transplantation

Two dimensionless indices of haemodynamic load on a partial liver graft, the published series used to test them, the code that reproduces every plot as a separate file, and the calculator.

```
zF = graft PVF per 100 g / donor PVF per 100 g     (donor reference from the same series if measured, otherwise 90 mL/min/100 g)
zP = (PVP − CVP) / 5 mmHg                          (CVP = 5 when not reported; flagged in the data)
```

Both equal 1 in a healthy donor. In 12 published living donor liver transplantation (LDLT) series, `zP` orders SFSS and mortality across groups (13 groups, Spearman ρ = 0.79, p = 0.001) and `zF` does not (ρ = 0.23, p = 0.46). The two coincide only when the graft keeps the donor's outflow resistance; `zP/zF` estimates that ratio.

Repository: https://github.com/Danpc11/liver_pressure_index

## Contents

```
sfss_calculator.html     the calculator (open in a browser; no server, no dependencies)
data/series.tsv          22 groups from 12 series, one row per group, source table/page for every value
data/patients.tsv        individual patients where the papers tabulate them (Ou 2010, Yamada 2008, Alim 2016)
data/cohort_template.tsv column layout for a future individual-patient cohort
data/interventions.tsv   pressure and flow changes produced by each inflow-modulation manoeuvre, with source
data/graph*.tsv, certificate_small.tsv   cached outputs of the long model campaigns
src/indices.py           computes zF, zP, Spearman correlations and the pooled logistic fit
src/plots.py             14 plots, each saved on its own as results/<name>.pdf and .png
src/build_app.py         regenerates sfss_calculator.html from the fitted coefficients and the series
src/app_template.html    the calculator with placeholders
src/model/               the perfusion-network model (2D, 3D, closed-form scaling) and its campaigns
run_all.sh               indices -> plots -> app
```

## Reproduce

```bash
pip install -r requirements.txt
./run_all.sh                       # everything, ~2 min
python src/plots.py risk_curve     # or any single plot
python src/plots.py                # prints the list of plots
```

Plots: `nomogram`, `risk_curve`, `series_pressure`, `series_flow`, `plane`, `resection`, `load_curve`, `interventions`, `network_2d`, `network_3d`, `scaling`, `allometry`, `shear_profile`, `certificate`. Each is a stand-alone figure with its own axes and legend, so they can be combined freely.

## Data

`data/series.tsv` columns: `study`, `group`, `n`, `source` (table or section and page in the original PDF), `GRWR`, `PVF_per_100g`, `donor_PVF_per_100g_ref`, `PVP`, `CVP`, `CVP_measured`, `zF`, `zP`, `outcome`, `outcome_type` (`SFSS_or_dysfunction` or `mortality_or_graft_loss`), `events`, `pct`, `notes`. Every value is a transcription from the publication; nothing is estimated except CVP where flagged.

Series: Troisi 2003 (Liver Transpl 9:S36), Troisi 2005 (Am J Transplant 5:1397), Ou 2010 (Transplant Proc 42:876), Vasavada 2014 (Int J Surg 12:177), Alim 2016 (Liver Transpl 22:1643), Chan 2011 (Liver Transpl 17:115), Yagi 2005 (Liver Transpl 11:68), Yagi 2006 (Transplantation 81:373), Wang 2014 (Surg Today 45:979), Osman 2017 (Hepatol Res 47:293), Ogura 2010 (Liver Transpl 16:718), Yamada 2008 (Am J Transplant 8:847).

## Adding a cohort

Fill `data/cohort_template.tsv` (one row per patient). `src/indices.py` shows the two functions `zF()` and `zP()`; a per-patient logistic fit replaces the pooled group fit in `fit_logistic()`, and `python src/build_app.py` then updates the calculator with the new coefficients and band.

## Calculator

`sfss_calculator.html` takes PVP, CVP, PVF, graft weight (or GRWR and recipient weight) and the donor flow reference, and returns both indices, the pooled SFSS risk with its bootstrap band, the position relative to the window 1–2, the outflow-resistance ratio `zP/zF`, and what each inflow-modulation manoeuvre achieved in the published series. It is a research tool calibrated on group means from published data; it is not a validated individual predictor.

## Model

`src/model/network.py` builds the liver as a graph (hilum → portal tree → hexagonal lobules → hepatic venous tree), solves flows and dissipation for prescribed lobule demand, finds the dissipation-optimal tree at fixed maintenance cost `Σ m_e^b`, and integrates the local shear set-point rule. `network3d.py` does the same on a hemisphere. `scaling.py` gives the closed form `D* ∝ F² C^(−2/b) N^((2−b)/b) L³`, the set-point `τ0 = √(D*/C)` and the allometric closure. `graph_scaling.py`, `network3d.py <R>` and `certify.py` are the campaigns whose outputs are cached in `data/`.

## Funding

DGAPA-PAPIIT IN234029; SECIHTI CBF-2025-G-789.

## Licence

[PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0): free to use, modify and share for research, teaching, personal and non-profit purposes (universities, public research organisations and hospitals included, whatever their funding source); any commercial use requires a separate licence from the authors. Values in `data/` are transcribed from the cited publications and remain the property of their authors. The calculator is a research tool, not a medical device.
