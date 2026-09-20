# liver-pressure-index

Code, data and manuscript for **"Portal pressure gradient, not portal flow, sets the risk of small-for-size syndrome: a perfusion-network model of the liver graft"** (submitted to *Journal of Hepatology*).

The repository has two layers. The clinical layer (Stage 0) reproduces every number and figure in the manuscript from the extraction table of the twelve published series in a few minutes. The physics layer (Stages 1–3) reproduces the network model that motivates the two indices; it takes one to two hours and its outputs are cached in `data/` so that Stage 0 runs without it.

## Repository layout

```
src/                  Python: model, campaigns, figures
manuscript/           make_manuscript.js (docx-js) and the generated Word file
app/                  sfss_calculator.html, self-contained risk calculator (research tool)
data/                 extraction table of the published series and cached model outputs
figures/              Figs. 1–4 and S1–S3 as vector PDF and 300-dpi PNG
assets/               liver pictograms used in Fig. 1C (AI-generated, disclosed in the legend)
notes/                working notes (Spanish) documenting each analysis step
run_all.sh            the whole pipeline, stage by stage
```

## Format

Both Word files follow the Journal of Hepatology guide: Times New Roman 12 pt, double spacing, numbered pages, section headings bold and subsections italic; figures use Liberation Sans (metrically identical to Arial), panel letters bold 14 pt, no text below 5 pt, vector PDF at print size.

## Quick start

```bash
pip install -r requirements.txt
npm install docx            # only for the Word manuscript
./run_all.sh                # Stage 0: series table, logistic fit, Figs. 1–4, manuscript
```

## The two indices

| index | definition | equals 1 in | data needed |
|---|---|---|---|
| `z_F` | graft PVF per 100 g ÷ donor PVF per 100 g | healthy donor | PVF, graft weight, donor PVF/100 g (90 if not measured) |
| `z_P` | (PVP − CVP) / 5 mmHg | healthy donor | PVP, CVP (5 if not measured, flagged) |

The two coincide only when the graft keeps the donor's outflow resistance (`z_P/z_F` estimates that resistance ratio). In the published series `z_P` orders outcomes (13 groups, Spearman 0.79, p = 0.001) and `z_F` does not (0.23, p = 0.46).

## Scripts, in the order they are run

### Stage 0 — clinical analysis (manuscript)

| script | input | output | what it does |
|---|---|---|---|
| `src/series_from_pdfs.py` | the twelve PDFs (values typed in, with page/table of origin in the `fuente` column) | `results/series_pdfs_z.tsv` | Builds the 22-group table with `z_F`, `z_P`, outcome counts. Every value is a transcription; nothing is estimated except CVP = 5 where flagged. |
| pooled logistic (inline in `run_all.sh`) | `series_pdfs_z.tsv` | `results/logit_zP.json` | Binomial logistic fit of SFSS on `z_P` for the five groups with event counts; 2000 bootstrap resamples for the CIs and the risk band. |
| `src/make_figs_jhep.py` | `series_pdfs_z.tsv`, `logit_zP.json`, `partB_states.pkl`, `liver_*.png` icons | `results/JHEP_Fig1-4.{pdf,png}` | Fig. 1 model (2D and 3D networks, four liver conditions, load curve); Fig. 2 nomogram, risk curve, effect of modulation; Fig. 3 outcome vs `z_P` and `z_F`; Fig. 4 flow–pressure plane with cirrhosis and resection. |
| `src/make_supp_figs.py` | cached scaling, allometry and certificate tables | `results/JHEP_FigS1-S3.{pdf,png}` | Supplementary figures. |
| `manuscript/make_manuscript.js`, `manuscript/make_supplement.js` | figures, `series_pdfs_z.tsv`, `closure_kruepunga.tsv` | `JHEP_manuscript.docx`, `JHEP_supplement.docx` | Word manuscript in the journal's structure (title page, structured abstract, Impact and implications, IMRaD, table, legends, figures). |
| `app/sfss_calculator.html` | – | – | Stand-alone calculator: enter PVP, CVP, PVF, graft weight or GRWR; returns both indices, the pooled risk with its band, the outflow-resistance ratio and what each manoeuvre achieved in the series. Research tool, not a validated score. |

### Stage 1 — network model of the liver (2D)

| script | what it does |
|---|---|
| `src/liver_canopy.py` | The model. `build_liver_graph` builds hilum → portal tree → hexagonal lobules → hepatic venous tree; `Loaded` solves flows and dissipation for prescribed lobule demand (with optional fluctuating demand `Q = μμᵀ + Cov`); `best_sinusoid_tree` finds the dissipation-optimal tree by coordinate descent with the exact closed-form allocation; `adapt` integrates the local shear set-point rule with pruning; `optimize` is the two-stage softmax/L-BFGS-B allocation used as a control. |
| `src/run_liver.py` | Campaign over maintenance exponent b ∈ {1, 3/4, 2/3, 1/2} and demand fluctuation σ ∈ {0, 1, 2, 4}: best tree, Lorente-type baseline, adapted rest point, sinusoids per lobule, cycle rank, efficiencies. Writes `partB_campaign.tsv`, `partB_states.pkl` (used by Fig. 1A). Also Part A: closed-form efficiency of generation geometries and inference of b from diameter/length ratios. |
| `src/make_figs.py` | Figures for the physics companion paper (network states, loops vs σ, efficiency vs b). |

### Stage 2 — scaling with size

| script | what it does |
|---|---|
| `src/scaling.py` | Closed form for the symmetric space-filling canopy: `D* ∝ F² C^(−2/b) N^((2−b)/b) L³`, set-point `τ0 = √(D*/C)`, vessel mass from `c_e = m_e^b`; exponent fits over 10³–10⁸ lobules and the allometric closures (portal pressure invariance, flow exponents, lobule size). |
| `src/graph_scaling.py R` | Exact tree optimum on the 2D graph at disk radius R; exponent of D*/F² with N. |
| `src/liver3d.py R` | Hemispheric 3D liver (hexagonal prisms, hilum at the flat face), fast incremental tree optimum; exponents in 3D. |
| `src/scaling_figs.py` | Compares graph exponents with the closed form. |

### Stage 3 — controls

| script | what it does |
|---|---|
| `src/certify.py` | Enumerates all 6⁷ = 279 936 spanning trees on a 7-lobule domain and shows the heuristic reaches the exact optimum at every b. |

## Data

`data/series_pdfs_z.tsv` columns: `estudio`, `grupo`, `n`, `fuente` (table and page in the original PDF), `GRWR`, `PVF_100g`, `ref_donante` (donor PVF per 100 g used for `z_F`), `PVP`, `CVP`, `CVP_medida` (True if reported), `desenlace`, `eventos`, `pct`, `z_F`, `z_P`, `notas`. Studies: Troisi 2003, Troisi 2005, Ou 2010, Vasavada 2014, Alim 2016, Chan 2011, Yagi 2005, Yagi 2006, Wang 2014, Osman 2017, Ogura 2010, Yamada 2008.

Other files in `data/` are cached outputs of Stages 1–3 (see the script table).

## Citation

Vázquez-Victorio G, Pérez-Calixto M, Escutia-Guadarrama L, Tovar H, Pérez-Calixto D. Portal pressure gradient, not portal flow, sets the risk of small-for-size syndrome: a perfusion-network model of the liver graft. 2026 (submitted).

Companion physics papers: Vázquez-Victorio et al., *Metabolic cost sets the shear stress profile of optimal vascular trees* (2026); Tovar et al., *Cost convexity controls the architecture of resistance-optimal networks* (2026).

## Funding

DGAPA-PAPIIT IN234029; SECIHTI CBF-2025-G-789.

## Note on AI assistance

Claude (Anthropic) was used for code-script assistance and language editing. All values were verified by the authors against the source publications and the source code.
