# The figures of the article

Every figure of the article is drawn by one script in one style:

```bash
python src/paper_figures.py           # the ten figures -> results/Figure_1 ... Figure_4, Figure_S1 ... Figure_S6 (.pdf and .png)
python src/paper_figures.py 2 S3      # only those
```

`src/figstyle.py` holds the style; `src/paper_figures.py` holds the figures. No number is typed into a figure: each
comes from `data/`, from `results/*.json` or from the model.

## Style

| | Rule |
|---|---|
| Size | 180 mm wide (two columns) or 88 mm (one column). PDF in vector with embedded TrueType fonts, PNG at 300 dpi. Axes are placed in millimetres, so the type has the same size in every figure. |
| Type | One sans family (Liberation Sans, Arial or Helvetica, whichever the machine has). 7 pt for axis labels and ticks, 6.5 pt inside a panel, 6 pt for second lines, 8 pt bold capitals for panel letters. |
| Text | In English, sentence case, units in parentheses. Ink or muted ink only: text never takes the colour of a series. Variables in italics, exponents as superscripts. |
| Frame | No title, no top or right spine, no grid, no legend frame. A legend never covers data; where a label can sit beside its mark it replaces the legend. |
| Lines | Axes 0.5 pt, data 1.0 to 1.2 pt, references 0.5 to 0.7 pt dashed. Markers that can overlap carry a white ring, so that the points stay countable. |
| Statistics | Each group or recipient is drawn, not only a mean. A test or a posterior probability that the legend quotes is also written in the panel. A scatter carries its trend: the fitted curve, the power law or the median. |

## Colour

One meaning per colour in the whole article:

| Colour | Meaning |
|---|---|
| purple `#b33fb0` | portal tree, inlet |
| orange `#eb6834` | sinusoids, lobules, total hepatic flow |
| blue `#2a78d6` | hepatic veins, outlet; also the single series of a panel that has no vessels |
| red `#e34948` | the one recipient singled out for an adverse outcome |
| blue ramp, light to dark | an ordered parameter: $b$, $h$, baseline gradient |
| blue tint | the range the network model reaches |
| warm grey | a reference window or an observed range |

The vessels follow the anatomical atlases: portal vein purple, hepatic veins blue. Green is not used for a vessel
because it is the colour of bile, and red is not because it is the colour of the hepatic artery. The three vascular
colours stay apart for every pair under the three common colour-vision deficiencies (checked with the validator of the
palette), and no panel uses more than three unordered colours. Figure 3 has five strata: the colour is the centre
(blue, orange, aqua) and the symbol the publication.

## What each figure reads

| Figure | Function | Reads |
|---|---|---|
| 1 | `figure_1` | the model itself: `model/network.py` and `model/network3d.py` recompute the optimum at $b=0.75$ for panels B and C |
| 2 | `figure_2` | `data/series.tsv`, `results/predictions.json` (the two intervals of panel B are `P2_discordance.uncertainty`) |
| 3 | `figure_3` | `results/bayes_main.json`, `results/spec_curve.json`, `data/series.tsv` |
| 4 | `figure_4` | `results/predictions.json` |
| S1 | `figure_s1` | `data/graph_scaling_all.tsv`, `data/graph3d_all.tsv`, `model/scaling.py` |
| S2 | `figure_s2` | `model/whole_graft.py` |
| S3 | `figure_s3` | `results/ipd_ecological.json` |
| S4 | `figure_s4` | `ipd/patient_level.py` |
| S5, S6 | `figure_s5`, `figure_s6` | `model/human.py` |

Figures 2, 3, 4 and S3 need the results of `run_all.sh`; the others run on a clean clone.

## Provenance

`READS` in `src/paper_figures.py` is the table above in code. When a figure is drawn, the script writes to
`results/figures_provenance.json` the SHA-256 of every file the figure reads, of `paper_figures.py` and `figstyle.py`,
and of the PDF and the PNG. `python src/manifest.py check` compares that record with the files on disk and names any
figure whose input, script or image is no longer the one recorded. The images carry no date, so the same inputs on the
same environment give the same bytes.

In Figure 2B each group has two intervals: the thick bar is the sampling error of the two group means, the thin bar
adds the uncertainty of the references (donor flow, normal gradient, assumed central venous pressure).
 The stand-alone plots of
`src/plots.py`, `src/fem/` and `src/ipd/` are working plots of the repository and are not figures of the article;
they share its three vascular colours.
