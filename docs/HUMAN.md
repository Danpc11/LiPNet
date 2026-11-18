# The network model in a human liver: units, growth and regrowth

`src/model/liver.py` works in the units of the optimisation. `src/model/human.py` puts the same optimum in
millimetres, mL/min, mmHg and Pa, sets it against what has been measured in people, and makes the lobule layer
(`src/fem/`) and the network describe one liver. Every number below is produced by

```bash
python src/model/human.py     # -> results/model_human.json
python tests/test_human.py
```

## 1. Calibration

| Quantity | Value | Source |
|---|---|---|
| liver | 1500 g adult, 1900 g for the cast, density 1.05 g/mL | the cast is one liver of about 1.9 kg [1] |
| portal flow | 90 mL/min per 100 g | donor reference of the manuscript; total hepatic flow is 100 to 130 mL/min per 100 g, one fifth to one third of it arterial [5] |
| arterial share | 0.25 | as in `liver.py` [5] |
| portocaval gradient | 5 mmHg | reference of the manuscript; about 5 mmHg [5] |
| viscosity | 3.5 mPa s | as in `src/fem/params.py` [2] |
| lobules | $3^{13}=1.6\times10^6$, 14 generations, branching 3, $d=3$ | volume over lobule gives $1.1\times10^6$; the hepatic veins of the cast extrapolate to 1.1 to 4.0 million terminal vessels [1] |
| branch length | 0.44 of the side of the territory supplied | portal generations 7 to 11 of the cast: 0.35 to 0.50 [1] |
| lobular drop | 1.0 mmHg at donor flow, imposed | section 4 |

The branch length is the one geometric number that had to change. The canopy made every branch as long as the side
of its territory, which gives a portal trunk of 122 mm and lengths two to three times those of the cast. With the
factor 0.44 the trunk of an adult liver is 50 mm long; the pooled length of the adult portal vein is 48.6 mm (35.6
to 61.6) [8]. The factor multiplies every length, so the shares of the pressure drop of `liver.py` are unchanged.

## 2. The exponent of the portal vein

Three readings of the same cast, none of them fitted to transplant data:

| Reading | $b$ |
|---|---|
| radius law of the optimum on the measured radii, lengths and counts of generations 7 to 11 | 0.89 (0.83 to 0.94 leaving one generation out) |
| absolute radii, trunk to generation 11, lobular drop imposed | 0.76 (rms error below 10% for 0.69 to 0.86) |
| absolute radii, lobule left to the optimum | 0.74 (0.67 to 0.82) |

The second reading uses no free scale: the liver volume, 90 mL/min per 100 g, 5 mmHg and 3.5 mPa s fix every radius.

| vessels in the generation | 1 | 216 | 432 | 1339 | 3715 | 10541 | 23976 |
|---|---|---|---|---|---|---|---|
| cast radius (mm) | 7.34 | 1.08 | 0.718 | 0.482 | 0.309 | 0.231 | 0.175 |
| model / cast, $b=0.5$ | 1.93 | 1.27 | 1.42 | 1.30 | 1.30 | 1.11 | 1.02 |
| model / cast, $b=0.75$ | 1.04 | 0.92 | 1.06 | 1.03 | 1.09 | 0.98 | 0.95 |
| model / cast, $b=1$ | 0.70 | 0.79 | 0.94 | 0.96 | 1.07 | 1.01 | 1.01 |

Neither end of $[0.5,1]$ reproduces the cast (rms error 0.33 and 0.17). With branches as long as their territory
no value of $b$ does. `B_HUMAN = 0.75` is used below; the range $[0.5,1]$ of the repository is not narrowed. The
trunk of the cast is 14.7 mm wide; ultrasound of the living portal vein reads 10.7 mm and puts it at $b=0.93$
(section 6).

## 3. The hepatic veins of the cast

Poiseuille drop of each measured generation at donor flow, mmHg:

| | generations | drop per generation | ratio per generation | wall shear (Pa) |
|---|---|---|---|---|
| portal vein | 7 to 11 | 0.14, 0.14, 0.18, 0.17, 0.13 | 0.98 | 0.79 to 0.99 (0.32 in the trunk) |
| hepatic vein | 6 to 8 | 0.42, 0.78, 1.35 | 1.79 | 2.4 to 5.4 |

The whole portal tree of the cast (17 to 18 generations) holds 1.9 to 2.0 mmHg, 0.9 to 1.0 mmHg if its radii were
20% larger in life. Three generations of the hepatic veins hold 2.6 mmHg, and their trend gives 56 to 99 mmHg for
the 13 to 14 generations of the tree, 27 to 48 mmHg with radii 20% larger. A liver with a gradient of 5 mmHg cannot
have these calibres at these counts. The exponent they give (0.46 here, 0.45 from the ratios quoted in `liver.py`)
is reported and no result uses it. The hepatic artery gives $b$ between 1.7 and 2.2, outside the cost law; the
artery enters the model only as a share of the flow.

## 4. One lobule for both layers

Until this release the two layers disagreed about the lobule. The lobule layer imposed 4 mmHg between the portal
triad and the central vein, which at the hepatic flow of the network is 4.9 mmHg, the whole gradient. The network
left its sinusoidal edge to the optimum, which at the size of a human liver gives 0.2 to 0.6 mmHg (4 to 12% for $b$
between 0.93 and 0.6; 21% on the explicit graph of 151 lobules).

Both now hold `LOBULE_DROP_MMHG = 1.0`, defined once in `src/fem/params.py` and imposed on the sinusoidal edge by
`human.liver`. The series of an adult liver is then 2.1 mmHg along the portal tree, 1.0 across the lobule and 1.9
along the hepatic veins:

| | portal tree | lobule | hepatic veins |
|---|---|---|---|
| model, lobular drop imposed, $b=0.75$ | 0.43 | 0.20 | 0.37 |
| model, lobule at the optimum, $b=0.75$ | 0.49 | 0.08 | 0.43 |
| Poiseuille on the human cast [1] | 0.38 to 0.40 | | |
| rat, micropuncture [6] | 0.47 | 0.44 | 0.09 |
| rat and dog, micropuncture [7] | | | 0.44 and 0.31 |
| lobule model of Debbaut et al. [3] | | 0.95 mmHg imposed | |

The permeability of the lobule layer is a calibration: $k_r=1.6\times10^{-13}$ m² carries 1 mL/min per g at 1 mmHg.
The permeability measured on the human cast, $1.56\times10^{-14}$ m² radially [2], carries a tenth of the
physiological flow at that drop in the lobule model of the same authors (157 of about 1500 mL/min) [3]. Darcy flow
depends on $k\,\Delta P/\mu$ only, so every flow, velocity, shear, transit time and concentration of the lobule
layer is identical to the previous release; its pressures are not.

The choice of 1 mmHg moves little. Over 0.35 mmHg (the optimum) to 2.2 mmHg (the rat by micropuncture) the portal
radii stay within 20% of the cast and the resistance load of a lobar graft with the arterial inflow stays
between 0.82 and 0.90.

## 5. The adult liver

Liver of 1500 g, lobular drop imposed:

| $b$ | portal trunk radius (mm) | mean velocity in the trunk (cm/s) | wall shear, trunk to terminal branch (Pa) | volume of the two trees (mL) | $z_R$ of a lobar graft | with the arterial inflow, $\varphi\le1$ |
|---|---|---|---|---|---|---|
| 0.60 | 9.9 | 7 | 0.10 to 2.8 | 201 | 0.997 | 0.86 to 0.90 |
| 0.67 | 8.4 | 10 | 0.17 to 2.4 | 165 | 0.994 | 0.85 to 0.90 |
| 0.75 | 7.0 | 15 | 0.29 to 2.0 | 136 | 0.989 | 0.85 to 0.89 |
| 0.85 | 5.8 | 21 | 0.50 to 1.6 | 117 | 0.976 | 0.84 to 0.88 |
| 0.93 | 5.2 | 27 | 0.73 to 1.2 | 109 | 0.960 | 0.82 to 0.86 |

Measured in adults: portal diameter 10.7 mm by ultrasound and 13.3 mm by CT (pooled) [8], 14.7 mm in the cast [1];
mean velocity 15 ± 4 and 16.5 ± 4.9 cm/s; flow 648 ± 186 and 899 ± 284 mL/min [10], below the 1350 mL/min that 90
mL/min per 100 g gives a liver of 1500 g. Liver volumes by body weight in section 6 follow the standard liver volume of Urata [9]. The calibration does not move the result of `docs/LIVER.md`, section 6:
the five observed resistance loads (0.31 to 0.74) lie below every value in the table.

## 6. Growth with the body

`growth(b, s)` scales the calibrated liver from 150 g to 1500 g at invariant gradient and invariant flow per gram.
The liver can grow by enlarging its lobules or by adding them; $s$ is the share of the growth in volume that is
enlargement (lobule volume $\propto V^{s}$, number of lobules $\propto V^{1-s}$). The sinusoidal bed keeps its
permeability, so the conductance of a lobule grows with its length.

Exponents with liver volume, for $b$ from 0.6 to 0.93:

| | $s=0$ (lobules added) | $s=0.5$ | $s=1$ (lobules enlarged) | measured |
|---|---|---|---|---|
| portal diameter | 0.34 to 0.40 | 0.35 to 0.38 | 0.35 | 0.47 (0.25 to 0.68) [10] |
| mean velocity in the portal trunk | 0.21 to 0.32 | 0.24 to 0.30 | 0.29 | 0.37 (−0.11 to 0.96) [10]; no difference between age groups [11] |
| lobular drop | 0 | 0.33 | 0.67 | not measured |
| velocity in the sinusoid | 0 | 0.17 | 0.33 | not measured |

The portal trunk does not tell how the liver grows: its diameter rises as $V^{0.34}$ to $V^{0.40}$ and its velocity
as $V^{0.21}$ to $V^{0.32}$ whatever $s$. At $b=0.75$ the diameter goes from 6.1 mm to 14.0 mm and the velocity
from 7 or 8 to 15 cm/s. The measurements available as ranges by body weight [10] contain these exponents; the pooled
portal diameter is 6.6 mm in children and 10.7 mm in adults [8]. In absolute value the ultrasound diameters sit on
the curve of $b=0.93$ and below that of $b=0.75$, which reproduces the cast:

| liver volume (mL) [9] | 155 | 334 | 551 | 784 | 1269 |
|---|---|---|---|---|---|
| portal diameter by ultrasound (mm) [8,10] | 3 to 4 | 4 to 7 | 4 to 8 | 6 to 9 | 10.7 |
| model, $b=0.75$ | 6.7 | 8.7 | 10.2 | 11.5 | 13.5 |
| model, $b=0.93$ | 4.9 | 6.4 | 7.5 | 8.5 | 9.9 |

The cast and the living vein give the two ends of the interval for the portal tree, 0.75 to 0.93. What $s$ changes is the lobule: if a newborn had
the lobules of an adult in number and they only enlarged ($s=1$), its lobular drop would be 0.22 mmHg and its
sinusoidal velocity 0.46 of the adult's; if lobules were added ($s=0$) both would be the adult's. In the mouse the
number of lobules doubles after birth and stops rising in the second week, and lobule size is final at four weeks
[12]; in the rat the liver lobe grows 12-fold while lobule area grows more than 5-fold and the number of lobules at
the surface by about 30% [13], which is $s$ close to one. No measurement of lobule size or number against age in
people was found.

## 7. Regrowth of a graft

After partial hepatectomy the rat liver regrows by enlarging its lobules; new lobules do not form [13,14]. In
people the implanted right lobe reaches 1.58 times its volume at two weeks (1.7 ± 0.3 in a second series), the
standard volume of the recipient, and then settles to 1.50, 1.44 and 1.31 times at one, three and six months [15,16].

`regeneration(rho, ...)` gives the loads of a graft whose volume has grown $\rho$-fold by lobule enlargement at
constant inflow. Lengths grow as $\rho^{1/3}$. The lobule keeps its permeability, so its resistance falls as
$\rho^{-1/3}$. The resistance of the trees grows as $\rho^{1/3}$ if their calibres do not change and falls as
$1/\rho$ if the calibres grow with the lengths. With the lobule holding one fifth of the drop,

$$\frac{z_F}{z_{F,0}}=\frac1\rho,\qquad \frac{z_P}{z_{P,0}}=0.8\,\rho^{\pm}+0.2\,\rho^{-1/3},\qquad \rho^{\pm}=\rho^{1/3}\ \text{or}\ \rho^{-1}.$$

| months | volume | flow load | pressure load, calibres unchanged | pressure load, calibres in proportion | resistance load, unchanged | in proportion |
|---|---|---|---|---|---|---|
| 0.5 | ×1.58 | ×0.63 | ×1.10 | ×0.68 | ×1.74 | ×1.07 |
| 1 | ×1.50 | ×0.67 | ×1.09 | ×0.71 | ×1.64 | ×1.06 |
| 3 | ×1.44 | ×0.69 | ×1.08 | ×0.73 | ×1.56 | ×1.06 |
| 6 | ×1.31 | ×0.76 | ×1.06 | ×0.79 | ×1.39 | ×1.04 |

Regrowth lowers the flow load by itself. It lowers the pressure load only if the vessels widen: with the calibres of
reperfusion the pressure load rises by 6 to 10%. A graft of one third that has regrown threefold with calibres in
proportion returns to the donor's flow per gram with a pressure load of 1.22 and a sinusoidal velocity 1.44 times
the donor's, because it carries the flow of a whole liver through a third of the lobules.

## 8. Limits

The anatomy is one cast of one liver, made after machine perfusion, with estimated vessel counts. Flows are the
references of the manuscript, viscosity is constant and vessels are rigid tubes. The share of the pressure drop by
compartment has been measured in the rat and the dog, not in people. The paediatric portal measurements are ranges
by body weight. The series of portal flow and pressure during the regrowth of human grafts [17] could not be read
and are not used. The growth of human lobules with age is not known to this note.

## Sources

1. Debbaut C, et al. J Anat 2014;224:509-517 (human liver cast, Table 1: radius, length and estimated number of vessels per generation). Macroscopic radii of generations 1 to 5 from [4].
2. Debbaut C, et al. Perfusion characteristics of the human hepatic microcirculation based on three-dimensional reconstructions and computational fluid dynamic analysis. J Biomech Eng 2012;134:011003.
3. Debbaut C, et al. A 3D porous media liver lobule model: the importance of vascular septa and anisotropic permeability for homogeneous perfusion. Comput Methods Biomech Biomed Engin 2014;17:1295-1310.
4. Debbaut C, et al. Constructing a model of the blood flow in the human liver using contrast-loaded corrosion casting and micro-CT imaging (conference abstract, Ghent University).
5. Lautt WW. Hepatic Circulation: Physiology and Pathophysiology, chapter 2 (NCBI Bookshelf NBK53069).
6. Shibayama Y. Kanzo 1975;16:145-152 (normal rat: portal vein 110, terminal portal venule 68, central vein 28, vena cava 20 mmH2O).
7. Bohlen HG, Maass-Moreno R, Rothe CF. Hepatic venular pressures of rats, dogs, and rabbits. Am J Physiol 1991;261:G539-G547.
8. Chau P, Yoon JS, Moses D, Pather N. A systematic review and meta-analysis of portal vein morphometry in pediatric and adult populations. Eur J Radiol 2023;168:111016.
9. Urata K, et al. Calculation of child and adult standard liver volume for liver transplantation. Hepatology 1995 (doi 10.1002/hep.1840210515); used in `urata_volume`.
10. López Barrio AM, et al. An Esp Pediatr 1996;44:45-49 (82 healthy children; adult values of Moriyasu and Ohnishi as quoted there).
11. Ohtsuka Y. J Jpn Soc Pediatr Surg 1994;30:1267-1277 (106 normal children, mean portal velocity 15.4 ± 3.5 cm/s).
12. Azizoglu DB, et al. Vascularization of neonatal liver lobules presages adult liver size. Nat Commun 2025;16:9989.
13. Papp V. Architectural changes in the rat liver during ontogenesis, regeneration and oncogenesis. Doctoral thesis, Semmelweis University, 2013.
14. Tsomaia K, et al. Liver structural transformation after partial hepatectomy and repeated partial hepatectomy in rats. World J Gastroenterol 2020;26:3899-3916.
15. Zhang Y, et al. Comparison of liver regeneration between donors and recipients after adult right lobe living-donor liver transplantation. Quant Imaging Med Surg 2022;12:3184-3192.
16. Byun SH, Yang HS, Kim JH. Liver graft hyperperfusion in the early postoperative period promotes hepatic regeneration 2 weeks after living donor liver transplantation. Medicine 2016;95:e5404.
17. Eguchi S, et al. Liver Transpl 2003;9:547-551. Jiang SM, et al. Liver Transpl 2009;15:1043-1049. Garcia-Valdecasas JC, et al. Liver Transpl 2003 (PMID 12783396). Yagi S, et al. Liver Transpl 2005 (doi 10.1002/lt.20317). Wu TJ, et al. Liver Transpl 2011 (doi 10.1002/lt.22326).
