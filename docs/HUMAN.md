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
| liver | 1500 g adult, 1900 g for the cast, density 1.05 g/mL | the cast is one liver of about 1.9 kg [1,2] |
| portal flow | 90 mL/min per 100 g | donor reference of the manuscript; five series measure it in the donor: 70, 81, 85, 91 and 118 (section 6) |
| arterial share | 0.25 | as in `liver.py`; 0.22 and 0.31 measured in donors (section 6) |
| portocaval gradient | 5 mmHg | reference of the manuscript; about 5 mmHg [6] |
| viscosity | 3.5 mPa s | as in `src/fem/params.py` [3] |
| lobules | $3^{13}=1.6\times10^6$, 14 generations, branching 3, $d=3$ | volume over lobule gives $1.1\times10^6$; a third as many moves $b$ by 0.04 (section 2) |
| branch length | 0.44 of the side of the territory supplied | portal generations 7 to 11 of the cast: 0.35 to 0.50 [1] |
| lobular drop | 1.0 mmHg at donor flow, imposed | section 4 |

The branch length is the one geometric number that had to change. The canopy made every branch as long as the side
of its territory, which gives a portal trunk of 122 mm. With the factor 0.44 the trunk of an adult liver is 50 mm
long; the pooled length of the adult portal vein is 48.6 mm (35.6 to 61.6) [9]. The factor multiplies every length,
so the shares of the pressure drop of `liver.py` are unchanged. The canopy shortens its branches by the same ratio
at every generation; in the cast the first four generations are all 32 to 35 mm long, so between generations 3 and 6
the model is 20 to 55% short (section 2).

## 2. The exponent of the portal vein

The cast gives radius, length and number of vessels for eleven portal generations: six measured on the whole cast
[2] and five on dissected samples [1]. Readings of $b$, none fitted to transplant data:

| Reading | $b$ |
|---|---|
| radius law of the optimum on the eleven measured generations | 0.73 (0.71 to 0.74 leaving one generation out) |
| the same on generations 1 to 6 alone | 0.77 (0.67 to 0.82) |
| the same on generations 7 to 11 alone | 0.89 (0.83 to 0.94) |
| absolute radii of the eleven generations, lobular drop imposed | 0.74 (rms error below 10% for 0.68 to 0.80) |
| absolute radii, lobule left to the optimum | 0.71 (0.67 to 0.77) |
| absolute radii, $5.3\times10^5$ lobules | 0.70 (0.64 to 0.76) |
| living portal trunk by ultrasound (section 7) | 0.93 |

The radius law uses the measured lengths and numbers and no scale. The absolute radii use the canopy and no free
scale: the liver volume, 90 mL/min per 100 g, 5 mmHg and 3.5 mPa s fix every radius.

| generation | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| vessels | 1 | 2 | 7 | 27 | 80 | 220 | 432 | 1339 | 3715 | 10541 | 23976 |
| cast radius (mm) | 7.3 | 6.3 | 4.0 | 2.1 | 1.5 | 1.1 | 0.718 | 0.482 | 0.309 | 0.231 | 0.175 |
| cast length (mm) | 35 | 35 | 33 | 32 | 27 | 16 | 8.1 | 5.2 | 3.2 | 2.7 | 1.5 |
| model length (mm) | 54 | 37 | 26 | 18 | 12 | 8.6 | 6.0 | 4.1 | 2.9 | 2.0 | 1.4 |
| model / cast radius, $b=0.5$ | 1.94 | 1.66 | 1.52 | 1.61 | 1.41 | 1.24 | 1.42 | 1.30 | 1.30 | 1.11 | 1.02 |
| model / cast radius, $b=0.75$ | 1.04 | 0.93 | 0.91 | 1.04 | 0.96 | 0.90 | 1.06 | 1.03 | 1.09 | 0.98 | 0.95 |
| model / cast radius, $b=1$ | 0.70 | 0.64 | 0.67 | 0.81 | 0.79 | 0.77 | 0.94 | 0.96 | 1.07 | 1.01 | 1.01 |

Neither end of $[0.5,1]$ reproduces the cast (rms error 0.37 and 0.25). With branches as long as their territory
no value of $b$ does. `B_HUMAN = 0.75` is used below; the range $[0.5,1]$ of the repository is not narrowed.

## 3. The hepatic veins and the hepatic artery of the cast

Poiseuille drop of each measured generation at donor flow, mmHg:

| | generations | drop per generation | sum | wall shear (Pa) |
|---|---|---|---|---|
| portal vein | 1 to 11 | 0.02, 0.02, 0.04, 0.12, 0.13, 0.09, 0.14, 0.14, 0.18, 0.17, 0.13 | 1.2 | 0.33 in the trunk to 0.99 |
| hepatic vein | 2 to 8 | 0.10, 0.14, 0.26, 0.43, 0.42, 0.78, 1.35 | 3.5 | 0.45 to 5.4 |

The whole portal tree of the cast (17 to 18 generations) holds 1.9 to 2.0 mmHg, 0.9 to 1.0 mmHg if its radii were
20% larger in life. Seven generations of the hepatic veins hold 3.5 mmHg, and their trend gives 56 to 100 mmHg for
the 13 to 14 generations of the tree, 27 to 48 mmHg with radii 20% larger. A liver with a gradient of 5 mmHg cannot
have these calibres at these counts. The authors of the cast met the same problem: with the trends extrapolated as
measured, their model put the sinusoidal pressure above the portal pressure, and they adjusted the lengths and
numbers of the unmeasured generations to reach 10 mmHg in the portal vein, 5 mmHg in the sinusoids and 3 mmHg at
the outlet [2]. The exponent of the hepatic veins (0.42 on the seven generations, 0.39 to 0.44 leaving one out) is
reported and no result uses it.

The hepatic artery follows Murray's law: $b=1.04$ (0.99 to 1.06) on its eleven measured generations. The artery
enters the model only as a share of the flow.

## 4. One lobule for both layers

Until release 0.4.4 the two layers disagreed about the lobule. The lobule layer imposed 4 mmHg between the portal
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
| Poiseuille on the human cast [1,2] | 0.37 to 0.40 | | |
| pressures taken from the literature by the authors of the cast [2] | 0.71 | | 0.29 |
| rat, micropuncture [7] | 0.47 | 0.44 | 0.09 |
| rat and dog, micropuncture [8] | | | 0.44 and 0.31 |
| lobule model of Debbaut et al. [4] | | 0.95 mmHg imposed | |

The permeability of the lobule layer is a calibration: $k_r=1.6\times10^{-13}$ m² carries 1 mL/min per g at 1 mmHg.
The permeability measured on the human cast, $1.56\times10^{-14}$ m² radially [3], carries a tenth of the
physiological flow at that drop in the lobule model of the same authors (157 of about 1500 mL/min) [4]. Darcy flow
depends on $k\,\Delta P/\mu$ only, so every flow, velocity, shear, transit time and concentration of the lobule
layer is identical to release 0.4.3; its pressures are not.

The choice of 1 mmHg moves little. Over 0.35 mmHg (the optimum) to 2.2 mmHg (the rat by micropuncture) the portal
radii stay within 19% of the cast and the resistance load of a lobar graft with the arterial inflow stays between
0.82 and 0.90.

In the adult biopsy a core of 16.4 mm² holds 21 portal tracts and 7 terminal hepatic venules of 67 µm, and 0.8 ± 0.5
complete portal triads per mm² [5]. The lobule layer uses a central vein of 60 µm.

## 5. The adult liver

Liver of 1500 g, lobular drop imposed:

| $b$ | portal trunk radius (mm) | mean velocity in the trunk (cm/s) | wall shear, trunk to terminal branch (Pa) | volume of the two trees (mL) | $z_R$ of a lobar graft | with the arterial inflow, $\varphi\le1$ |
|---|---|---|---|---|---|---|
| 0.60 | 9.9 | 7 | 0.10 to 2.8 | 201 | 0.997 | 0.86 to 0.90 |
| 0.67 | 8.4 | 10 | 0.17 to 2.4 | 165 | 0.994 | 0.85 to 0.90 |
| 0.75 | 7.0 | 15 | 0.29 to 2.0 | 136 | 0.989 | 0.85 to 0.89 |
| 0.85 | 5.8 | 21 | 0.50 to 1.6 | 117 | 0.976 | 0.84 to 0.88 |
| 0.93 | 5.2 | 27 | 0.73 to 1.2 | 109 | 0.960 | 0.82 to 0.86 |

Measured in adults: portal diameter 10.7 mm by ultrasound and 13.3 mm by CT (pooled) [9], 14.6 mm in the cast [2];
mean velocity 15 ± 4 and 16.5 ± 4.9 cm/s; flow 648 ± 186 and 899 ± 284 mL/min [11], below the 1350 mL/min that 90
mL/min per 100 g gives a liver of 1500 g. In donors the pressure of the right portal vein is 8 ± 1.8 mmHg [20].
The calibration does not move the result of `docs/LIVER.md`, section 6: the six observed resistance loads (0.31 to
0.74) lie below every value in the table.

## 6. Two inflows, measured

`docs/LIVER.md`, section 5, left the arterial flow per lobule of the graft, $\varphi$, unobserved. Two series
measured portal and arterial flow with a transit-time flowmeter in the same tissue, the right lobe in the donor
before retrieval and the graft after reperfusion. `two_inflows_measured` puts their values in the relation of
`liver.py`.

| mL/min per 100 g | Wu 2011 [18], 64 recipients, means | García-Valdecasas 2003 [19], 22 recipients, medians |
|---|---|---|
| donor, portal and arterial | 84.7 and 38.2 | 69.6 and 20.2 |
| graft, portal and arterial | 275.5 and 25.3 | 243.7 and 16.0 |
| flow load $z_F$ | 3.25 | 3.50 |
| arterial share of the donor | 0.31 | 0.22 |
| arterial flow per gram, graft over donor ($\varphi$) | 0.66 | 0.79 |
| total flow per gram, graft over donor | 2.45 | 2.89 |
| portal share of the graft's inflow | 0.92 | 0.94 |
| gradient after reperfusion (mmHg) | 10.6 ± 4.5 | not measured |
| $z_R$ observed | 0.65 | |
| $z_R$ of the model, lobar graft, $b$ 0.6 to 0.93 | 0.82 to 0.86 | 0.86 to 0.90 |
| $z_R$ of the model with all the resistance after the junction | 0.75 | 0.83 |

The arterial flow per gram of the graft falls to 0.66 and 0.79 of the donor's while the portal flow per gram rises
3.3 and 3.5 fold. Wu 2011 is the sixth group of `data/series.tsv` with both loads, read against its own donor
reference. Its resistance load, 0.65, is inside the range of the other five and below the model at any split of the
resistance between the trees.

## 7. Growth with the body

`growth(b, s)` scales the calibrated liver from 150 g to 1500 g at invariant gradient and invariant flow per gram.
The liver can grow by enlarging its lobules or by adding them; $s$ is the share of the growth in volume that is
enlargement (lobule volume $\propto V^{s}$, number of lobules $\propto V^{1-s}$). The sinusoidal bed keeps its
permeability, so the conductance of a lobule grows with its length.

Exponents with liver volume, for $b$ from 0.6 to 0.93:

| | $s=0$ (lobules added) | $s=0.5$ | $s=1$ (lobules enlarged) | measured |
|---|---|---|---|---|
| portal diameter | 0.34 to 0.40 | 0.35 to 0.38 | 0.35 | 0.47 (0.25 to 0.68) [11] |
| mean velocity in the portal trunk | 0.21 to 0.32 | 0.24 to 0.30 | 0.29 | 0.37 (−0.11 to 0.96) [11]; no difference between age groups [12] |
| lobular drop | 0 | 0.33 | 0.67 | not measured |
| velocity in the sinusoid | 0 | 0.17 | 0.33 | not measured |

The portal trunk does not tell how the liver grows: its diameter rises as $V^{0.34}$ to $V^{0.40}$ and its velocity
as $V^{0.21}$ to $V^{0.32}$ whatever $s$. At $b=0.75$ the diameter goes from 6.1 mm to 14.0 mm and the velocity
from 7 or 8 to 15 cm/s. The measurements available as ranges by body weight [11] contain these exponents; the pooled
portal diameter is 6.6 mm in children and 10.7 mm in adults [9]. In absolute value the ultrasound diameters sit on
the curve of $b=0.93$ and below that of $b=0.75$, which reproduces the cast:

| liver volume (mL) [10] | 155 | 334 | 551 | 784 | 1269 |
|---|---|---|---|---|---|
| portal diameter by ultrasound (mm) [9,11] | 3 to 4 | 4 to 7 | 4 to 8 | 6 to 9 | 10.7 |
| model, $b=0.75$ | 6.7 | 8.7 | 10.2 | 11.5 | 13.5 |
| model, $b=0.93$ | 4.9 | 6.4 | 7.5 | 8.5 | 9.9 |

The cast and the living vein give the two ends of the interval for the portal tree, 0.73 to 0.93. The standard
liver volume is that of Urata, $706.2\times$ body surface area $+\,2.4$ mL, from 96 subjects aged one month to 27
years, in whom liver volume over body weight falls with age until about 16 years [10].

What $s$ changes is the lobule: if a newborn had the lobules of an adult in number and they only enlarged ($s=1$),
its lobular drop would be 0.22 mmHg and its sinusoidal velocity 0.46 of the adult's; if lobules were added ($s=0$)
both would be the adult's. In the mouse the number of lobules doubles after birth and stops rising in the second
week, and lobule size is final at four weeks [13]; in the rat the liver lobe grows 12-fold while lobule area grows
more than 5-fold and the number of lobules at the surface by about 30% [14], which is $s$ close to one. No
measurement of lobule size or number against age in people was found.

## 8. Regrowth of a graft

After partial hepatectomy the rat liver regrows by enlarging its lobules; new lobules do not form [14,15]. In adult
recipients the graft reaches 1.6 to 1.9 times its implanted volume in two weeks, about the standard volume of the
recipient:

| volume over implanted volume | recipients | 2 weeks | 1 month | 2 to 3 months | 6 months |
|---|---|---|---|---|---|
| Yagi 2005, portal pressure below 20 mmHg [21] | 15 | 1.61 | 1.54 | 1.49 | 1.52 |
| Yagi 2005, portal pressure of 20 mmHg or more [21] | 17 | 1.94 | 1.99 | 2.02 | 1.47 |
| Eguchi 2003, cirrhosis [22] | 9 | 1.66 | | | |
| Byun 2016 [17] | 32 | 1.7 | | | |
| Zhang 2022 [16] | 62 | 1.58 | 1.50 | 1.44 | 1.31 |
| Jiang 2009 [23] | 18 | | 2.19 | | |
| García-Valdecasas 2003 [19] | 22 | | | 1.81 | |
| Olthoff 2015 [24] | 308 | | | 1.56 | |
| Wu 2011 [18] | 64 | | | 2.09 | |

The portal inflow does not stay at its value of reperfusion. It falls from 1700 to 1199 mL/min in ten days [18];
the portal velocity is 106, 79, 57 and 59 cm/s on days 1, 3, 5 and 7, against 24 to 32 cm/s in the donors [25], and
its peak falls from 95 cm/s one hour after reperfusion to 59 cm/s at three months, against 48 cm/s in the donors
[20].

`regeneration(rho, inflow=...)` gives the loads of a graft whose volume has grown $\rho$-fold by lobule enlargement
with $h$ times the inflow of reperfusion. Lengths grow as $\rho^{1/3}$. The lobule keeps its permeability, so its
resistance falls as $\rho^{-1/3}$. The resistance of the trees grows as $\rho^{1/3}$ if their calibres do not change
and falls as $1/\rho$ if the calibres grow with the lengths. With the lobule holding one fifth of the drop,

$$\frac{z_F}{z_{F,0}}=\frac h\rho,\qquad \frac{z_P}{z_{P,0}}=h\left(0.8\,\rho^{\pm}+0.2\,\rho^{-1/3}\right),\qquad \rho^{\pm}=\rho^{1/3}\ \text{or}\ \rho^{-1}.$$

The first week, with the measured growth and the measured inflow (`first_week`):

| at one week, relative to reperfusion | |
|---|---|
| graft volume [21] | ×1.57 |
| portal inflow, day 10 [18] | ×0.70 |
| flow load | ×0.45 |
| model, portal pressure with calibres unchanged | ×0.78 |
| model, portal pressure with calibres in proportion | ×0.48 |
| measured portal pressure, 12 hours to day 7 [21] | 24.2 to 22.0 mmHg and 16.2 to 14.8 mmHg, ×0.91 |
| the same as a gradient above a caval pressure of 8 mmHg | ×0.83 to ×0.86 |
| resistance of the graft implied by the measurements | ×1.18 to ×1.30 |
| resistance of the graft in the model | ×1.10 unchanged, ×0.68 in proportion |

The flow load falls by more than half in a week. The portal pressure falls by a tenth, which excludes vessels that
widen with the graft and sits at or above the bound of unchanged calibres: over the first week the resistance of the
graft rises while its volume grows by half.

At constant inflow, for reference:

| months | volume [16] | flow load | pressure load, calibres unchanged | pressure load, calibres in proportion |
|---|---|---|---|---|
| 0.5 | ×1.58 | ×0.63 | ×1.10 | ×0.68 |
| 1 | ×1.50 | ×0.67 | ×1.09 | ×0.71 |
| 3 | ×1.44 | ×0.69 | ×1.08 | ×0.73 |
| 6 | ×1.31 | ×0.76 | ×1.06 | ×0.79 |

A graft of one third that has regrown threefold with calibres in proportion returns to the donor's flow per gram
with a pressure load of 1.22 and a sinusoidal velocity 1.44 times the donor's, because it carries the flow of a
whole liver through a third of the lobules.

## 9. Limits

The anatomy is one cast of one liver, made after machine perfusion, with estimated vessel counts beyond the sixth
generation. Flows are the references of the manuscript, viscosity is constant and vessels are rigid tubes. The share
of the pressure drop by compartment has been measured in the rat and the dog, not in people. The paediatric portal
measurements are ranges by body weight. The first week joins three cohorts: the inflow of one, the volume and the
portal pressure of another, read from its figure as means with standard deviations of 3 to 7 mmHg and without the
caval pressure. The growth of human lobules with age is not known to this note.

## Sources

1. Debbaut C, et al. J Anat 2014;224:509-517 (human liver cast, Table 1: radius, length and estimated number of vessels of generations 5 to 13).
2. Debbaut C, et al. From vascular corrosion cast to electrical analog model for the study of human liver hemodynamics and perfusion. IEEE Trans Biomed Eng 2011;58:25-35 (Tables II to IV: generations 1 to 6 of the same cast).
3. Debbaut C, et al. Perfusion characteristics of the human hepatic microcirculation based on three-dimensional reconstructions and computational fluid dynamic analysis. J Biomech Eng 2012;134:011003.
4. Debbaut C, et al. A 3D porous media liver lobule model: the importance of vascular septa and anisotropic permeability for homogeneous perfusion. Comput Methods Biomech Biomed Engin 2014;17:1295-1310.
5. Crawford AR, Lin XZ, Crawford JM. The normal adult human liver biopsy: a quantitative reference standard. Hepatology 1998;28:323-331.
6. Lautt WW. Hepatic Circulation: Physiology and Pathophysiology, chapter 2 (NCBI Bookshelf NBK53069).
7. Shibayama Y. Kanzo 1975;16:145-152 (normal rat: portal vein 110, terminal portal venule 68, central vein 28, vena cava 20 mmH2O).
8. Bohlen HG, Maass-Moreno R, Rothe CF. Hepatic venular pressures of rats, dogs, and rabbits. Am J Physiol 1991;261:G539-G547.
9. Chau P, Yoon JS, Moses D, Pather N. A systematic review and meta-analysis of portal vein morphometry in pediatric and adult populations. Eur J Radiol 2023;168:111016.
10. Urata K, et al. Calculation of child and adult standard liver volume for liver transplantation. Hepatology 1995;21:1317-1321; used in `urata_volume`.
11. López Barrio AM, et al. An Esp Pediatr 1996;44:45-49 (82 healthy children; adult values of Moriyasu and Ohnishi as quoted there).
12. Ohtsuka Y. J Jpn Soc Pediatr Surg 1994;30:1267-1277 (106 normal children, mean portal velocity 15.4 ± 3.5 cm/s).
13. Azizoglu DB, et al. Vascularization of neonatal liver lobules presages adult liver size. Nat Commun 2025;16:9989.
14. Papp V. Architectural changes in the rat liver during ontogenesis, regeneration and oncogenesis. Doctoral thesis, Semmelweis University, 2013.
15. Tsomaia K, et al. Liver structural transformation after partial hepatectomy and repeated partial hepatectomy in rats. World J Gastroenterol 2020;26:3899-3916.
16. Zhang Y, et al. Comparison of liver regeneration between donors and recipients after adult right lobe living-donor liver transplantation. Quant Imaging Med Surg 2022;12:3184-3192.
17. Byun SH, Yang HS, Kim JH. Liver graft hyperperfusion in the early postoperative period promotes hepatic regeneration 2 weeks after living donor liver transplantation. Medicine 2016;95:e5404.
18. Wu TJ, et al. Impact of portal venous hemodynamics on indices of liver function and graft regeneration after right lobe living donor liver transplantation. Liver Transpl 2011;17:1035-1045.
19. García-Valdecasas JC, et al. Changes in portal vein flow after adult living-donor liver transplantation: does it influence postoperative liver function? Liver Transpl 2003;9:564-569.
20. Gondolesi GE, et al. Venous hemodynamics in living donor right lobe liver transplantation. Liver Transpl 2002;8:809-813.
21. Yagi S, et al. Impact of portal venous pressure on regeneration and graft damage after living-donor liver transplantation. Liver Transpl 2005;11:68-75 (pressures of Fig. 1 read from the figure).
22. Eguchi S, et al. Relationship between portal venous flow and liver regeneration in patients after living donor right-lobe liver transplantation. Liver Transpl 2003;9:547-551.
23. Jiang SM, et al. Role of splanchnic hemodynamics in liver regeneration after living donor liver transplantation. Liver Transpl 2009;15:1043-1049.
24. Olthoff KM, et al. Liver regeneration after living donor transplantation: adult-to-adult living donor liver transplantation cohort study. Liver Transpl 2015;21:79-88.
25. Sugimoto H, et al. Normal hepatic hemodynamics during early postoperative period in recipients with adult live donor liver transplantation. Transplant Proc 2007;39:160-165.
