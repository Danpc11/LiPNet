# The network model applied to the liver

This note states what the cost-constrained network of `src/model/` implies for a human liver and for a partial
graft, and what the extraction table and the fits support, with nothing added to the model. Every number is
produced by

```bash
python src/model/liver.py     # -> results/model_liver.json
python src/audit.py           # -> results/audit.json
python tests/test_liver.py
```

## 1. Wall shear stress is not the scale of the design

The optimum minimises $\sum_e f_e^2/w_e$ at fixed $\sum_e m_e^b$, with $w_e=r_e^4/l_e$ and $m_e=r_e^2 l_e$. Its
stationarity condition equalises $D_e/c_e$ over the edges, and that common value is $\tau_0^2=D^*/C$. The Poiseuille
wall shear of an edge is $\tau_e=f_e/r_e^3$, and

$$\tau_e=\tau_0\,m_e^{(b-1)/2}.$$

The two coincide only at $b=1$ (Murray's law). For $b<1$ shear rises from the trunk to the sinusoid: by a factor of
7 at $b=0.75$ and of 73 at $b=0.5$ in a liver of $1.6\times10^6$ lobules. $\tau_0$ is a scale of the whole design,
not the shear held at any wall.

The allometric closure (invariant portocaval pressure drop across body mass) therefore has three exponents, not one:

| Quantity | Exponent with body mass, over $b\in[0.5,1]$ and the four scaling inputs |
|---|---|
| $\tau_0=\sqrt{D^*/C}$ | −0.077 to +0.071 |
| wall shear in the sinusoid | −0.077 to −0.003 |
| wall shear in the trunk | −0.298 to −0.033 |
| vascular mass | 0.89 to 1.08 |

Sinusoidal shear is within 0.08 of invariant in every case; trunk shear is not for $b<1$. The closure is compared
with one observation, the scaling of hepatic blood volume ($M^{0.86}$, 0.80 to 0.92). Five of the 28 combinations
fall inside it, and they are the same family: portal flow $\propto M^{0.77}$, lobule size $\propto M^{0.10}$,
$b\ge0.75$.

## 2. The exponent of a human liver

The optimum fixes the diameter ratio between successive generations,

$$\frac{r_{i+1}}{r_i}=n^{-1/(2+b)}\,\lambda^{(1-b)/(2(2+b))},$$

with $n$ the branching number and $\lambda$ the length ratio. Inverting it for the human casts already quoted in
`network.py`:

| Tree | diameter ratio | length ratio | $n$ | model at $b=0.5$ | model at $b=1$ | $b$ that reproduces it |
|---|---|---|---|---|---|---|
| portal vein | 0.70 | 0.72 | 2.80 | 0.641 | 0.709 | 0.92 |
| hepatic vein | 0.59 | 0.64 | 3.22 | 0.599 | 0.677 | 0.45 |
| hepatic artery | 0.74 | 0.66 | 2.76 | 0.639 | 0.713 | 1.22 |

One exponent does not describe both venous trees. A model with one exponent per tree is not closed: the costs
$m^b$ of the two trees then carry different units, and their joint optimum depends on one further number, the
budget of one tree relative to the other, which the model does not supply. The exponent of the hepatic vein is
reported and not used: by Poiseuille the same radii hold more than the gradient of the whole liver
([HUMAN.md](HUMAN.md), section 3). On the eleven measured generations of the cast the portal vein gives 0.73 and the
hepatic artery 1.04. Results are reported across $b$, and where they depend on the split between the
trees, over every split.

## 3. Where the pressure drop sits in a human liver

Closed-form canopy, $3^{13}=1.6\times10^6$ lobules, $d=3$:

| $b$ | portal tree | sinusoids | hepatic venous tree | two trunks | exponent of $D^*$ with $N$ at this size | asymptote $(2-b)/b$ |
|---|---|---|---|---|---|---|
| 1 | 0.483 | 0.034 | 0.483 | 0.069 | 1.20 | 1.00 |
| 0.92 | 0.476 | 0.048 | 0.476 | 0.045 | 1.30 | 1.17 |
| 0.75 | 0.456 | 0.088 | 0.456 | 0.013 | 1.71 | 1.67 |
| 0.5 | 0.424 | 0.151 | 0.424 | 0.001 | 3.00 | 3.00 |

The sinusoids carry 3 to 15% of the drop. The scaling $D^*\propto N^{(2-b)/b}$ is the large-$N$ limit for $b<1$; at
$b=1$ the local exponent is 1.20 at human size, and on the explicit graphs (37 to 1,478 lobules) it is 1.36 to 3.30.
What the explicit graphs verify is the closed-form canopy evaluated at the same size.

## 4. A partial graft is a sub-network

A graft is the sub-network below one first-order branch. Its vessels keep their conductance, so with $z_F$ times the
donor flow in every lobule the drop across it is $z_F$ times the donor drop, less the trunks it left behind:

$$z_R(\text{graft})=1-\frac{\Delta P_{\text{removed generations}}}{\Delta P_{\text{donor}}}.$$

| | $z_R$ of a lobar graft |
|---|---|
| canopy, $b=1$ to $0.5$ | 0.93 to 1.00 |
| explicit 3D graphs (160 and 293 lobules), one first-order portal branch | 0.90 to 0.96 |

On the explicit graphs the graft keeps the hepatic venous edges that drain it at the calibre they had in the whole
liver, which lowers the outlet term to 0.78 to 0.91 of its donor value.

## 5. Two inflows

The hepatic artery joins the portal blood at the sinusoid. The portal tree carries portal flow; the sinusoids and
the hepatic veins carry portal plus arterial flow. The clinical flow load is a portal load, so with $a$ the
arterial-to-portal flow ratio of the donor (one third for an arterial share of 25%), $\varphi$ the arterial flow per
lobule of the graft relative to the donor, and $s_{in}$ the share of the portal tree in the donor drop,

$$z_P=s_{in}z_F+(1-s_{in})\frac{z_F+\varphi a}{1+a}.$$

| arterial flow per lobule | $z_R$ of a lobar graft at $z_F$ 3.3 to 4.1, $b=1$ to $0.5$ | for any split between the trees |
|---|---|---|
| none ($\varphi=0$) | 0.81 to 0.86 | 0.75 to 1 |
| donor value ($\varphi=1$) | 0.84 to 0.90 | 0.81 to 1 |
| scaling with the portal vein ($\varphi=z_F$) | 0.93 to 1.00 | the network value |

The series of the table report no arterial flow. Two others measured it in the donor and in the graft: $\varphi$ is 0.66
and 0.79, and the model then gives 0.82 to 0.90 against an observed 0.65 ([HUMAN.md](HUMAN.md), section 6).

## 6. The model against the six paired groups

| | $z_R$ |
|---|---|
| network, lobar graft | 0.93 to 1.00 |
| explicit graphs | 0.90 to 0.96 |
| with the arterial inflow, $\varphi\le1$ | 0.81 to 0.90 |
| with the shear feedback of the lobule layer (`src/fem/network_link.py`) | 0.91 to 0.92 |
| observed | 0.31, 0.54, 0.55, 0.65, 0.72, 0.74 |

In physical units, with the lobular drop shared with the lobule layer, a lobar graft is at 0.96 to 1.00 and at 0.82 to
0.90 with the arterial inflow ([HUMAN.md](HUMAN.md), section 5). At the reference values all six observed values lie below the lowest value the model reaches; with the donor references varied together this holds for every group in 34 of 60 combinations, and in none with a normal gradient of 4 mmHg (THEORY.md, P2). The model has the sign of the departure and
not its size. Widening the outlet anastomosis in `whole_graft.py` lowers $z_R$ to 0.95 at most, because that
anastomosis holds 5% of the resistance.

## 7. What the rigid network leaves out: the passive relations of the normal liver

The network is rigid: its calibres are those of the optimum and do not change with pressure. `src/model/distension.py`
asks how much of the distance to the observations each mechanism it leaves out can carry. Nothing is fitted to the six
groups; every number is the model's own or a published one.

The rigid network with its two inflows gives 0.81 on average in its lowest configuration and 0.90 in its highest, and
the groups report 0.58: the resistance would have to fall by a further 28 to 35%. (No figure is computed against the
midpoint of that range, which is not a configuration of the network.)

| Mechanism | Acts on | Lowest $z_R$ it can give, relative to the rigid value |
|---|---|---|
| Dilation or recruitment of sinusoids, unlimited | the 3 to 15% of the drop the sinusoids carry ($b$ from 1 to 0.5) | 0.85 to 0.97 |
| Outlet anastomosis, unlimited widening | the 5% of the resistance it holds | 0.95 |
| Uniform widening of a share $s$ of the bed | $s$ | closes the distance only if $s\geq0.28$ to $0.35$; with the whole bed, a radius 8 to 11% larger |

So neither the sinusoids nor the outlet can account for the distance: it has to come from vessels that hold more than
a quarter of the resistance, which in this network are the two trees.

### Three published relations

The venous resistance of the normal liver is passively distensible. Three relations measured on the whole livers of
cats say by how much; all were read in the original papers.

- **Cube law.** Resistance is proportional to $1/P_d^3$ (index of contractility, $IC=R\,P_d^3$), with $P_d$ the mean
  of the pressures on either side of the resistance [L3, L2], here CVP + $G/2$ for a gradient $G$. It was derived in
  isolated perfused livers over 10 to 50 mL/min/kg, a third to 1.7 times the control flow, with outflow pressures up to
  15 mmHg, and confirmed in vivo, separately for the presinusoidal and the postsinusoidal resistance sites, whose
  indices add in series [L3].
- **Compliance.** Hepatic blood volume is linear in intrahepatic pressure. Raising the hepatic venous pressure by
  9.4 mmHg added 32 mL/100 g to 27 [L4]; "changing intrahepatic pressure by 8 mmHg will double hepatic blood volume"
  [L2]; a compliance of 2.6 to 3.1 mL/mmHg per 100 g on a volume of 34 mL/100 g doubles it in 11 to 13 mmHg [L8]. The
  doubling pressure is taken over 8 to 13 mmHg. With the length of the vessels fixed and Poiseuille flow, resistance
  goes as volume$^{-2}$: this takes the resistance vessels to distend as the capacitance vessels do, an assumption
  and an upper bound.
- **Pressure-flow.** Intrahepatic pressure, the mean of the portal and the lobar venous pressures, is linear in total
  hepatic flow, measured up to 240 mL/min/100 g of portal flow, about 2.2 times the control total flow; nonlinear fits
  were worse [L8]. The intercept is 2.2 mmHg at zero flow and the slope 0.046 mmHg per mL/min/100 g, of which 0.015 is
  the outflow cannula raising the caval pressure with flow, so 0.030 relative to the caval pressure, at a control flow of
  139. On the portocaval gradient used here, which exceeds the intrahepatic-to-caval difference by half the 1 to 2 mmHg
  between portal and lobar veins (taken proportional to flow), a share $f$ = 0.30 to 0.35 of the control gradient does
  not scale with flow, and $G/g_0=f+(1-f)\,q$ at a flow load $q$. The slope with the cannula (which would give 0.26) is
  not a gradient and is not used. This is the relation measured at the highest flows, and it stops distending:
  resistance tends to the slope as flow rises.

The first two keep distending at the flows of a graft, the third does not. Their authors note that in the normal liver
the changes of resistance are notable below 10 mmHg and minor above [L7], and a graft is at 14 to 22 mmHg.

For a group with flow load $z_F$ and central venous pressure CVP each relation is solved for the gradient $G$ from
$G = 5\,z_F\,z_R^{\text{rigid}}\,d(G)$, with the donor at 5 mmHg; the pressure-flow relation is applied to the
effective load $z_F\,z_R^{\text{rigid}}$. The inputs are the flow load, the rigid network and the central venous
pressure; the portal pressure of the recipient is the output and is not used. The rigid value runs from 0.81 (with the
two inflows) to 1.00 (the network alone).

### The relations against measurements they were not derived from

Doubling the flow of a donor liver raises its gradient by 1.9 mmHg with the cube law and 1.7 to 2.3 mmHg with the compliance; the source
gives about 2 mmHg [L6]. The pressure-flow relation gives 3.3 to 3.5 mmHg.

Ricciardi et al. [L5] measured portal pressure and resistance at two flows in pig livers in situ, isolated on a bypass
circuit, and after 2 h of cold ischaemia. Each relation is anchored at the lower flow of each state, with the gradient
and caval pressure that the printed pressure and resistance imply, and gives the resistance and the portal pressure at
the higher flow (`published_check`). The caval pressure implied at the two flows agrees within 0.3 mmHg in every state.

| Pig liver | Flow raised | Resistance, measured (± SE) | Cube law | Compliance | Pressure-flow |
|---|---:|---:|---:|---:|---:|
| in situ (9 livers) | 2.75-fold | 0.54 ± 0.07 | 0.48 | 0.53 to 0.59 | 0.78 to 0.81 |
| isolated, perfused (7 livers) | 2.67-fold | 0.59 ± 0.05 | 0.48 | 0.57 to 0.63 | 0.78 to 0.81 |
| after 2 h of cold ischaemia (6 livers) | 2.67-fold | 0.59 ± 0.08 | 0.48 | 0.50 to 0.55 | 0.78 to 0.81 |

| Pig liver | Portal pressure at the higher flow, measured (± SE) | Error, cube law | Error, compliance | Error, pressure-flow |
|---|---:|---:|---:|---:|
| in situ | 9.6 ± 0.2 mmHg | −0.8 | +0.2 to +1.1 | +4.5 to +5.0 |
| isolated, perfused | 8.1 ± 0.5 mmHg | −1.6 | −0.4 to +0.4 | +2.6 to +3.0 |
| after 2 h of cold ischaemia | 15.8 ± 1.3 mmHg | −2.8 | −2.3 to −1.0 | +5.5 to +6.3 |

In a normal pig liver the compliance relation holds: it predicts the portal pressure at the higher flow within 1.2 mmHg.
The cube law gives too little resistance and the pressure-flow relation too much. After cold ischaemia the resistance
was twice as high at every flow and fell with flow by the same factor, but at its higher pressure the relations
expected more distension than occurred: the ischaemic liver distended less than a normal liver at the same pressure.

### The six groups

| Group | $z_F$ | CVP (mmHg) | Rigid network | Cube law | Compliance | Pressure-flow | Observed (95% sampling interval) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Chan 2011 | 3.93 | 8.0 | 0.81 to 1.00 | 0.47 to 0.52 | 0.42 to 0.52 | 0.62 to 0.78 | 0.31 (0.25 to 0.37) |
| Wang 2014, splenectomy | 4.12 | 3 to 9, assumed | 0.81 to 1.00 | 0.39 to 0.52 | 0.40 to 0.50 | 0.61 to 0.77 | 0.54 (0.48 to 0.60) |
| Yagi 2006 | 3.28 | 5.2 | 0.81 to 1.00 | 0.48 to 0.53 | 0.46 to 0.57 | 0.64 to 0.79 | 0.55 (0.44 to 0.69) |
| Wu 2011 | 3.25 | 4.4 | 0.81 to 1.00 | 0.47 to 0.52 | 0.47 to 0.57 | 0.64 to 0.79 | 0.65 (0.57 to 0.75) |
| Wang 2014, no splenectomy | 3.74 | 3 to 9, assumed | 0.81 to 1.00 | 0.41 to 0.54 | 0.43 to 0.53 | 0.62 to 0.78 | 0.72 (0.64 to 0.80) |
| Ishizaki 2012 | 3.34 | 9.1 | 0.81 to 1.00 | 0.52 to 0.58 | 0.46 to 0.56 | 0.63 to 0.79 | 0.74 (0.60 to 0.92) |

At these flow loads the relations that keep distending give 0.39 to 0.58 and the one that stops gives 0.61 to 0.79.
Each group is placed against four regions with one rule, the point estimate and the 95% sampling interval reported
apart (`bracket`, `_place`):

| Region | Points inside | Points below | Intervals reaching it | Intervals below |
|---|---:|---:|---:|---:|
| rigid network | 0 | 6 | 1 | 5 |
| keeps distending | 1 | 1 | 3 | 1 |
| stops distending | 3 | 3 | 4 | 2 |
| all relations | 5 | 1 | 5 | 1 |

All six points lie below the rigid network, and five intervals (that of Ishizaki 2012 reaches it). One point lies inside
the relations that keep distending (Yagi 2006), three inside the one that stops (Wu 2011, Ishizaki 2012, Wang 2014
without splenectomy), and Chan 2011 lies below every region: its gradient was 1 mmHg above the donor's at four times the
flow. This is a description: a sampling interval that reaches a region is not a probability that its mechanism acts.

### The three relations are one composition (`distension.consistency`)

The three relations use one definition of pressure (the portocaval gradient over the donor's) and of flow (the portal
flow load), and each acts on the gradient ratio the rigid network gives, $q=z_F z_R^{\text{rigid}}$. Checked in code:
with the whole donor network at $z_F=1$ every relation returns $z_P=z_F=z_R=1$ at every CVP and parameter; two networks
with the same $q$ give the same gradient, so architecture and arterial inflow enter once; over $z_F$ from 0.25 to 20,
$z_P$ rises and $z_R$ falls with no jump. At the highest flow the local exponent of $z_P$ on $z_F$ is 0.32 to 0.38 (cube law) and 0.37 to 0.41 (compliance), still falling, and 0.97 to 0.98 (pressure-flow), whose resistance tends to $(1-f)$ times the rigid one.

### Which vessels distend (`src/model/compartments.py`)

The whole-liver form multiplies the resistance of the rigid network by one factor, as if every vessel distended as the
whole cat liver did. `compartments.py` applies the relations vessel by vessel instead. The cube law was validated for
each resistance site with its own distending pressure [L3], so each generation $i$ of the network keeps its own index,
$R_i = IC_i/P_{d,i}^3$, with $P_{d,i}$ the mean pressure in that generation as the pressure profile of the network gives
it; with the compliance relation each generation changes its volume with its own change of pressure. The vessels that do
not distend keep their rigid resistance, $R/R_0=(1-s)+s\,d$, and $s$ is fixed by the network, not by the groups:

| Distending vessels | $s$ | Relations that keep distending, over the six groups | Groups within all relations | Below |
|---|---:|---:|---:|---|
| whole liver, published form | 1 | 0.39 to 0.58 | 5 | Chan 2011 |
| every vessel at its own pressure | 1 | 0.39 to 0.65 | 5 | Chan 2011 |
| pre- and postsinusoidal sites (all but the sinusoids) | 0.80 to 0.97 | 0.39 to 0.70 | 5 | Chan 2011 |
| postsinusoidal site, the hepatic venous tree (basal resistance site of the cat) | 0.37 to 0.48 | 0.67 to 0.92 | 4 | Chan 2011; Wang 2014, splenectomy |
| presinusoidal site, the portal tree | 0.42 to 0.48 | 0.46 to 0.76 | 5 | Chan 2011 |

The rigid limit is the network of P2 exactly. A higher flow raises the pressure of the vessels upstream of the sinusoids
by most of the rise of the gradient and that of the hepatic veins, which sit close to the caval pressure, by little: the
relations that keep distending reach the lower observations only if the portal tree distends.

The distending pressure of the whole-liver form was also replaced by CVP − $P_{ext}$ + $\lambda G$, $\lambda$ from 0.25
to 1 and an external pressure of 0 or 2 mmHg: the span of the relations that keep distending moves between 0.34 and 0.77,
and in every variant five groups lie within all the relations and Chan 2011 below them. The placement rests on which
vessels distend, not on the choice of distending pressure.

### Every uncertainty together (`distension.joint`)

Measurement (the two group means), references (donor flow, normal gradient, assumed CVP), network (eight networks,
arterial flow per lobule from none to the donor's) and relations (doubling pressure, flow-independent share) drawn
together, 20,000 draws per group, the observed ratio and each model from the same draw. The ranges are published or
plausible ranges drawn uniformly, not measured distributions, so the shares of draws below are a sensitivity analysis
over the joint range, one group at a time, and not probabilities that a mechanism acts; the deterministic ranges above
remain the primary statement.

| Group | Observed $z_R$, median (95% range) | Draws below the rigid network | Below every relation | Within | Above every relation |
|---|---:|---:|---:|---:|---:|
| Chan 2011 | 0.31 (0.23 to 0.42) | >99% | >99% | <1% | <1% |
| Wang 2014, splenectomy | 0.49 (0.31 to 0.74) | >99% | 32% | 58% | 9% |
| Yagi 2006 | 0.55 (0.39 to 0.77) | >99% | 23% | 66% | 11% |
| Wu 2011 | 0.65 (0.49 to 0.87) | 97% | <1% | 61% | 39% |
| Wang 2014, no splenectomy | 0.66 (0.45 to 0.97) | 90% | 2% | 50% | 47% |
| Ishizaki 2012 | 0.74 (0.53 to 1.03) | 80% | <1% | 31% | 69% |

### Level and order (`distension.heterogeneity`, `distension.wang_contrast`)

The relations depend on the flow load and CVP of each group and on parameters common to all, so they set the level of
$z_R$, not its order. At central parameters the standard deviation of $\log z_R$ over the six groups is
0.016 to 0.075 against 0.325 observed, and the rank correlation with the observations 0.49 (P = 0.33).
The one contrast inside a cohort, the two arms of Wang 2014, solved from the flow loads only (4.12 and 3.74) with
CVP, network and parameters common to both arms in each draw:

| Model | Ratio of $z_R$, splenectomy to none (95% range) | Share of the log contrast |
|---|---:|---:|
| observed | 0.74 (0.61 to 0.88) | 100% |
| rigid network | 1.00 (0.99 to 1.00) | 0% |
| cube law | 0.95 (0.88 to 1.02) | 18% |
| compliance | 0.95 (0.88 to 1.02) | 18% |
| pressure-flow | 0.99 (0.97 to 1.00) | 4% |

From the flow loads alone the relations give the direction and at most a fifth of the contrast. The rest is what the flow
load of the arms does not carry; the arms were not randomised, so these data cannot say whether it is the splenectomy or
the recipients chosen for it.

What this supports and what it does not. The passive relations of the normal liver, extrapolated to the flow of a graft,
span the resistance loads the cohorts report, and the rigid network does not; the lower ones only if the portal tree
distends. That bounds the mechanism and does not identify it. The relations come from normal cat livers, mostly below
10 mmHg and at a third of the flow of a graft or less, and disagree with each other there by as much as the groups
differ. A graft differs by species, ischaemia and reperfusion, endothelial injury and vascular tone, a partial graft is
not a whole liver, and the pig livers after cold ischaemia distended less than the relations give. The relations are
read as the largest plausible passive contribution, not as a prediction for a graft. Which relation a graft follows
needs flow and pressure in the same graft at two or more flows.

[L2] Lautt WW, Legare DJ. Passive autoregulation of portal venous pressure: distensible hepatic resistance. *Am J Physiol* 1992;263:G702–G708.
[L3] Lautt WW, Greenway CV, Legare DJ. Index of contractility: quantitative analysis of hepatic venous distensibility. *Am J Physiol* 1991;260:G325–G332.
[L4] Lautt WW, Greenway CV. Hepatic venous compliance and role of liver as a blood reservoir. *Am J Physiol* 1976;231:292–295.
[L5] Ricciardi R, Foley DP, Quarfordt SH, et al. Alterations in intrahepatic hemodynamics of the harvested porcine liver. *J Gastrointest Surg* 2001;5:490–498.
[L6] Lautt WW, Ming Z. Hepatic hemodynamics in portal hypertension. In: Sanyal AJ, Shah VH, eds. *Portal Hypertension: Pathobiology, Evaluation, and Treatment*. Totowa: Humana Press; 2005. p. 85–97.
[L7] Inglés AC, Legare DJ, Lautt WW. Distensibility of portacaval shunts in portal hypertensive cats: index of contractility model. *Can J Physiol Pharmacol* 1994;72:687–692.
[L8] Greenway CV, Seaman KL, Innes IR. Norepinephrine on venous compliance and unstressed volume in cat liver. *Am J Physiol* 1985;248:H468–H476.

## 8. What the table supports

- Arithmetic of `data/series.tsv`: no disagreement between percentages and counts, or between stored gradients and
  the pressures they come from.
- The donor flow reference is measured in four series (81, 85, 91 and 118 mL/min/100 g) and set to 90 elsewhere. No
  series reports a donor gradient; the normal gradient of 5 mmHg is assumed for all 22 groups with a pressure load.
  The table does not record when the flow was measured relative to the pressure.
- Sampling intervals of $z_R$ with the dispersions the papers report (largest reported value where none is given)
  stay below one: 0.25 to 0.37, 0.44 to 0.69, 0.48 to 0.60, 0.64 to 0.80, 0.60 to 0.92 and, for Wu 2011, 0.57 to
  0.75.
- The reference values decide the level. A group reaches $z_R=1$ if the donor gradient is at or below 1.5, 2.7, 2.7,
  3.6, 3.7 and 3.3 mmHg respectively; at the largest measured donor reference (118) the two highest groups are at
  0.94 and 0.97.
- Without any reference: the flow load spans 1.27-fold across the six groups and the pressure load 2.23-fold;
  their rank correlation is +0.09; 12% of the variance of $\log z_P$ is in $\log z_F$. In Wang 2014 the arm with
  10% more flow has 12% less portal pressure.

## 9. What the outcome model supports

| stratum | groups defined by | span of the predictor | slope on its own (SE) | events |
|---|---|---|---|---|
| Osman 2017 | final portal pressure | 5.0 mmHg | +2.00 (1.10) | 7/76 |
| Wang 2014 | splenectomy | 2.3 mmHg | +1.78 (0.76) | 39/276 |
| Uemura 2016 | graft size | 0.37 mmHg | +4.51 (7.35) | 30/221 |
| Ogura 2010 | final portal pressure | 3.7 mmHg | +2.52 (0.72) | 20/129 |
| Yagi 2005 | portal pressure, days 1 to 3 | 8.3 mmHg | +1.37 (0.69) | 8/32 |

Uemura 2016 separates its groups by 0.37 mmHg against a within-group standard deviation of 3 to 4 mmHg and carries
no information on the slope: the primary estimate is +1.55 (+0.33 to +2.79) with it and +1.53 (+0.30 to +2.77)
without it, on 46 events in 352 recipients from two series. The same model fitted to the flow load gives +0.84
(−0.59 to +2.36) on 48 events, mostly from other cohorts; the two slopes are not compared, because the two fits
share one series and a joint model would need cohorts with both loads and an outcome. The network
model enters this analysis only through the normalisation of the predictor: it contains no relation between load
and outcome.
