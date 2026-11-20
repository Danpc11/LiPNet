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

## 7. What the rigid network leaves out: a passively distensible bed

The network is rigid: its calibres are those of the optimum and do not change with pressure. `src/model/distension.py`
asks how much of the distance to the observations each mechanism it leaves out can carry. Nothing is fitted to the six
groups; every number is the model's own or a published one.

The rigid network with its two inflows gives 0.85 on average and the groups report 0.58: the resistance would have to
fall by a further 32%.

| Mechanism | Acts on | Lowest $z_R$ it can give, relative to the rigid value |
|---|---|---|
| Dilation or recruitment of sinusoids, unlimited | the 3 to 15% of the drop the sinusoids carry ($b$ from 1 to 0.5) | 0.85 to 0.97 |
| Outlet anastomosis, unlimited widening | the 5% of the resistance it holds | 0.95 |
| Uniform widening of a share $s$ of the bed | $s$ | closes the distance only if $s\geq0.32$; with the whole bed, a radius 10% larger |

So neither the sinusoids nor the outlet can account for the distance: it has to come from vessels that hold at least a
third of the resistance, which in this network are the two trees.

### Two published relations

The venous resistance of the normal liver is passively distensible. Two relations measured on the livers of cats say by
how much.

- **Cube law.** Resistance is proportional to $1/P_d^3$ (index of contractility, $IC=R\,P_d^3$). The distending pressure
  $P_d$ is the mean of the pressures above and below the resistance site [L2, L3, L6, L7], here CVP + $G/2$ for a
  gradient $G$.
- **Compliance.** Hepatic blood volume doubles for about 8 mmHg of intrahepatic pressure, which the portal venous
  pressure estimates, linearly over the physiological range [L1, L4, L8]. With the length of the vessels fixed and
  Poiseuille flow, resistance goes as volume$^{-2}$. This takes the resistance vessels to distend as the capacitance
  vessels do, which is an assumption and an upper bound.

For a group with flow load $z_F$ and central venous pressure CVP the gradient $G$ of a bed that follows the relation $d$
solves

$$G = 5\,z_F\,z_R^{\text{rigid}}\,d(G),$$

with the donor at a gradient of 5 mmHg. The inputs are the flow load, the rigid network and the central venous
pressure; the portal pressure of the recipient is the output and is not used. The rigid value runs from 0.81 (with the
two inflows) to 1.00 (the network alone), because the published relations were measured with the hepatic artery in
place and may already contain its part.

### The relations against measurements they were not derived from

Doubling the flow of a donor liver raises the gradient by 1.9 mmHg with the cube law and 1.7 mmHg with the
compliance; the source gives about 2 mmHg for a doubling of portal flow [L6].

Ricciardi et al. [L5] measured portal pressure and resistance at two flows in pig livers in situ, isolated on a
bypass circuit, and after 2 h of cold ischaemia. With the published starting gradient and caval pressure:

| Pig liver | Flow raised | Resistance, measured | Cube law | Compliance |
|---|---:|---:|---:|---:|
| in situ (9 livers) | 2.75-fold | 0.54 | 0.48 | 0.53 |
| isolated, perfused (7 livers) | 2.67-fold | 0.59 | 0.48 | 0.57 |
| after 2 h of cold ischaemia (6 livers) | 2.67-fold | 0.59 | 0.48 | 0.50 |

A normal pig liver lies at the distensible bed or just above it (0.54 and 0.59 against 0.48 to 0.57). After cold
ischaemia the resistance still falls by the same factor when flow is raised, but it is twice as high at every flow:
relative to the isolated liver at low flow, the cold liver at 2.7 times the flow is at 1.18.

### The six groups

| Group | $z_F$ | CVP (mmHg) | Rigid network | Distensible bed, cube law | Distensible bed, compliance | Observed (95% sampling interval) | Position |
|---|---:|---:|---:|---:|---:|---:|---|
| Chan 2011 | 3.93 | 8.0 | 0.81 to 1.00 | 0.47 to 0.52 | 0.42 to 0.45 | 0.31 (0.25 to 0.37) | below both |
| Wang 2014, splenectomy | 4.12 | 3 to 9, assumed | 0.81 to 1.00 | 0.39 to 0.52 | 0.40 to 0.44 | 0.54 (0.48 to 0.60) | at the bed |
| Yagi 2006 | 3.28 | 5.2 | 0.81 to 1.00 | 0.48 to 0.53 | 0.46 to 0.50 | 0.55 (0.44 to 0.69) | at the bed |
| Wu 2011 | 3.25 | 4.4 | 0.81 to 1.00 | 0.47 to 0.52 | 0.47 to 0.51 | 0.65 (0.57 to 0.75) | between the two |
| Wang 2014, no splenectomy | 3.74 | 3 to 9, assumed | 0.81 to 1.00 | 0.41 to 0.54 | 0.43 to 0.47 | 0.72 (0.64 to 0.80) | between the two |
| Ishizaki 2012 | 3.34 | 9.1 | 0.81 to 1.00 | 0.52 to 0.58 | 0.46 to 0.50 | 0.74 (0.60 to 0.92) | between the two |

The two relations give nearly the same bed: 0.46 to 0.53 and 0.44 to 0.48 on average, and 0.39 to 0.58 over the groups.
Three groups lie between the bed and the rigid network with their whole sampling interval. Two lie at the bed: their
values, 0.54 and 0.55, are just above its upper end and their intervals reach it. Chan 2011 lies below both: its gradient
was 1 mmHg above the donor's at four times the flow, which no passive relation gives.

What this supports and what it does not. The observed resistance loads run from what a fully distensible bed gives to
short of the rigid network, so the missing resistance is of the size passive distension of the venous trees can supply.
It does not show that distension is the mechanism. The relations come from cats, in the normal range of pressure, where
their authors note that the changes of resistance are notable below 10 mmHg and minor above [L7]; a graft is at 14 to
22 mmHg. And a liver after cold ischaemia has twice the resistance of the same liver before it [L5], which moves a graft
up from the bed. The share of the bed that would have to follow the relation to give each observed value is
0.12 to 0.91 in the five groups at or above the bed; it is a description, not an estimate.

[L1] Lautt WW. Regulatory processes interacting to maintain hepatic blood flow constancy: vascular compliance, hepatic arterial buffer response, hepatorenal reflex, liver regeneration, escape from vasoconstriction. *Hepatol Res* 2007;37:891–903.
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
