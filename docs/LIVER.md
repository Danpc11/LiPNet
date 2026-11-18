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
budget of one tree relative to the other, which the model does not supply. Results are reported at each exponent,
and where they depend on the split between the trees, over every split.

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
| each tree at its cast exponent, any split | 0.95 to 1.00 |
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

The series report no arterial flow, so $\varphi$ is not observed.

## 6. The model against the five paired groups

| | $z_R$ |
|---|---|
| network, lobar graft | 0.93 to 1.00 |
| explicit graphs | 0.90 to 0.96 |
| with the arterial inflow, $\varphi\le1$ | 0.81 to 0.90 |
| with the shear feedback of the lobule layer (`src/fem/network_link.py`) | 0.91 to 0.92 |
| observed | 0.31, 0.54, 0.55, 0.72, 0.74 |

All five observed values lie below the lowest value the model reaches. The model has the sign of the departure and
not its size. Widening the outlet anastomosis in `whole_graft.py` lowers $z_R$ to 0.95 at most, because that
anastomosis holds 5% of the resistance.

## 7. What the table supports

- Arithmetic of `data/series.tsv`: no disagreement between percentages and counts, or between stored gradients and
  the pressures they come from.
- The donor flow reference is measured in three series (81, 91 and 118 mL/min/100 g) and set to 90 elsewhere. No
  series reports a donor gradient; the normal gradient of 5 mmHg is assumed for all 21 groups with a pressure load.
  The table does not record when the flow was measured relative to the pressure.
- Sampling intervals of $z_R$ with the dispersions the papers report (largest reported value where none is given)
  stay below one: 0.25 to 0.37, 0.44 to 0.69, 0.48 to 0.60, 0.64 to 0.80, 0.60 to 0.92.
- The reference values decide the level. A group reaches $z_R=1$ if the donor gradient is at or below 1.5, 2.7, 2.7,
  3.6 and 3.7 mmHg respectively; at the largest measured donor reference (118) the two highest groups are at 0.94
  and 0.97.
- Without any reference: the flow load spans 1.26-fold across the five groups and the pressure load 2.23-fold;
  their rank correlation is −0.10; 10% of the variance of $\log z_P$ is in $\log z_F$. In Wang 2014 the arm with
  10% more flow has 12% less portal pressure.

## 8. What the outcome model supports

| stratum | groups defined by | span of the predictor | slope on its own (SE) | events |
|---|---|---|---|---|
| Osman 2017 | final portal pressure | 5.0 mmHg | +2.00 (1.10) | 7/76 |
| Wang 2014 | splenectomy | 2.3 mmHg | +1.78 (0.76) | 39/276 |
| Uemura 2016 | graft size | 0.37 mmHg | +4.51 (7.35) | 30/221 |
| Ogura 2010 | final portal pressure | 3.7 mmHg | +2.52 (0.72) | 20/129 |
| Yagi 2005 | portal pressure, days 1 to 3 | 8.3 mmHg | +1.37 (0.69) | 8/32 |

Uemura 2016 separates its groups by 0.37 mmHg against a within-group standard deviation of 3 to 4 mmHg and carries
no information on the slope: the primary estimate is +1.55 (+0.33 to +2.79) with it and +1.53 (+0.30 to +2.77)
without it, on 46 events in 352 recipients from two series. The same model fitted to the flow load gives +0.86
(−0.58 to +2.36) on 48 events; the probability that the pressure slope exceeds the flow slope is 0.77. The network
model enters this analysis only through the normalisation of the predictor: it contains no relation between load
and outcome.
