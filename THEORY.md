# LiPNet: theory

The model implemented in `src/`. Symbols are listed where they first appear.

## 1. The perfusion network

The liver is a directed graph: a source at the portal vein feeds an inflow tree, each of $N$ lobules is crossed by one
sinusoidal path, and an outflow tree drains through the hepatic veins to the inferior vena cava. Edge $e$ carries
flow $f_e$ and has radius $r_e$, length $l_e$, conductance $w_e=r_e^4/l_e$ and vascular mass $m_e=r_e^2 l_e$.

The network is the one that minimises dissipation at a fixed maintenance budget,

$$\min\ D=\sum_e \frac{f_e^2}{w_e}\qquad\text{subject to}\qquad C=\sum_e m_e^b,$$

with $b$ the maintenance exponent ($b=1$: cost proportional to vessel volume; $b<1$: large vessels cheaper). Results
are computed over $b\in[0.5,1]$. For a fixed topology the optimal conductances are analytic; the tree is found by a
greedy construction refined by local swaps (`model/network.py`, `model/network3d.py`) and checked by enumeration on
seven lobules (`model/certify.py`).

At the optimum $D_e/c_e$ is equal on every edge, and its common value defines the scale of the design

$$\tau_0=\sqrt{D^*/C}.$$

The Poiseuille wall shear stress of an edge is

$$\tau_e=\tau_0\,m_e^{(b-1)/2},$$

so $\tau_0$ is a wall shear stress only at $b=1$. In the limit of many lobules, for $b<1$,

$$D^*\propto F^2\,C^{-2/b}\,N^{(2-b)/b}\,L^3,$$

with $F$ the total flow and $L$ the lobule spacing (`model/scaling.py`). On the symmetric canopy with branching number
$n$ and length ratio $\lambda$, the radius ratio between successive generations is

$$\frac{r_{i+1}}{r_i}=n^{-1/(2+b)}\,\lambda^{(1-b)/(2(2+b))}.$$

The pressure drop across a lobule is local:

$$\Delta P = f\,R.$$

## 2. Normalised loads

A partial graft keeps a fraction $g$ of the network and receives $h$ times the donor inflow, so each surviving lobule
carries $h/g$ times its donor flow. In clinical quantities, each equal to 1 in a healthy donor:

$$z_F=\frac{\text{graft PVF per }100\,\text{g}}{\text{donor PVF per }100\,\text{g}},\qquad
z_P=\frac{\text{PVP}-\text{CVP}}{5\ \text{mmHg}},\qquad
z_R=\frac{R_{\text{graft}}}{R_{\text{donor}}}.$$

Since $R=\Delta P/Q$ per 100 g,

$$z_P = z_F\,z_R.$$

The identity holds by definition. Flow and pressure carry the same information only while $z_R=1$. The donor
references are 90 mL/min/100 g of portal flow (or the value a series measured in its donors) and a portocaval
gradient of 5 mmHg; where a series reports no CVP it is taken as 5 mmHg (`indices.py`).

## 3. A partial graft in the network

A lobar graft is the sub-network below one first-order branch. Its vessels keep their conductance, so

$$z_R(\text{graft})=1-\frac{\Delta P_{\text{removed generations}}}{\Delta P_{\text{donor}}}$$

(`model/liver.py`, `lobar_graft`). The hepatic artery joins the portal blood at the sinusoid: the portal tree carries
portal flow, the sinusoids and hepatic veins carry portal plus arterial flow. With $s_{in}$ the share of the portal
tree in the donor drop, $a$ the arterial-to-portal flow ratio of the donor (1/3 for an arterial share of 25%) and
$\varphi$ the arterial flow per lobule of the graft relative to the donor,

$$z_P=z_R^{\text{net}}\left[s_{in}z_F+(1-s_{in})\frac{z_F+\varphi a}{1+a}\right]$$

(`two_inflows`). $\varphi=0$: no arterial flow; $\varphi=1$: the donor's; $\varphi=z_F$: arterial flow scaling with
portal flow.

## 4. Physical units

`model/human.py` evaluates the canopy for an adult liver: 1500 g, density 1.05 g/mL, portal flow 90 mL/min/100 g,
arterial share 0.25, gradient 5 mmHg, blood viscosity 3.5 mPa s, $3^{13}$ lobules in 14 generations with branching
number 3 in three dimensions, branch length 0.44 times the side of the territory supplied, and a lobular drop of
1 mmHg at donor flow shared with the lobule model (`src/fem/params.py`). Conductance is Poiseuille,
$\pi r^4/(8\mu l)$.

## 5. Whole grafts

A whole graft has $g=1$, so $z_F=h$. The two anastomoses are resistances in series with the lobular bed, and the
splanchnic bed drains in parallel through collaterals $R_{coll}$:

$$R = r_{in}d_{in}^{-4} + r_{lob} + r_{out}d_{out}^{-4},$$

with $r_{in}$, $r_{out}$ the shares of the normal graft resistance held by each anastomosis and $d$ the diameter
relative to normal (`model/whole_graft.py`). Narrowing the inlet lowers flow and sinusoidal pressure together
($z_R\le1$) and raises the pressure of the recipient portal trunk; narrowing the outlet raises the pressure every
lobule drains against at nearly constant flow ($z_R>1$).

## 6. Passive distensibility

The network is rigid. `model/distension.py` adds three published relations of the normal liver, each a factor $d$ on
the resistance of the rigid network, with $G$ the portocaval gradient and $g_0=5$ mmHg:

| Relation | Form | Parameter |
|---|---|---|
| cube law | $R\propto P_d^{-3}$, $P_d=\text{CVP}+G/2$ | none |
| compliance | blood volume linear in intrahepatic pressure, $R\propto V^{-2}$: $d=\left(1+\dfrac{G-g_0}{p}\right)^{-2}$ | doubling pressure $p$, 8 to 13 mmHg |
| pressure-flow | $G/g_0=f+(1-f)\,q$ | flow-independent share $f$, 0.30 to 0.35 |

For the first two the gradient solves

$$G = g_0\, z_F\, z_R^{\text{rigid}}\, d(G),$$

whose right side minus $G$ decreases strictly with $G$, so the root is unique. The third is applied to the effective
load $q=z_F z_R^{\text{rigid}}$. Every relation acts on the load of the rigid network and returns $z_P=z_F=z_R=1$ for
the donor.

**Vessel by vessel** (`model/compartments.py`). Generation $i$ of the network keeps its own index,
$R_i=IC_i/P_{d,i}^3$, with $P_{d,i}$ the mean pressure in that generation from the pressure profile of the network,
or changes its volume with its own change of pressure (compliance). Only the vessels of the chosen sites distend:

$$\frac{R}{R_0}=(1-s)+s\,\bar d,$$

with $s$ the share of the resistance in the distending vessels (all vessels; pre- and postsinusoidal sites; hepatic
venous tree; portal tree) taken from the network. The coupled system pressure → calibre → resistance → pressure is
solved as a fixed point for the factor of each generation, with a convergence check on the original equations.

## 7. The lobule

`src/fem/` solves one hexagonal lobule (circumradius 0.55 mm, height 2 mm, central vein 30 µm) by finite elements.

Perfusion, anisotropic Darcy flow:

$$\nabla\cdot\mathbf q=0,\qquad \mathbf q=-\frac{\mathbf K}{\mu}\nabla p,\qquad
\mathbf K=k_r\,\mathbf e_r\mathbf e_r^{\top}+k_t\,\mathbf e_\theta\mathbf e_\theta^{\top}.$$

Transport, with $g(c)$ the advected content (free plus haemoglobin-bound for oxygen,
$g(p)=\alpha p + c_{Hb}\,p^n/(p_{50}^n+p^n)$) and $R(c)$ zonated Michaelis-Menten kinetics:

$$\nabla\cdot\big(\mathbf q\,g(c)\big)-\nabla\cdot(D\nabla c)+R(c)=0.$$

Sinusoid and hepatocyte as two phases exchanging through a membrane of permeability $\beta$:

$$\nabla\cdot(\mathbf q\,g(c_v))-\nabla\cdot(\phi D_v\nabla c_v)+\beta\,(c_v-c_t/K_p)=0,\qquad
-\nabla\cdot((1-\phi)D_t\nabla c_t)-\beta\,(c_v-c_t/K_p)+R(c_t)=0.$$

Nitric oxide is solved on one cylinder around a sinusoid with shear-driven endothelial release
$q=q_0(\tau_w/\tau_{ref})^{\text{gain}}$; the resulting calibre change closes the loop shear → NO → permeability and
gives the resistance of the lobule as a function of flow, $z_R^{\text{lob}}(z_F)=R(z_F Q_0)/R(Q_0)$, which
`fem/network_link.py` returns to the network.

## 8. Outcome model

For group $g$ in stratum $s$ (a study and an outcome definition), with $n_g$ recipients and $y_g$ events,

$$y_g\sim\text{Binomial}\big(n_g,\ \text{logit}^{-1}(\alpha_s+\beta_s\,(z_{P,g}-\bar z_{P,s}))\big),$$

$$\beta_s\sim\mathcal N(\mu_\beta,\tau_\beta^2),\qquad \alpha_s\sim\mathcal N(0,2.5^2),\qquad \mu_\beta\sim\mathcal N(0,1.5^2),\qquad \tau_\beta\sim\text{HalfNormal}(0.5),$$

fitted by NUTS (`bayes.py`). Centring within stratum separates the within-study slope from the level of each study.

With the within-group spread of the load, each patient carries $x\sim\mathcal N(\mu_g,\sigma_g^2)$ and the group
probability is the average of the individual curve,

$$\bar p_g=\int \text{logit}^{-1}\big(\alpha_s+\beta_s(x-\bar z_{P,s})\big)\,\mathcal N(x;\mu_g,\sigma_g^2)\,dx,$$

evaluated by Gauss-Hermite quadrature (`ipd/ecological.py`); $\sigma_g=0$ returns the model above.

A centre recalibrates the level with its own rate and mean load, keeping the slope:
$\alpha=\text{logit}(\text{rate})-\beta\,\bar z_P$ (`indices.baseline_from_rate`, `indices.recalibrate`).
