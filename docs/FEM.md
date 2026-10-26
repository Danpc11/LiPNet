# The lobule as a finite-element problem

LiPNet treats the liver as a graph: a portal tree, one edge per lobule, a hepatic venous
tree. That abstraction answers the question the clinical analysis asks, which is what turns
portal flow into portal pressure, and it is deliberately blind to everything inside the
lobule. This layer opens that edge up.

It does three things the graph cannot do:

1. it computes the lobule's resistance from its geometry instead of assigning it;
2. it gives a resistance that depends on flow, through shear-mediated vasodilation, where
   the graph edge is a rigid Poiseuille tube;
3. it carries oxygen, ethanol and drugs across the lobule, so that a change in perfusion
   can be read as a change in delivery and in zonal exposure rather than only as a change
   in pressure.

Everything here is in SI units and lives in `src/fem/`.

---

## The model

### Geometry

One hexagonal lobule of Kiernan: circumradius $R = 0.55$ mm from the centre to a portal
triad, prism height $H = 2$ mm, a central vein of radius 30 µm. The domain is the hexagon
minus the central vein disc.

The mesh is structured, built directly in polar coordinates mapped onto the hexagon,

$$r(i,j) = r_{cv} + s_i\,\bigl(R_{hex}(\theta_j) - r_{cv}\bigr), \qquad s_i = (i/n_r)^{\gamma},$$

with $R_{hex}(\theta) = a / \cos\bigl((\theta \bmod \pi/3) - \pi/6\bigr)$ and apothem
$a = R\cos(\pi/6)$. No Delaunay step, no degenerate elements, exact periodicity in
$\theta$, and with $n_t$ a multiple of six the portal triads are mesh nodes rather than
interpolated points. The grading exponent $\gamma > 1$ refines towards the central vein,
where the converging flow is steepest.

Two inlet conditions are available. `perimeter` feeds the whole hexagon boundary, which is
the usual porous-lobule condition and the microanatomically honest one: interlobular
venules run along every edge. `triads` feeds six windows at the corners, which is the
point-source idealisation the LiPNet graph uses, and reproduces its stagnant zones at the
face midpoints. Both give nearly the same lobule resistance; they differ in the periportal
corners.

### Perfusion

Anisotropic Darcy flow,

$$\nabla\cdot\mathbf{q} = 0, \qquad \mathbf{q} = -\frac{\mathbf{K}}{\mu}\nabla p, \qquad
\mathbf{K} = k_r\,\mathbf{e}_r\mathbf{e}_r^{\!\top} + k_t\,\mathbf{e}_\theta\mathbf{e}_\theta^{\!\top},$$

with $p = p_{portal}$ at the inlet, $p = p_{central}$ on the central vein and no flux on the
rest of the boundary, which is the symmetry plane shared with the neighbouring lobule. The
tensor is rotated at each element centroid, so the sinusoids are more permeable along the
portal-central axis than across it. $\mathbf{q}$ is superficial; a red cell moves at
$\mathbf{q}/\phi$.

The flow is read off the Dirichlet reactions, $\mathbf{A}p$ restricted to the constrained
nodes, which balances inflow against outflow to machine precision rather than to
quadrature error.

### Transport

$$\nabla\cdot\bigl(\mathbf{q}\,g(c)\bigr) - \nabla\cdot(D\nabla c) + R(c) = 0$$

solved by damped Newton with SUPG stabilisation, $\tau = \frac{h}{2|u|}(\coth Pe - 1/Pe)$.

$g$ is the advected content. For a free solute it is the identity. For oxygen it is not:
blood carries free O₂ and the much larger haemoglobin-bound pool,

$$g(p) = \alpha p + c_{Hb}\,\frac{p^{n}}{p_{50}^{n} + p^{n}},$$

while diffusion acts only on the free pool. The capacity $dg/dp$ is tens of times $\alpha$
near $p_{50}$ and falls at high $pO_2$, which is what flattens the periportal part of the
profile and steepens the pericentral part. Treating oxygen as a diluted species misses
this, and gets the shape of the gradient wrong even after the endpoints are fitted.

### Metabolism

| species | kinetics | zonation |
|---|---|---|
| O₂ | Michaelis-Menten, $K_m \approx 1$ mmHg, $V_{max} \approx 3$ µmol min⁻¹ g⁻¹ | uniform by default |
| ethanol | ADH ($K_m \approx 0.5$ mM) + CYP2E1 ($K_m \approx 10$ mM) | CYP2E1 pericentral |
| acetaldehyde | ALDH2 ($K_m \approx 1$ µM), source = the ethanol sink | follows ethanol |
| drug | one CYP, Michaelis-Menten, with a reactive-metabolite yield | pericentral (paracetamol geometry) |

A zonation weight redistributes enzyme capacity along $\zeta$ (0 periportal, 1 pericentral)
with unit lobule mean, so sharpness can be swept at constant $V_{max}$.

### Sinusoid and hepatocyte: `two_phase.py`

The single-continuum model puts the solute and the enzyme in the same place, so whatever
is in the blood is immediately available to be metabolised. For oxygen and ethanol that is
close enough, because the sinusoidal endothelium is fenestrated and has no basement
membrane, so the barrier is a micron of water. For anything whose uptake is transporter
limited it is wrong in a way that matters: it makes the lobule look as if it clears the
drug as fast as the enzyme can work, when it clears it as fast as the membrane can deliver
it.

Two coupled fields, with the metabolism moved to where the enzymes are:

$$\nabla\cdot(\mathbf q\,g(c_v)) - \nabla\cdot(\phi D_v \nabla c_v) + \beta\,(c_v - c_t/K_p) = 0$$
$$-\nabla\cdot((1-\phi) D_t \nabla c_t) - \beta\,(c_v - c_t/K_p) + R(c_t) = 0$$

$\beta$ is not a free parameter. It is $\beta = P\,a_v$ with $a_v = 4\phi/d_{sin}
\approx 7.5\times10^{4}$ m⁻¹ the sinusoidal surface per unit volume and $P$ the
permeability of the endothelium plus the space of Disse. Free diffusion through a micron of
water gives $P\sim10^{-3}$ m s⁻¹ and $\beta\sim10^{2}$ s⁻¹; a transporter-limited drug can
sit at $P\sim10^{-7}$ m s⁻¹.

The membrane need not be uniform. Periportal sinusoids have more and smaller fenestrae,
pericentral ones fewer and larger, and `zonal_permeability` redistributes $P$ along $\zeta$
at constant lobule mean, so a sweep in sharpness isolates the zonation from a net gain or
loss of surface. The result is mild and has a clear structure:

| extraction | membrane periportal | uniform | membrane pericentral |
|---|---|---|---|
| **enzyme periportal** | **0.158** | 0.157 | 0.146 |
| **enzyme uniform** | 0.155 | **0.159** | 0.155 |
| **enzyme pericentral** | 0.144 | 0.154 | **0.160** |

The diagonal wins: a membrane matched to the enzyme clears more than a mismatched one. But
the penalty for full mismatch is about 10%, against the five-fold effect of the
permeability level itself. Zonal fenestration is a second-order correction, which is worth
knowing before anyone builds an argument on it.

The regime is decided by an exchange Damköhler number built on the transit time,
$Da_x = \beta\,\tau_{transit}$, with $\beta$ divided by the tissue capacity when the
unknown is a partial pressure rather than a concentration. Sweeping the membrane:

| $P$ [m/s] | $10^{-3}$ | $10^{-5}$ | $10^{-6}$ | $10^{-7}$ | $10^{-8}$ |
|---|---|---|---|---|---|
| $Da_x$ | 645 | 6.4 | 0.64 | 0.064 | 0.0064 |
| extraction | 0.208 | 0.208 | 0.202 | 0.154 | 0.038 |
| hepatocyte / blood | 1.000 | 0.995 | 0.948 | 0.609 | 0.103 |

The curve is flat for four decades and then falls off a cliff at $Da_x \approx 1$. Above it
the answer is $V_{max}$ and the membrane is invisible; below it the answer is $\beta$ and
$V_{max}$ is irrelevant. Nothing between the two regimes needs fitting to see where the
boundary is.

For oxygen $Da_x \approx 1300$ and the blood-to-hepatocyte gap is 0.25 mmHg, so the
two-continuum model reproduces the single-continuum extraction to 0.6% and the gap closes
like $1/P$. That reduction is the test that the block solver is right, and it is in the
suite.

### Nitric oxide

NO is deliberately **not** a lobule-scale field. With haemoglobin in the lumen the
scavenging rate is $k[\mathrm{Hb}] \sim 5\times10^{7}\,\mathrm{M^{-1}s^{-1}} \times
10^{-2}\,\mathrm{M} \sim 5\times10^{5}\,\mathrm{s^{-1}}$, giving a diffusion length of
about 0.18 µm; in the perivascular tissue it is 5.7 µm. Both are below one element of any
lobule mesh anyone would run. An advection-diffusion-reaction field for NO on this mesh
would resolve nothing.

What is solved instead is one cylinder around a single sinusoid: red-cell core, plasma free
layer, endothelium with shear-driven eNOS release $q = q_0(\tau_w/\tau_{ref})^{\,\text{gain}}$,
and tissue out to $b = a/\sqrt{\phi}$ with zero flux by symmetry. Every consumption term is
first order, so the problem is **linear in the production**: one tridiagonal solve gives a
geometric factor $\chi = [\mathrm{NO}](r_{target})/q$, and the lobule field is then
$\chi q_0 (\tau_w/\tau_{ref})^{\text{gain}}$ element by element. One radial solve for the
whole lobule, not one per element.

The sub-model reproduces the thing that makes NO signalling possible at all: the plasma
free layer. Without it, the wall concentration collapses by eight orders of magnitude into
the red-cell core.

---

## Two numerical points worth knowing

Both of these were found by tests in `tests/test_fem.py` rather than by inspection, and
both change the answer, not just its accuracy.

**The advection must be conservative.** Written as $u\cdot\nabla g$, the operator is
equivalent to the conservative form only if the discrete velocity is pointwise divergence
free. A P1 Darcy flux is not: it is divergence free weakly, with jumps across faces. In the
non-conservative form the mass lost scales with the *magnitude* of $g$ rather than with its
variation, and for oxygen $g \approx 8$ mol m⁻³ while its variation across the lobule is
about 2. The lobule balance then got *worse* under refinement, from 9% to 19%, while the
field itself looked perfectly plausible.

**The boundary term must use nodal fluxes, not edge integrals.** Per element,
$-\int_K \nabla\phi_i\cdot u + \oint_{\partial K}\phi_i\, u\cdot n = 0$ exactly, but summing
that over elements leaves the jumps of $u$ across interior faces. An edge integral
therefore preserves a constant only to $O(h)$, with the whole error piled up at the central
vein where the flux jumps most: 11% on a coarse mesh. The nodal Darcy outflow
$q_n = -(\mathbf{A}p)$ is by definition the quantity that cancels $-\int\nabla w_i\cdot u$
node by node, so with it a uniform inlet and no reaction give the constant solution to
$10^{-14}$, and

$$\text{influx} = \text{outflux} + \text{consumption}$$

becomes an identity rather than an approximation. The geometric boundary integral is kept
as an independent $O(h)$ check and reported as `balance_error`.

---

## What it reproduces

Nothing below was fitted. The permeability was chosen to put the perfusion at
1 mL min⁻¹ g⁻¹ and the porosity to put the transit time in range; everything else follows.

| quantity | model | reported range |
|---|---|---|
| perfusion | 1.05 mL min⁻¹ g⁻¹ | ~1 |
| sinusoidal transit | 8.6 s | 5-10 s |
| sinusoidal velocity | 54 µm s⁻¹ | 20-200 |
| wall shear stress | 0.19 Pa (mean), 1.0 Pa (p90) | 0.1-0.5 |
| $pO_2$, periportal → pericentral | 56 → 33 mmHg | 60-65 → 30-35 |
| haemoglobin saturation | 0.91 → 0.63 | — |
| hepatic O₂ extraction | 33% | 20-40% |
| ethanol first pass | 7.6%, of which 14% via CYP2E1 | — |
| acetaldehyde cleared within the lobule | 99.9% | blood levels stay in the low µM |

The reactive drug metabolite and the acetaldehyde peak both land at $\zeta = 1$, which is
where centrilobular necrosis is.

Perivascular NO comes out **higher pericentrally** (3.3 nM) than periportally (0.65 nM).
That is not an input: flow converges towards the central vein, so wall shear rises along the
sinusoid, and eNOS follows it. It is a prediction of the geometry.

---

## What goes back to the network

### The shape factor

LiPNet's absolute conductance scale comes from the cost budget, not from physiology, so
only a *relative* correction is transferable. The reference is an equal-area radial annulus
with the same axial permeability,

$$R_{annulus} = \frac{\mu}{k_r}\,\frac{\ln(b/r_{cv})}{2\pi H}, \qquad \pi b^2 = A_{hex},$$

and $\kappa = R_{FEM}/R_{annulus}$ is a pure number set by the hexagonal outline, the
anisotropy and the size of the central vein. For the default lobule $\kappa = 0.996$: the
hexagon is almost exactly its equal-area annulus. That is a real result, and a slightly
boring one. It says the single-edge idealisation is not costing the network model anything
geometric, which is worth knowing before anyone proposes that it is.

`coupling.apply_to_graph` rescales the sinusoidal edges of a graph from
`model.network.build_liver_graph`, either through the length (so the cost weight
$a_e = L^{3b/2}$ follows) or as a pure conductance multiplier.

### The non-linear part

This is where the layer earns its place. Drive the lobule at a prescribed flow instead of a
prescribed pressure, close the shear → NO → calibre → permeability loop at each flow, and
the lobule's own resistance load is

$$z_R^{lobule}(z_F) = \frac{R(z_F Q_0)}{R(Q_0)}.$$

| $z_F$ | 0.5 | 1 | 1.5 | 2 | 3 | 4 |
|---|---|---|---|---|---|---|
| $z_R$ | 1.21 | 1.00 | 0.85 | 0.74 | 0.59 | 0.56 |

With the gain set to zero this is identically 1, which is the rigid Poiseuille edge the
graph currently uses, so the two models are nested and the test suite checks it. Mean wall
shear scales exactly with $z_F$, as it must at prescribed flow whatever the calibre does.

The interpretation needs care. The published $z_R < 1$ values come from grafts with
reconstructed or enlarged venous outflow, and that is a macroscopic effect this model says
nothing about. What the lobule adds is an independent, same-signed contribution: a
hyperperfused graft dilates its own sinusoidal bed, so part of the measured drop in
resistance load is microcirculatory and would appear even with the outflow untouched. The
model gives the *shape* of $R(Q)$ rather than assuming linearity, and it bounds how much of
the effect the microcirculation can carry. The clamp on the permeability multiplier, 0.6 to
1.8, is what keeps the curve from running away; it is the least constrained parameter here
and the one a sensitivity sweep should hit first.

---

## What the network does to all this

`network_link.py` asks the question the clinical index forces: the lobule's $z_R$ is
measured across the whole graft, and the sinusoidal bed is one resistance in a series that
also holds the portal tree, the hepatic venous tree and two anastomoses. How much of the
effect survives?

First, how much of the resistance is sinusoidal at all. Taking the dissipation-optimal
network of `src/model/network.py` and splitting its dissipation by edge type:

| lobules | 37 | 73 | 151 |
|---|---|---|---|
| sinusoidal share | 0.247 | 0.220 | 0.205 |

It falls like $N^{-0.13}$. Extrapolated to the $10^{6}$ lobules of a real liver that is
about 0.07, which is below everything that has been measured. That looked like a problem
with the model. It is not: it is a constraint on the cost exponent.

Running the same campaign for each $b$ (`share_vs_cost_exponent`):

| $b$ | 37 | 73 | 151 | 301 | slope | extrapolated to $10^6$ |
|---|---|---|---|---|---|---|
| 0.50 | 0.308 | 0.288 | 0.284 | 0.274 | **-0.053** | **0.18** |
| 0.75 | 0.247 | 0.220 | 0.206 | 0.188 | -0.126 | 0.068 |
| 1.00 | 0.196 | 0.165 | 0.141 | 0.120 | -0.232 | 0.018 |

The drift with size is not a fixed property of the model, it is set by $b$. At $b = 1$ the
share collapses with liver size; at $b = 0.5$ it is close to scale invariant, and that is
the only branch that survives extrapolation inside the measured range. Since the measured
shares in cat, rabbit and rat livers, which differ by orders of magnitude in lobule count,
are all of the same order, approximate scale invariance is itself an observation, and it
points at the low end of $b$.

That turns a loose end into a fifth prediction: **the sinusoidal fraction of hepatic
vascular resistance is an independent observable that constrains the vascular maintenance
exponent**, with no clinical data involved. It is worth checking against the cost-convexity
work directly, since $b$ is exactly the parameter at issue there.

The campaign above covers 37 to 301 lobules. It should be pushed further before the
constraint is quoted as a number: the lever arm is one decade and the extrapolation is
four. Note also that this is a property of the cost-constrained
optimum: the optimiser assigns every conductance under a budget, so this is what the model
*thinks* the share should be, not what a vascular cast would show.

### What the measurements say

The share is measurable, and it has been measured, with a long disagreement about the
answer.

| source | method | sinusoidal share |
|---|---|---|
| Lautt et al. 1986 (cat); Lautt et al. 1987 (dog) | localisation of the resistance site | ~0 |
| Nakata et al. 1997 (rabbit) | servo-null micropuncture | 0.08 - 0.26 |
| Shibayama & Nakata, isolated perfused rat liver | portal / hepatic venous occlusion | ~0.33 |
| **this network, 37 to 151 lobules** | dissipation at the optimum | **0.21** |

Lautt's group put essentially the whole portal-to-caval drop in a short postsinusoidal
segment of the hepatic veins, with portal and lobar venous pressures indistinguishable, so
for them the sinusoidal share is nil. The rabbit micropuncture, which measures pressures
vessel by vessel rather than inferring them, finds a gradient of under 1 mmHg across the
sinusoids out of a total of 3.8, and concludes that most of the transhepatic resistance
cannot be attributed to any one vascular segment. The perfused rat gives roughly half
presinusoidal, a third sinusoidal and a fifth postsinusoidal.

Two things follow for this model.

The first is a validation that was not aimed for: the cost-constrained optimum lands at
0.21, inside the measured range and close to the micropuncture, without the share having
been used anywhere in the optimisation. The network model predicts an observable it was
never fitted to.

The second is that **the half-of-all-resistance scenario is not supported in a normal
liver.** It is a cirrhotic pattern, and the grafts in `data/series.tsv` come from healthy
living donors; it is the recipient who is cirrhotic, not the graft. So the last column of
the table above should be read as a sensitivity bound and not as a live possibility, and
the conclusion of this section hardens rather than softens.

### Where the resistance sits inside the lobule

Henriksen and Lassen (Liver 1988;8:88-94) modelled the pressure profile along a converging
sinusoidal network and concluded that the main hydrodynamic resistance is located
downstream, at the outlet, and that the precise shape of the microcirculatory unit does not
much matter. This FEM reproduces that without being asked to. Of the total intralobular
pressure drop:

| $\zeta$ | 0 to 0.5 | 0.5 to 0.75 | 0.75 to 1 |
|---|---|---|---|
| share of the drop | 27% | 26% | **47%** |

Half the drop happens in the last quarter, next to the central vein, for a purely
geometric reason: the flow converges, so the cross-section available falls as $1/r$.

This also goes some way towards reconciling the two schools above. If the resistance
concentrates at the lobular outlet and in the hepatic venules immediately beyond it, then a
catheter placed one way calls it sinusoidal and placed the other way calls it
postsinusoidal, and both groups may be describing the same region. That is an argument the
model can make on its own evidence.

Then the series calculation of `whole_graft.py`, repeated with the sinusoidal part made
flow dependent and solved as a fixed point, since a lower resistance admits more flow and
more flow dilates further:

| $z_F$ | lobule $z_R$ | graft $z_R$, share 0.21 | graft $z_R$, share 0.07 | graft $z_R$, share 0.5 |
|---|---|---|---|---|
| 1.5 | 0.85 | 0.97 | 0.99 | 0.93 |
| 2 | 0.74 | 0.95 | 0.98 | 0.88 |
| 3 | 0.59 | 0.91 | 0.97 | 0.80 |
| 4 | 0.56 | 0.91 | 0.97 | 0.79 |

The transmitted fraction is exactly $\text{share}\times r_{lob}/(r_{lob}+r_{out})$, which
the test suite checks against the solver.

**So the answer is mostly negative, and that is the useful part.** The lobule effect is
real, it has the right sign, and it is too small to matter: a graft at four times the
donor's flow per gram reaches a measured $z_R$ of 0.91, not the 0.31 to 0.74 of the
published series. Even if the sinusoids carried half of all hepatic resistance, twice the
model's own estimate, it would only reach 0.79. Sinusoidal autoregulation cannot account
for the clinical observation, which leaves the outflow reconstruction as the explanation
the data already pointed to. What the lobule adds is a small, same-signed background term
that should be subtracted before attributing the rest, not a competing mechanism.

The escape route is now closed. The sinusoidal share was the weakest number in the chain,
and the measurements bracket it between 0 and a third in a normal liver, with this model's
own 0.21 sitting in the middle. There is no version of the healthy graft in which the
sinusoidal bed carries enough of the resistance to deliver a $z_R$ of 0.31.

### The sign is wrong at high flow, and the data say so

The feedback implemented here is adaptive dilation: more shear, more eNOS, wider bed. The
measurements in the regime that matters say the opposite. In the porcine small-for-size
model of Kelly et al. (*Liver Transpl* 2010), with grafts at 19 to 25% of the recipient's
standard liver volume, portal flow, pressure gradient **and resistance** were all
significantly higher than in whole-liver recipients after reperfusion, with endothelial
injury markers rising first and hepatocellular ones later. Resistance up, not down. There
are clinical case reports of the same thing, a transient rise in intrahepatic portal
resistance in a small-for-size graft, severe enough to reverse portal flow.

So the sign flips somewhere. `Feedback.injury_shear` puts that fork in the model: above a
wall shear stress threshold the calibre is lost instead of gained, measured against the
baseline like the NO term so the reference state stays its own reference. With a threshold
at 0.4 Pa, about twice the baseline mean:

| $z_F$ | 1 | 1.5 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| lobule $z_R$, dilation only | 1.00 | 0.85 | 0.74 | 0.59 | 0.56 |
| lobule $z_R$, with injury | 1.00 | 1.30 | 1.32 | 1.36 | 1.40 |
| graft $z_R$, dilation only | 1.00 | — | 0.95 | 0.92 | 0.91 |
| graft $z_R$, with injury | 1.00 | — | 1.06 | 1.07 | 1.08 |

Neither the threshold nor the slope is known, so the branch is off by default and exists to
let the fork be explored, not to assert a number. What the fork does establish is that the
two branches straddle 1 and both stay within 9% of it at the graft level, which is the next
section's point.

### Against the five published groups

`predict_series` runs the prediction over the groups in `data/series.tsv` that report both
indices:

| group | outflow | $z_F$ | $z_R$ observed | $z_R$ predicted |
|---|---|---|---|---|
| Chan 2011 | reconstructed | 3.93 | 0.306 | 0.914 |
| Yagi 2006 | reconstructed | 3.28 | 0.549 | 0.918 |
| Wang 2014 | reconstructed | 4.12 | 0.539 | 0.913 |
| Wang 2014 | standard | 3.74 | 0.716 | 0.915 |
| Ishizaki 2012 | reconstructed | 3.34 | 0.742 | 0.917 |

The flow load barely varies across them, from 3.3 to 4.1, so the lobule model predicts the
same $z_R$ to three decimals for all five. The observed spread is 0.436; the predicted
spread is 0.004. A factor of a hundred.

This is the cleanest form of the result. Whatever separates a $z_R$ of 0.31 from one of
0.74 in these series, it is not the flow through the lobules, because the flow through the
lobules is nearly the same in all of them. It also means the data cannot discriminate
between the dilation and the injury branch above: both sit within 9% of 1 while the data
range over 44%.

---

## Running it

```bash
python src/fem/run_fem.py perfusion      # Darcy, resistance, transit, shear
python src/fem/run_fem.py transport      # O2, ethanol, acetaldehyde, drug
python src/fem/run_fem.py twophase       # sinusoid/hepatocyte exchange: membrane vs enzyme
python src/fem/run_fem.py no             # the NO wall sub-model
python src/fem/run_fem.py coupling       # shape factor and z_R(z_F)
python src/fem/run_fem.py network        # how much of it survives the rest of the series
python src/fem/run_fem.py figures
python src/fem/run_fem.py all
python tests/test_fem.py
```

Options: `--nr --nt` (mesh), `--inlet perimeter|triads`, `--gain` (NO feedback), `--quiet`.
Results land in `results/fem_*.json`, figures in `results/fem_*.pdf` and `.png`.

No new dependencies: NumPy, SciPy and Matplotlib only, all already required. The default
mesh is 3528 nodes and the whole pipeline runs in under a minute on one core; the coupling
sweep is the slow part, at a few seconds per flow point.

---

## References for the numbers used here

Values taken from the literature, by section. Nothing in this layer was fitted to the
clinical series in `data/series.tsv`.

Lobule geometry, sinusoidal permeability and porosity: Debbaut et al., *J Biomech Eng*
2012 and *IEEE Trans Biomed Eng* 2014. Porous-medium formulations of the lobule: Rani
et al. 2006; Bonfiglio et al., *Am J Physiol* 2010.

Oxygen gradient and metabolic zonation: Jungermann & Kietzmann, *Hepatology* 2000;
Kietzmann, *Redox Biol* 2017. Ethanol kinetics, ADH and CYP2E1: Lieber, *Alcohol Res
Health* 2004.

Nitric oxide near a vessel wall, red-cell scavenging and the plasma free layer: Vaughn,
Kuo & Liao, *Am J Physiol* 1998; Buerk, *Annu Rev Biomed Eng* 2001.

Portal hyperperfusion raising graft resistance: Kelly et al., *Liver Transpl* 2010
(porcine small-for-size model); Kiuchi/Ito case reports of transient intrahepatic
resistance rise with hepatofugal flow.

Localisation of hepatic vascular resistance: Lautt, Greenway, Legare & Weisman,
*Am J Physiol* 1986;251:G375 (cat); Lautt et al., *Can J Physiol Pharmacol* 1987 (dog);
Nakata et al., *Am J Physiol* 1997;272:H2826 (rabbit servo-null micropuncture);
Shibayama & Nakata, isolated perfused rat liver, segmental resistances by occlusion.
Intralobular pressure profile and outlet-dominated resistance: Henriksen & Lassen,
*Liver* 1988;8:88-94.

Portal hyperperfusion, sinusoidal shear and outflow reconstruction in small-for-size
grafts: see the series already collected in `data/series.tsv` and the reviews cited in
`THEORY.md`.

---

## Limits of this layer

The lobule is 2D. The prism is handled by multiplying by the height $H$, which assumes
nothing varies along the central vein axis. That is fine for resistance and for zonal
averages and wrong near the ends of the lobule.

The sinusoidal bed is a porous continuum. The permeability and porosity are effective
parameters calibrated against the perfusion and the transit time, not measured. Anything
that depends on an individual sinusoid, including the shear used to drive eNOS, goes
through the Poiseuille estimate $\tau_w = 8\mu v/d$ and inherits whatever error is in the
assumed diameter.

The two-phase model shares the lobule geometry with the single-phase one, so the exchange
area enters only through $a_v = 4\phi/d_{sin}$. Fenestration density, the width of the
space of Disse and anything else that varies zonally are folded into one number, $P$, which
is taken as uniform. A zonal $P$ would be straightforward to add and is probably the first
refinement worth making, since fenestration is known to differ between zone 1 and zone 3.

The NO sub-model assumes the lumen is well mixed and the production is uniform along the
sinusoid. It gives the right scales and the right dependence on shear, not a quantitative
concentration; and since the feedback uses only $NO/NO_{ref}$, $q_0$ cancels out of
everything reported here.

Nothing is calibrated against liver data from a graft. The $z_R(z_F)$ curve is a prediction
of the mechanism, not a fit to the clinical series, and should be read that way.
