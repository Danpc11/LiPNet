"""Physical parameters of the hepatic lobule, in SI units, with the source of each value.

Everything the FEM layer needs lives here so that a calibration or a sensitivity sweep
touches one file. Pressures are stored in Pa; the helpers mmHg() and to_mmHg() convert.

Scales that the defaults are meant to reproduce (these are the acceptance targets used in
tests/test_fem.py):

    perfusion            ~1 mL/min per g of liver        -> Q_lobule ~ 2.7e-11 m^3/s
    sinusoidal transit   ~5-10 s
    sinusoidal shear     ~0.1-0.5 Pa
    pO2 periportal       60-65 mmHg  ->  pericentral 30-35 mmHg
    hepatic O2 extraction 20-40 %

References for the ranges (not reproduced here, see docs/FEM.md):
  Debbaut et al. 2012, 2014 (sinusoidal permeability, lobule geometry)
  Rani et al. 2006; Bonfiglio et al. 2010 (porous-medium models of the lobule)
  Kietzmann 2017 (metabolic zonation); Jungermann & Kietzmann 2000 (O2 gradient)
  Lieber 2004 (ADH / CYP2E1 kinetics); Vaughn, Kuo & Liao 1998 (NO wall model)
"""
from dataclasses import dataclass

MMHG = 133.322387415          # Pa per mmHg


def mmHg(x):
    """mmHg -> Pa."""
    return x * MMHG


def to_mmHg(x):
    """Pa -> mmHg."""
    return x / MMHG


@dataclass
class Geometry:
    """Hexagonal lobule of Kiernan: six portal triads at the corners, one central vein."""
    R: float = 0.55e-3          # circumradius, corner (portal triad) to centre [m]
    r_cv: float = 3.0e-5        # radius of the central vein lumen [m]
    H: float = 2.0e-3           # lobule height (prism axis) [m]
    inlet: str = 'perimeter'    # 'perimeter' (interlobular venules) or 'triads' (six corners)
    triad_arc: float = 0.35     # half-width of the portal inlet window, used when inlet='triads' [rad]
    nr: int = 48                # radial elements
    nt: int = 72                # angular elements (multiple of 6: corners land on nodes)
    grade: float = 1.35         # radial grading exponent (>1 refines towards the central vein)


@dataclass
class Flow:
    """Darcy (superficial) flow through the sinusoidal bed."""
    mu: float = 3.5e-3          # apparent blood viscosity in the sinusoids [Pa s]
    k_r: float = 4.0e-14        # permeability along the portal-central axis [m^2]
    k_t: float = 1.3e-14        # permeability transverse to it [m^2]
    phi: float = 0.15           # sinusoidal volume fraction (porosity) [-]
    d_sin: float = 8.0e-6       # sinusoid diameter, used only for the wall shear estimate [m]
    p_portal: float = mmHg(7.0)   # portal triad pressure [Pa]
    p_central: float = mmHg(3.0)  # central vein pressure [Pa]


@dataclass
class Oxygen:
    """O2 as pO2 [mmHg]. Advection carries total O2 (free + haemoglobin-bound)."""
    D: float = 2.0e-9           # diffusivity in tissue [m^2/s]
    alpha: float = 1.3e-3       # solubility [mol m^-3 mmHg^-1]
    c_hb: float = 9.0           # carrying capacity of haemoglobin [mol/m^3]
    p50: float = 27.0           # Hill half-saturation [mmHg]
    hill: float = 2.7           # Hill coefficient [-]
    vmax: float = 5.0e-2        # maximal consumption [mol m^-3 s^-1] (~3 umol/min/g)
    km: float = 1.0             # Michaelis constant [mmHg]
    p_in: float = 65.0          # inlet pO2 at the portal triads [mmHg]
    zonation: float = 0.0       # >0 makes consumption periportal (beta-oxidation gradient)


@dataclass
class Ethanol:
    """Ethanol [mol/m^3 = mM] with the two-enzyme kinetics that sets the zonal pattern."""
    D: float = 1.0e-9           # diffusivity [m^2/s]
    vmax_adh: float = 2.4e-2    # ADH capacity [mol m^-3 s^-1]
    km_adh: float = 0.5         # ADH Michaelis constant [mol/m^3]
    vmax_cyp: float = 6.0e-3    # CYP2E1 capacity [mol m^-3 s^-1]
    km_cyp: float = 10.0        # CYP2E1 Michaelis constant [mol/m^3]
    cyp_zonation: float = 3.0   # pericentral sharpness of CYP2E1 [-]
    c_in: float = 20.0          # inlet ethanol [mol/m^3] (~0.9 g/L, legally drunk)


@dataclass
class Acetaldehyde:
    """Acetaldehyde: produced by whatever oxidises ethanol, cleared by ALDH2."""
    D: float = 1.1e-9
    vmax_aldh: float = 6.0e-2   # ALDH2 capacity [mol m^-3 s^-1]
    km_aldh: float = 1.0e-3     # ALDH2 Michaelis constant [mol/m^3] (~1 uM, very avid)
    c_in: float = 0.0


@dataclass
class Drug:
    """A generic small molecule cleared by a pericentral CYP (the paracetamol geometry)."""
    D: float = 0.7e-9
    vmax: float = 8.0e-3        # [mol m^-3 s^-1]
    km: float = 1.0             # [mol/m^3]
    zonation: float = 3.0       # >0 pericentral, <0 periportal, 0 uniform
    c_in: float = 1.0           # inlet concentration [mol/m^3]
    toxic_yield: float = 0.05   # fraction of the flux that becomes a reactive metabolite


@dataclass
class NitricOxide:
    """Wall sub-model. NO is NOT transported on the lobule mesh: with haemoglobin
    scavenging at ~1e5 1/s its diffusion length is well under a micron, so a lobule-scale
    advection-diffusion field for NO is meaningless. See docs/FEM.md."""
    D: float = 3.3e-9           # diffusivity [m^2/s]
    q0: float = 1.0e-8          # basal eNOS flux at the endothelium [mol m^-2 s^-1];
    # only ratios matter, the feedback uses NO / NO_reference
    tau_ref: float = 0.3        # reference wall shear stress [Pa]
    shear_gain: float = 1.0     # eNOS sensitivity to shear [-] (q = q0 (tau/tau_ref)^gain)
    k_tissue: float = 1.0e2     # tissue consumption [1/s]
    k_rbc: float = 1.0e5        # effective scavenging in the RBC core [1/s]
    free_layer: float = 0.5e-6  # RBC-free plasma layer thickness [m]
    wall: float = 0.5e-6        # endothelial layer thickness [m]
    r_target: float = 2.0e-6    # distance outside the wall where NO is read (stellate cells) [m]
    n_radial: int = 240         # nodes of the 1D wall solve


@dataclass
class Feedback:
    """Shear -> NO -> sinusoidal calibre -> permeability. Small gains by default."""
    gain: float = 0.35          # relative permeability change per unit relative NO change
    k_min: float = 0.6          # clamps on the permeability multiplier
    k_max: float = 1.8
    tol: float = 1e-4
    max_iter: int = 40
    relax: float = 0.5
    # Injury branch. The dilation above is adaptive: more shear, more eNOS, wider bed.
    # The porcine small-for-size data say the opposite happens at extreme flow, with
    # portal resistance significantly HIGHER in small grafts, endothelial injury markers
    # rising within minutes of reperfusion and an endothelin-to-NO imbalance. So above
    # some wall shear stress the sign flips. Neither the threshold nor the slope is
    # known, so this is off by default (injury_shear = 0) and exists to let the fork be
    # explored, not to assert a number.
    injury_shear: float = 0.0   # wall shear stress above which injury starts [Pa]
    injury_gain: float = 1.0    # relative calibre loss per unit of relative overshoot


@dataclass
class Params:
    geom: Geometry = None
    flow: Flow = None
    o2: Oxygen = None
    etoh: Ethanol = None
    acet: Acetaldehyde = None
    drug: Drug = None
    no: NitricOxide = None
    fb: Feedback = None

    def __post_init__(self):
        if self.geom is None: self.geom = Geometry()
        if self.flow is None: self.flow = Flow()
        if self.o2 is None: self.o2 = Oxygen()
        if self.etoh is None: self.etoh = Ethanol()
        if self.acet is None: self.acet = Acetaldehyde()
        if self.drug is None: self.drug = Drug()
        if self.no is None: self.no = NitricOxide()
        if self.fb is None: self.fb = Feedback()


def default():
    return Params()
