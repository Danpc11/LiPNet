"""Whole-liver transplantation in the same model: what the anastomoses do.

A whole graft keeps the donor parenchyma intact, so g = 1 and the flow per lobule is set only by the recipient's
splanchnic inflow, zF = h. What can still change is the two anastomoses, which the model adds as series
resistances around the lobular bed:

    portal trunk --[R_in]-- graft portal tree --[R_lob]-- central veins --[R_out]-- cava (CVP)

With Poiseuille conductance a narrowing to a fraction d of the original diameter multiplies that resistance by
d^-4: 0.8 gives 2.4x, 0.7 gives 4.2x, 0.5 gives 16x. Write the baseline shares of the total graft resistance as
r_in and r_out (a few per cent each in a normal anastomosis) and r_lob = 1 - r_in - r_out.

The splanchnic bed is not a perfect flow source: part of the inflow can escape through collaterals, so the graft
receives F = dP / (R_in + R_lob + R_out) in parallel with a collateral path R_coll. That is what makes inlet and
outlet stenosis behave differently, which is the prediction:

  * inlet (portal anastomosis) stenosis   -> flow falls and sinusoidal pressure falls with it. The graft is
    protected, the recipient is not: pressure rises upstream, in the portal trunk, and flow can be stolen by
    collaterals. zP/zF stays at or below 1 while both indices drop.
  * outlet (hepatic vein or caval) stenosis -> flow is maintained by the driving pressure but every lobule now
    drains against a higher pedestal, so sinusoidal pressure rises at nearly constant flow. zP/zF climbs above 1.

So the two failures are separable at the bedside with the two numbers already measured, and they are separable in
opposite directions, which a single threshold on flow or on pressure cannot do.
"""
import numpy as np


def whole_graft(d_in=1.0, d_out=1.0, h=1.0, r_in=0.05, r_out=0.05, r_coll=12.0, measure_pvp_after_anastomosis=True):
    """Indices of a whole graft whose inlet and outlet are narrowed to fractions d_in and d_out of their diameter.

    h          recipient splanchnic inflow relative to the donor's (1 = normal, 1.5-2 = hyperdynamic cirrhosis)
    r_in,r_out share of the total graft resistance taken by each anastomosis when it is not narrowed
    r_coll     resistance of the portosystemic collateral path, in units of the normal graft resistance
    measure_pvp_after_anastomosis  True if PVP is measured in the graft portal vein (distal to the anastomosis),
               which is the usual site; False if measured in the recipient portal trunk, proximal to it.

    Returns the flow per lobule (zF), the sinusoidal pressure relative to normal (zP), their ratio, the fraction
    of the inflow that escapes through collaterals, and the pressure in the recipient portal trunk.
    """
    r_lob = 1.0 - r_in - r_out
    R_in, R_out = r_in * d_in ** -4, r_out * d_out ** -4
    R = R_in + r_lob + R_out                                   # graft path, in units of the normal graft resistance
    # The splanchnic bed offers h units of flow at a pressure set by the graft and the collaterals in parallel.
    # Everything is normalised to the same liver with normal anastomoses, so that d_in = d_out = 1 gives zF = h
    # exactly, whatever collateral resistance is assumed; the collaterals then only act through what a rising
    # graft resistance diverts.
    par = lambda RR: 1 / (1 / RR + 1 / r_coll)
    f0 = h * par(1.0) / 1.0                                    # flow with normal anastomoses, the reference
    dP_trunk = h * par(R)
    f_abs = dP_trunk / R
    f = h * f_abs / f0                                         # flow per lobule, relative to the donor
    dP_trunk = dP_trunk / (h * par(1.0)) * h                   # trunk pressure on the same scale
    steal = 1 - f_abs / f0
    # pressure that the lobules actually see: their own drop plus the pedestal the outlet imposes
    zP_true = f * (r_lob + R_out) / (r_lob + r_out)     # the lobule's own drop plus the pedestal from the outlet
    zP_meas = zP_true if measure_pvp_after_anastomosis else dP_trunk
    return dict(zF=float(f), zP=float(zP_meas), zP_over_zF=float(zP_meas / f) if f > 0 else np.nan,
                collateral_steal=float(steal), portal_trunk_pressure=float(dP_trunk),
                R_in=float(R_in), R_out=float(R_out), R_total=float(R))


def sweep(which='out', d=np.linspace(1.0, 0.4, 25), **kw):
    """Indices along a narrowing of the inlet ('in') or the outlet ('out')."""
    rows = []
    for x in d:
        r = whole_graft(**{('d_in' if which == 'in' else 'd_out'): float(x)}, **kw)
        r['diameter_ratio'] = float(x); r['narrowed'] = which; rows.append(r)
    return rows


if __name__ == '__main__':
    print('whole graft, normal anastomoses, recipient inflow 1.5x donor:')
    print('  ', {k: round(v, 3) for k, v in whole_graft(h=1.5).items()})
    for which, lab in (('out', 'outlet (hepatic vein / caval)'), ('in', 'inlet (portal anastomosis)')):
        print(f'\n{lab} narrowed, recipient inflow 1.5x donor:')
        print('   d     zF     zP    zP/zF   steal   trunk P')
        for r in sweep(which, d=[1.0, 0.9, 0.8, 0.7, 0.6, 0.5], h=1.5):
            print(f"  {r['diameter_ratio']:.1f}  {r['zF']:5.2f}  {r['zP']:5.2f}  {r['zP_over_zF']:5.2f}  "
                  f"{100*r['collateral_steal']:5.1f}%  {r['portal_trunk_pressure']:5.2f}")
