"""Command line for the lobule FEM. Run from the repository root:

    python src/fem/run_fem.py perfusion      # Darcy solve, resistance, transit, shear
    python src/fem/run_fem.py transport      # O2, ethanol, acetaldehyde, drug
    python src/fem/run_fem.py no             # the NO wall sub-model
    python src/fem/run_fem.py twophase       # sinusoid/hepatocyte exchange, membrane vs enzyme
    python src/fem/run_fem.py coupling       # shape factor and z_R(z_F) for the network
    python src/fem/run_fem.py network        # how much of it survives the rest of the series
    python src/fem/run_fem.py figures        # every figure
    python src/fem/run_fem.py all

Numbers go to results/fem_*.json, figures to results/fem_*.pdf and .png.
Options: --nr --nt (mesh), --inlet perimeter|triads, --gain (NO feedback), --quiet.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import numpy as np
from fem import params as P
from fem import mesh as mesh_mod, darcy, transport as tr, no_wall, coupling, two_phase as tp
from fem import network_link as nl

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
OUT = os.path.join(ROOT, 'results')


def dump(name, obj, quiet=False):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name + '.json')
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=True)
    if not quiet:
        print(json.dumps(obj, indent=1, sort_keys=True))
    print('wrote', os.path.relpath(path, ROOT))


def build(args):
    p = P.default()
    if args.nr:
        p.geom.nr = args.nr
    if args.nt:
        p.geom.nt = args.nt
    if args.inlet:
        p.geom.inlet = args.inlet
    if args.gain is not None:
        p.fb.gain = args.gain
    m = mesh_mod.build(p.geom)
    return p, m


def cmd_perfusion(args):
    p, m = build(args)
    sol = darcy.solve(m, p.flow)
    out = sol.summary()
    out['mass_balance_error'] = sol.mass_error
    out['nodes'] = m.n
    out['elements'] = int(len(m.area))
    out['inlet'] = m.inlet_mode
    out['lobule_volume_mm3'] = m.volume() * 1e9
    out['perfusion_mL_min_per_g'] = sol.Q * 6e7 / (m.volume() * 1e6)
    dump('fem_perfusion', out, args.quiet)
    return p, m, sol


def cmd_transport(args):
    p, m, sol = cmd_perfusion(args)
    o2 = tr.oxygen(m, sol, p.o2)
    et = tr.ethanol(m, sol, p.etoh)
    ac = tr.acetaldehyde(m, sol, et, p.acet)
    dr = tr.drug(m, sol, p.drug)
    out = {}
    for r in (o2, et, ac, dr):
        out.update(r.summary())
    out['pO2_saturation_in'] = o2.sat_in
    out['pO2_saturation_out'] = o2.sat_out
    out['ethanol_cyp2e1_fraction'] = et.cyp_fraction
    out['acetaldehyde_peak_zeta'] = ac.peak_zeta
    out['drug_reactive_peak_zeta'] = dr.reactive_peak_zeta
    dump('fem_transport', out, args.quiet)
    return p, m, sol, o2, et, ac, dr


def cmd_no(args):
    p, m = build(args)
    sol = darcy.solve(m, p.flow)
    out = no_wall.summary(m, sol, p.no)
    rc, cq = no_wall.radial_profile(p.no, p.flow)
    a = 0.5 * p.flow.d_sin
    out['no_lumen_over_wall'] = float(np.interp(0.0, rc, cq) / max(np.interp(a, rc, cq), 1e-300))
    out['no_diffusion_length_tissue_um'] = float(np.sqrt(p.no.D / p.no.k_tissue) * 1e6)
    out['no_diffusion_length_lumen_um'] = float(np.sqrt(p.no.D / p.no.k_rbc) * 1e6)
    dump('fem_no', out, args.quiet)
    return out


def cmd_twophase(args):
    """Sinusoid and hepatocyte as two continua: how much of the clearance is the enzyme
    and how much is the membrane."""
    p, m = build(args)
    sol = darcy.solve(m, p.flow)
    out = tp.permeability_sweep(m, sol, p.drug,
                                [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])
    o2 = tp.oxygen(m, sol, p.o2)
    out.update(o2.summary())
    out['zonation_matrix'] = tp.zonation_matrix(m, sol, p.drug, 1e-7)
    out['oxygen_single_phase_extraction'] = tr.oxygen(m, sol, p.o2).extraction
    out['drug_single_phase_extraction'] = tr.drug(m, sol, p.drug).extraction
    dump('fem_two_phase', out, args.quiet)
    return out


def cmd_network(args):
    """How much of the lobule's z_R survives the rest of the series."""
    p, _ = build(args)
    sc = nl.share_scaling((3.0, 4.5, 6.5))
    load = coupling.resistance_load(p)
    out = dict(share_by_size=[(r['n_lobules'], r['sinusoidal_share']) for r in sc['rows']],
               share_log_slope=sc['log_slope'],
               share_extrapolated_1e6_lobules=sc['extrapolated_at_1e6'])
    out['at_model_share'] = nl.sweep(load, sc['rows'][-1]['sinusoidal_share'])
    out['at_extrapolated_share'] = nl.sweep(load, sc['extrapolated_at_1e6'])
    out['if_sinusoids_carried_half'] = nl.sweep(load, 0.5)
    try:
        out['against_the_published_groups'] = nl.predict_series(
            load, sc['rows'][-1]['sinusoidal_share'])
    except Exception as exc:                       # pandas or the series file missing
        out['against_the_published_groups'] = 'unavailable: %s' % exc
    dump('fem_network_link', out, args.quiet)
    return out


def cmd_coupling(args):
    p, _ = build(args)
    out = coupling.shape_factor(p)
    out.update(coupling.resistance_load(p))
    dump('fem_coupling', out, args.quiet)
    return out


def cmd_figures(args):
    from fem import figures
    p, m, sol, o2, et, ac, dr = cmd_transport(args)
    no_c, _ = no_wall.field(m, sol, p.no)
    figures.fields(m, sol, o2, et, ac, dr, no_c)
    figures.profiles(m, o2, et, ac, dr)
    figures.no_profile(p.no, p.flow)
    sweep = tp.permeability_sweep(m, sol, p.drug, [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8])
    figures.exchange_sweep(sweep, tr.drug(m, sol, p.drug).extraction)
    path = os.path.join(OUT, 'fem_coupling.json')
    load = json.load(open(path)) if os.path.exists(path) else coupling.resistance_load(p)
    figures.resistance_curve(load)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('what', nargs='*', default=['all'])
    ap.add_argument('--nr', type=int)
    ap.add_argument('--nt', type=int)
    ap.add_argument('--inlet', choices=['perimeter', 'triads'])
    ap.add_argument('--gain', type=float)
    ap.add_argument('--quiet', action='store_true')
    args = ap.parse_args(argv)
    what = args.what
    if 'all' in what:
        what = ['perfusion', 'transport', 'no', 'twophase', 'coupling', 'network',
                'figures']
    jobs = dict(perfusion=cmd_perfusion, transport=cmd_transport, no=cmd_no,
                twophase=cmd_twophase, coupling=cmd_coupling, network=cmd_network,
                figures=cmd_figures)
    for w in what:
        if w not in jobs:
            raise SystemExit('unknown target %r; choose from %s' % (w, ', '.join(sorted(jobs))))
        print('\n=== %s ===' % w)
        jobs[w](args)


if __name__ == '__main__':
    main()
