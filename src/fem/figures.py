"""Figures for the lobule FEM. Same convention as src/plots.py: one function per figure,
results/<name>.pdf and .png.

    python src/fem/run_fem.py figures
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'font.size': 8,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
ROOT = os.path.join(os.path.dirname(__file__), '..', '..')
OUT = os.path.join(ROOT, 'results')


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved', name)


def _tri(mesh):
    return Triangulation(mesh.xy[:, 0] * 1e6, mesh.xy[:, 1] * 1e6, mesh.tri)


def _map(ax, mesh, c, title, cmap, elementwise=False, units=''):
    t = _tri(mesh)
    if elementwise:
        h = ax.tripcolor(t, facecolors=c, cmap=cmap, shading='flat')
    else:
        h = ax.tripcolor(t, c, cmap=cmap, shading='gouraud')
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=8)
    cb = plt.colorbar(h, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    cb.set_label(units, fontsize=6)
    return h


def fields(mesh, sol, o2, et, ac, dr, no_c, name='fem_lobule_fields'):
    """One page: pressure, speed, shear, pO2, ethanol, acetaldehyde, reactive metabolite, NO."""
    from .params import to_mmHg
    fig, axes = plt.subplots(2, 4, figsize=(11, 5.6))
    _map(axes[0, 0], mesh, to_mmHg(sol.p), 'pressure', 'viridis', units='mmHg')
    _map(axes[0, 1], mesh, sol.interstitial * 1e6, 'sinusoidal velocity', 'magma',
         elementwise=True, units='um/s')
    _map(axes[0, 2], mesh, sol.wall_shear, 'wall shear stress', 'inferno',
         elementwise=True, units='Pa')
    _map(axes[0, 3], mesh, no_c * 1e6, 'perivascular NO', 'cividis',
         elementwise=True, units='nM')
    _map(axes[1, 0], mesh, o2.c, 'pO2', 'RdYlBu_r', units='mmHg')
    _map(axes[1, 1], mesh, et.c, 'ethanol', 'YlGnBu', units='mM')
    _map(axes[1, 2], mesh, ac.c * 1e3, 'acetaldehyde', 'YlOrBr', units='uM')
    _map(axes[1, 3], mesh, dr.reactive * 1e6, 'reactive metabolite', 'Reds',
         units='umol m-3 s-1')
    fig.suptitle('Hepatic lobule, finite elements: perfusion and transport', fontsize=9)
    fig.tight_layout()
    save(fig, name)


def profiles(mesh, o2, et, ac, dr, name='fem_zonal_profiles'):
    """Everything against the zonation coordinate, periportal (0) to pericentral (1)."""
    fig, axes = plt.subplots(1, 4, figsize=(10, 2.4))
    for ax, (res, lab, scale, unit) in zip(axes, [(o2, 'pO2', 1.0, 'mmHg'),
                                                  (et, 'ethanol', 1.0, 'mM'),
                                                  (ac, 'acetaldehyde', 1e3, 'uM'),
                                                  (dr, 'drug', 1.0, 'mM')]):
        z, p = mesh.radial_profile(res.c * scale)
        ax.plot(z, p, color='#1f5fbf', lw=1.4)
        ax.set_xlabel('zeta   (0 periportal, 1 pericentral)')
        ax.set_ylabel('%s [%s]' % (lab, unit))
        ax.set_xlim(0, 1)
    axes[0].axhspan(0, 10, color='#c62828', alpha=0.12)
    axes[0].text(0.03, 11, 'hypoxic', color='#c62828', fontsize=6)
    fig.tight_layout()
    save(fig, name)


def resistance_curve(load, name='fem_resistance_load'):
    """The lobule's own z_R against the flow load z_F: what the shear-NO loop can
    contribute to the resistance load that LiPNet measures clinically."""
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    zf, zr = np.array(load['z_F']), np.array(load['z_R'])
    ax.axhline(1.0, color='0.6', lw=0.8, ls='--')
    ax.plot(zf, zr, 'o-', color='#c62828', lw=1.4, ms=3.5)
    ax.plot(zf, np.ones_like(zf), '-', color='0.6', lw=0.8)
    ax.set_xlabel('flow load  $z_F$')
    ax.set_ylabel('lobule resistance load  $z_R$')
    ax.set_title('shear-NO autoregulation, gain %.2f' % load['no_feedback_gain'], fontsize=8)
    ax.text(zf[-1], 1.02, 'rigid Poiseuille edge', color='0.5', fontsize=6, ha='right')
    fig.tight_layout()
    save(fig, name)


def no_profile(no, flow, name='fem_no_wall'):
    """The NO sub-model: why it cannot live on the lobule mesh."""
    from . import no_wall
    rc, cq = no_wall.radial_profile(no, flow)
    a = 0.5 * flow.d_sin
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.semilogy(rc * 1e6, np.maximum(cq, 1e-12) / max(cq.max(), 1e-30), color='#2e7d32', lw=1.4)
    ax.axvspan(0, (a - no.free_layer) * 1e6, color='#c62828', alpha=0.10)
    ax.axvspan((a - no.free_layer) * 1e6, a * 1e6, color='0.85')
    ax.axvspan(a * 1e6, (a + no.wall) * 1e6, color='#1f5fbf', alpha=0.18)
    ax.text(0.2, 0.5, 'red cells', fontsize=6, color='#c62828')
    ax.set_xlabel('radius from the sinusoid axis [um]')
    ax.set_ylabel('[NO] / peak')
    ax.set_title('NO around one sinusoid', fontsize=8)
    fig.tight_layout()
    save(fig, name)


def exchange_sweep(sweep, single_phase_extraction=None, name='fem_two_phase'):
    """Clearance against the exchange Damkohler number: where the membrane takes over
    from the enzyme."""
    da = np.array(sweep['Da_exchange'])
    ext = np.array(sweep['extraction'])
    ratio = np.array(sweep['hepatocyte_over_blood'])
    fig, ax = plt.subplots(figsize=(3.6, 2.6))
    ax.axvline(1.0, color='0.6', lw=0.8, ls='--')
    ax.semilogx(da, ext, 'o-', color='#1f5fbf', lw=1.4, ms=3.5, label='extraction')
    ax.semilogx(da, ratio, 's-', color='#2e7d32', lw=1.2, ms=3,
                label='hepatocyte / blood')
    if single_phase_extraction is not None:
        ax.axhline(single_phase_extraction, color='#c62828', lw=0.9, ls=':',
                   label='single continuum')
    ax.set_xlabel(r'exchange Damkohler number  $\beta\,\tau_{transit}$')
    ax.set_ylabel('fraction')
    ax.text(1.15, 0.05, 'membrane  |  enzyme', color='0.5', fontsize=6)
    ax.legend(fontsize=6, loc='center left')
    fig.tight_layout()
    save(fig, name)
