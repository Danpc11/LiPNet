"""Supplementary figures of the model in a human liver (src/model/human.py; docs/HUMAN.md).

    python src/human_figures.py     -> results/human_liver_values, results/human_liver_growth (.pdf and .png)

human_liver_values   A the portal radii of the cast against the model with no fitted scale; B Poiseuille drop of each measured
                     generation of the cast; C where the portocaval drop sits in the model, in micropuncture and on the cast
human_liver_growth   A portal diameter against liver volume; B lobular drop for growth by adding or enlarging lobules;
                     C graft volume after transplantation in the published series; D portal pressure at one week
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import human

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); OUT = f'{ROOT}/results'
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)     # 'Arial' is a fallback where Liberation Sans is absent
plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'font.size': 7, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.minor.width': 0.35, 'ytick.minor.width': 0.35, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5,
                     'ytick.minor.size': 1.5, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e6e5e0'
BLUE, AQUA, ORANGE, NEUTRAL = '#2a78d6', '#1baf7a', '#eb6834', '#d9d8d2'
mm = 1 / 25.4


def save(fig, name):
    os.makedirs(OUT, exist_ok=True); fig.savefig(f'{OUT}/{name}.pdf'); fig.savefig(f'{OUT}/{name}.png', dpi=300); plt.close(fig); print('saved', name)


def human_liver_values():
    RAMP = {0.5: '#86b6ef', 0.75: '#2a78d6', 1.0: '#0d366b'}
    fig = plt.figure(figsize=(180 * mm, 62 * mm))
    BOX = [(0.062, 0.175, 0.225, 0.70), (0.365, 0.175, 0.175, 0.70), (0.742, 0.175, 0.243, 0.70)]

    # A. absolute portal radii of the eleven measured generations, no fitted scale
    ax = fig.add_axes(BOX[0]); g, r, l, N = human.cast_arrays('portal vein'); x = np.logspace(0, np.log10(6e4), 80)
    for b, col in RAMP.items():
        ax.plot(x, human.radius_at(human.liver(b, liver_g=human.CAST_LIVER_G), x), color=col, lw=1.2, solid_capstyle='round', zorder=2, label=f'Model, b = {b:g}')
    ax.plot(N, r, 'o', ms=4.2, mfc=INK, mec='white', mew=0.7, zorder=3, label='Human cast')
    ax.legend(loc='lower left', frameon=False, fontsize=6.5, handlelength=1.4, labelspacing=0.3, borderaxespad=0.3)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.6, 8e4); ax.set_ylim(0.1, 25)
    ax.set_xticks([1, 10, 100, 1e3, 1e4]); ax.set_xticklabels(['1', '10', '10$^2$', '10$^3$', '10$^4$']); ax.set_yticks([0.1, 1, 10]); ax.set_yticklabels(['0.1', '1', '10'])
    ax.set_xlabel('Vessels in the generation'); ax.set_ylabel('Portal vein radius (mm)')

    # B. Poiseuille on the measured generations
    ax = fig.add_axes(BOX[1])
    for name, col, mk, lab, xy, ha in (('portal vein', BLUE, 'o', 'Portal vein', (9.6, 0.26), 'center'), ('hepatic vein', ORANGE, 's', 'Hepatic veins', (5.9, 1.10), 'right')):
        p = human.cast_poiseuille(name); ax.plot(p['generations'], p['drop_mmHg'], '-', color=col, lw=1.2, zorder=2)
        ax.plot(p['generations'], p['drop_mmHg'], mk, ms=4.2, mfc=col, mec='white', mew=0.7, zorder=3); ax.text(xy[0], xy[1], lab, color=INK, fontsize=6.5, ha=ha, va='bottom')
    ax.set_xlim(0.4, 11.6); ax.set_ylim(0, 1.5); ax.set_xticks([1, 3, 5, 7, 9, 11]); ax.set_yticks([0, 0.5, 1.0, 1.5]); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_xlabel('Generation'); ax.set_ylabel('Pressure drop per generation (mmHg)')

    # C. where the portocaval drop sits
    ax = fig.add_axes(BOX[2]); c1, c0 = human.liver(0.75), human.liver(0.75, lobule_drop=None)
    pv = 0.5 * sum(human.cast_tree_drop('portal vein', G) for G in human.CAST_GENERATIONS['portal vein']) / 5 * 100
    sh = lambda c: [(c['share_inlet'] * 100, BLUE), (c['share_sinusoid'] * 100, AQUA), (c['share_outlet'] * 100, ORANGE)]
    rows = [('Model, lobule imposed', sh(c1)), ('Model, lobule at optimum', sh(c0)), ('Rat (Shibayama 1975)', [(42 / 90 * 100, BLUE), (40 / 90 * 100, AQUA), (8 / 90 * 100, ORANGE)]),
            ('Rat (Bohlen 1991)', [(56, NEUTRAL), (44, ORANGE)]), ('Dog (Bohlen 1991)', [(69, NEUTRAL), (31, ORANGE)]), ('Human cast (Poiseuille)', [(pv, BLUE), (100 - pv, NEUTRAL)])]
    for k, (lab, segs) in enumerate(rows):
        y = len(rows) - 1 - k; left = 0
        for v, col in segs:
            ax.barh(y, v, left=left, height=0.62, color=col, edgecolor='white', linewidth=1.0, zorder=2)
            if v >= 9 and col != NEUTRAL: ax.text(left + v / 2, y, f'{v:.0f}', ha='center', va='center', fontsize=6.5, color='white' if col != AQUA else INK, zorder=3)
            left += v
        ax.text(-2.5, y, lab, ha='right', va='center', fontsize=6.5, color=INK)
    ax.set_xlim(0, 100); ax.set_ylim(-0.6, len(rows) - 0.4); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xticks([0, 25, 50, 75, 100]); ax.set_xlabel('Share of the portocaval pressure drop (%)')
    fig.legend(handles=[Patch(color=BLUE, label='Portal tree'), Patch(color=AQUA, label='Sinusoids'), Patch(color=ORANGE, label='Hepatic veins'), Patch(color=NEUTRAL, label='Not resolved')],
               loc='lower right', bbox_to_anchor=(0.99, 0.885), ncol=4, frameon=False, fontsize=6.5, handlelength=0.9, handleheight=0.8, columnspacing=0.9, handletextpad=0.4, borderaxespad=0)
    for a, letter in zip(fig.axes, 'ABC'): fig.text(a.get_position().x0 - (0.055 if letter != 'C' else 0.158), 0.965, letter, fontsize=9, fontweight='bold', va='top', ha='left', color=INK)
    save(fig, 'human_liver_values')


def human_liver_growth():
    RAMP = {0.6: '#86b6ef', 0.75: '#2a78d6', 0.93: '#0d366b'}
    fig = plt.figure(figsize=(180 * mm, 118 * mm))
    BOX = [(0.075, 0.585, 0.37, 0.355), (0.60, 0.585, 0.37, 0.355), (0.075, 0.095, 0.37, 0.355), (0.72, 0.095, 0.25, 0.355)]
    sizes = np.geomspace(human.NEWBORN_LIVER_G, human.ADULT_LIVER_G, 25)

    # A. portal diameter against liver volume
    ax = fig.add_axes(BOX[0])
    for b, col in RAMP.items():
        g = human.growth(b, 1.0, sizes); ax.plot(sizes / human.DENSITY, [r['portal_diameter_mm'] for r in g['rows']], color=col, lw=1.2, label=f'Model, b = {b:g}', zorder=2)
    V = human.measured_growth()['liver_volume_mL']
    for k, (v, (w, d, vel)) in enumerate(zip(V, human.PORTAL_BY_WEIGHT)):
        ax.plot([v, v], d, color=INK, lw=1.6, solid_capstyle='round', zorder=3, label='Children, range by body weight' if k == 0 else None)
    ax.plot([human.urata_volume(70)], [human.PORTAL_DIAMETER_MM['adults'][0]], 's', ms=4, mfc=INK, mec='white', mew=0.7, zorder=4, label='Adults, pooled ultrasound')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(120, 1700); ax.set_ylim(2.5, 60)
    ax.set_xticks([150, 300, 600, 1200]); ax.set_xticklabels(['150', '300', '600', '1200']); ax.set_yticks([3, 5, 10, 20]); ax.set_yticklabels(['3', '5', '10', '20']); ax.minorticks_off()
    ax.set_xlabel('Liver volume (mL)'); ax.set_ylabel('Portal vein diameter (mm)')
    ax.legend(loc='upper left', frameon=False, fontsize=6.3, handlelength=1.4, labelspacing=0.25, borderaxespad=0.1)

    # B. lobular drop against liver volume
    ax = fig.add_axes(BOX[1])
    for s, ls, lab in ((0.0, '-', 'Lobules added (s = 0)'), (0.5, (0, (4, 1.5)), 's = 0.5'), (1.0, (0, (1, 1.2)), 'Lobules enlarged (s = 1)')):
        g = human.growth(0.75, s, sizes); y = [r['lobule_drop_mmHg'] for r in g['rows']]; ax.plot(sizes / human.DENSITY, y, color=INK, lw=1.2, ls=ls, zorder=2)
        xy = {0.0: (150, 1.04), 0.5: (170, 0.62), 1.0: (330, 0.14)}[s]; ax.text(xy[0], xy[1], lab, fontsize=6.3, color=INK, ha='left', va='bottom')
    ax.set_xscale('log'); ax.set_xlim(120, 1700); ax.set_ylim(0, 1.25); ax.set_xticks([150, 300, 600, 1200]); ax.set_xticklabels(['150', '300', '600', '1200']); ax.minorticks_off()
    ax.set_yticks([0, 0.5, 1.0]); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.set_xlabel('Liver volume (mL)'); ax.set_ylabel('Pressure drop across the lobule (mmHg)')

    # C. graft volume after transplantation, every series
    ax = fig.add_axes(BOX[2]); names = []
    for a, n, t, x in human.GRAFT_VOLUME:
        if a not in names: names.append(a)
    for a in names:
        pts = sorted((t, x, n) for b_, n, t, x in human.GRAFT_VOLUME if b_ == a); t = [p[0] for p in pts]; x = [p[1] for p in pts]
        ax.plot(t, x, '-', color=BLUE, lw=0.8, alpha=0.6, zorder=2); ax.scatter(t, x, s=[6 + 0.12 * p[2] for p in pts], color=BLUE, edgecolor='white', linewidth=0.6, zorder=3)
    ax.axhline(1, color=GRID, lw=0.5, zorder=0); ax.set_xscale('log'); ax.set_xlim(10, 240); ax.set_ylim(0.9, 2.4)
    ax.set_xticks([14, 30, 60, 90, 180]); ax.set_xticklabels(['14', '30', '60', '90', '180']); ax.minorticks_off(); ax.set_yticks([1.0, 1.5, 2.0])
    ax.set_xlabel('Days after transplantation'); ax.set_ylabel('Graft volume / implanted volume')
    for n, yy in ((15, 2.33), (60, 2.2), (300, 2.03)): ax.scatter([150], [yy], s=6 + 0.12 * n, color=BLUE, edgecolor='white', linewidth=0.6, zorder=3); ax.text(162, yy, f'n = {n}', fontsize=6.3, va='center', ha='left', color=INK)

    # D. portal pressure at one week
    ax = fig.add_axes(BOX[3]); w = human.first_week()
    rows = [('Model, calibres unchanged', w['unchanged']['pressure_load'], None, ORANGE), ('Model, calibres grow with the graft', w['proportional']['pressure_load'], None, BLUE)]
    for k, key in (('Measured, portal pressure ≥ 20 mmHg', 'high'), ('Measured, portal pressure < 20 mmHg', 'low')):
        rows.append((k, w['observed_gradient_if_caval_pressure_8'][key], w['observed_portal_pressure'][key], INK))
    for i, (lab, lo, hi, col) in enumerate(rows):
        y = len(rows) - 1 - i
        if hi is None: ax.plot([lo], [y], 'o', ms=5, mfc=col, mec='white', mew=0.7, zorder=3)
        else: ax.plot([lo, hi], [y, y], color=col, lw=2.2, solid_capstyle='round', zorder=3)
        ax.text(-0.04, y, lab, fontsize=6.3, ha='right', va='center', color=INK)
    ax.axvline(1, color=GRID, lw=0.5, zorder=0); ax.set_xlim(0, 1.1); ax.set_ylim(-0.6, len(rows) - 0.4); ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_xticks([0, 0.5, 1.0]); ax.set_xlabel('Portal pressure at one week /\nportal pressure at reperfusion')
    for a, letter, dx in zip(fig.axes, 'ABCD', (0.06, 0.06, 0.06, 0.24)):
        bb = a.get_position(); fig.text(bb.x0 - dx, bb.y1 + 0.045, letter, fontsize=9, fontweight='bold', va='top', ha='left', color=INK)
    save(fig, 'human_liver_growth')


if __name__ == '__main__':
    human_liver_values(); human_liver_growth()
