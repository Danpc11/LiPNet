"""Data panels of the main figures of the article, drawn from the table and the cached results.

    python src/main_figures.py            -> results/main_figure_1DE, main_figure_2, main_figure_3, main_figure_4 (.pdf and .png)

Figure 1D-E   load per lobule against graft size; the flow-pressure plane with the identity zP = zF zR
Figure 2A-B   the six groups that report both loads in that plane; their resistance load with its sampling interval
Figure 3A-C   slope per stratum; outcome against pressure load with the fitted curve of each stratum; specification curve
Figure 4A-B   clinical thresholds on the pressure-load scale; pressure load of the remnant after resection

The schematics of Figure 1A-C are not produced here (network_2d and network_3d of plots.py draw the networks).
Every number comes from data/series.tsv through indices.compute, from audit.py and from results/*.json.
"""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import indices, audit

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); OUT = f'{ROOT}/results'
import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)     # 'Arial' is a fallback where Liberation Sans is absent
plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'font.size': 7, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
INK, MUTED, GRID, BAND = '#0b0b0b', '#52514e', '#e6e5e0', '#f0efec'
BLUE, ORANGE, AQUA, YELLOW, MAGENTA = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4'
MM = 1 / 25.4
MODEL_ZR = (0.81, 1.0)                       # resistance load the network reaches for a partial graft (docs/LIVER.md, section 6)
SHORT = {'Chan 2011': 'Chan 2011', 'Yagi 2006': 'Yagi 2006', 'Ishizaki 2012': 'Ishizaki 2012', 'Wu 2011': 'Wu 2011'}

series = lambda: indices.compute(pd.read_csv(f'{ROOT}/data/series.tsv', sep='\t'))
load = lambda name: json.load(open(f'{OUT}/{name}.json'))


def save(fig, name):
    os.makedirs(OUT, exist_ok=True); fig.savefig(f'{OUT}/{name}.pdf'); fig.savefig(f'{OUT}/{name}.png', dpi=300); plt.close(fig); print('saved', name)


def letter(fig, ax, s, dx=0.055, dy=0.05):
    b = ax.get_position(); fig.text(b.x0 - dx, b.y1 + dy, s, fontsize=9, fontweight='bold', va='top', ha='left', color=INK)


def label_of(r):
    return r.study + (', ' + ('splenectomy' if r.group == 'splenectomy' else 'no splenectomy') if r.study == 'Wang 2014' else '')


def plane(ax, lim=5.0):
    """The flow-pressure plane: donor diagonal, the wedge the network reaches, the window between one and two."""
    x = np.array([0, lim]); ax.fill_between(x, MODEL_ZR[0] * x, MODEL_ZR[1] * x, color=BLUE, alpha=0.16, lw=0, zorder=1)
    ax.plot(x, x, color=INK, lw=0.7, ls=(0, (4, 2)), zorder=2); ax.plot([1], [1], 'o', ms=4, mfc=INK, mec='white', mew=0.6, zorder=5)
    ax.text(4.05, 4.05 * np.mean(MODEL_ZR) - 0.02, 'network model', fontsize=6, rotation=42.5, ha='center', va='center', color=INK, zorder=3)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect('equal'); ax.set_xticks(range(0, int(lim) + 1)); ax.set_yticks(range(0, int(lim) + 1))
    ax.set_xlabel('Flow load, z$_F$'); ax.set_ylabel('Pressure load, z$_P$')


def figure_1():
    fig = plt.figure(figsize=(180 * MM, 68 * MM)); a = fig.add_axes((0.065, 0.17, 0.33, 0.72)); b = fig.add_axes((0.60, 0.17, 0.272, 0.72))
    grwr = np.linspace(0.4, 2.4, 200)                                  # a whole liver is 2% of body weight: g = GRWR / 2
    for h, col in ((1.0, '#86b6ef'), (1.5, BLUE), (2.0, '#0d366b')):
        z = h / (grwr / 2); a.plot(grwr, z, color=col, lw=1.2, zorder=3); a.text(2.44, z[-1], f'h = {h:g}', fontsize=6.5, color=INK, ha='left', va='center', clip_on=False)
    a.axhspan(1, 2, color=BAND, lw=0, zorder=0); a.text(0.46, 1.5, 'window', fontsize=6.5, color=MUTED, ha='left', va='center')
    a.set_xlim(0.4, 2.4); a.set_ylim(0, 7); a.set_xticks([0.5, 1.0, 1.5, 2.0]); a.set_yticks(range(0, 8))
    a.set_xlabel('Graft size, GRWR (%)'); a.set_ylabel('Load per lobule, z = h / g')
    plane(b); b.fill_between([1, 2], 1, 2, color=BAND, lw=0, zorder=0)
    b.text(1.12, 0.72, 'donor', fontsize=6.5, color=INK); b.text(2.9, 3.2, 'z$_R$ = 1', fontsize=6.5, rotation=45, ha='center', va='center', color=INK)
    b.text(0.25, 4.3, 'z$_R$ > 1\npressure above flow', fontsize=6.5, color=INK, ha='left', va='center')
    b.text(4.85, 1.2, 'z$_R$ < 1\npressure below flow', fontsize=6.5, color=INK, ha='right', va='center')
    letter(fig, a, 'D', 0.055, 0.07); letter(fig, b, 'E', 0.06, 0.07); save(fig, 'main_figure_1DE')


def figure_2():
    d = series(); p = d.dropna(subset=['zF', 'zP']); flow_only = d[d.zF.notna() & d.zP.isna() & d.in_main_analysis]
    A = {(g['study'], g['group']): g for g in audit.paired(d)['groups']}
    fig = plt.figure(figsize=(180 * MM, 76 * MM)); a = fig.add_axes((0.065, 0.16, 0.36, 0.74)); b = fig.add_axes((0.70, 0.16, 0.285, 0.74))
    plane(a); a.fill_between([1, 2], 1, 2, color=BAND, lw=0, zorder=0)
    a.plot(flow_only.zF, np.zeros(len(flow_only)), '|', ms=7, mew=0.9, color=MUTED, clip_on=False, zorder=4)
    order = {(r.study, r.group): k + 1 for k, r in enumerate(sorted(p.itertuples(), key=lambda r: r.zR))}; left = {'Yagi 2006', 'Wu 2011', 'Ishizaki 2012'}
    for _, r in p.iterrows():
        std = not r.outflow.startswith('reconstructed'); ms = 3.2 + r.n ** 0.5 / 3.2; dx = (0.09 + ms / 95) * (-1 if r.study in left else 1)
        a.plot(r.zF, r.zP, 'o', ms=ms, mfc='white' if std else BLUE, mec=BLUE, mew=1.0, zorder=6)
        a.text(r.zF + dx, r.zP, str(order[(r.study, r.group)]), fontsize=6.5, color=INK, ha='right' if r.study in left else 'left', va='center', zorder=7)
    a.text(0.72, 1.28, 'donor', fontsize=6.5, color=INK, ha='right'); a.text(2.9, 3.2, 'z$_R$ = 1', fontsize=6.5, rotation=45, ha='center', va='center')
    a.text(2.9, 0.22, 'flow load only', fontsize=6.5, color=MUTED, ha='center', va='bottom')
    a.legend(handles=[Line2D([], [], marker='o', ls='', ms=4.5, mfc=BLUE, mec=BLUE, label='Reconstructed outflow'), Line2D([], [], marker='o', ls='', ms=4.5, mfc='white', mec=BLUE, mew=1.0, label='Standard outflow')],
             loc='upper left', fontsize=6.5, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.3)
    rows = sorted(p.itertuples(), key=lambda r: r.zR); b.fill_between(MODEL_ZR, -0.6, len(rows) - 0.35, color=BLUE, alpha=0.16, lw=0, zorder=0); b.plot([1, 1], [-0.6, len(rows) - 0.35], color=INK, lw=0.7, ls=(0, (4, 2)), zorder=1)
    for k, r in enumerate(rows):
        y = len(rows) - 1 - k; g = A[(r.study, r.group)]; lo, hi = g['zR_ci_sampling']; std = not r.outflow.startswith('reconstructed')
        b.plot([lo, hi], [y, y], color=BLUE, lw=1.2, solid_capstyle='round', zorder=2); b.plot(r.zR, y, 'o', ms=5, mfc='white' if std else BLUE, mec=BLUE, mew=1.0, zorder=3)
        b.text(-0.04, y + 0.13, f'{k + 1}  ' + label_of(r), fontsize=6.5, color=INK, ha='right', va='center'); b.text(r.zR, y + 0.2, f'{r.zR:.2f}', fontsize=6.5, color=INK, ha='center', va='bottom')
        b.text(-0.04, y - 0.2, f'n = {int(r.n)}, z$_F$ = {r.zF:.2f}', fontsize=6, color=MUTED, ha='right', va='center')
    b.text(np.mean(MODEL_ZR), len(rows) - 0.28, 'network model', fontsize=6.5, color=INK, ha='center', va='bottom')
    b.set_xlim(0, 1.2); b.set_ylim(-0.6, len(rows) + 0.25); b.set_yticks([]); b.spines['left'].set_visible(False); b.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]); b.set_xlabel('Resistance load, z$_R$ = z$_P$ / z$_F$')
    letter(fig, a, 'A', 0.05, 0.06); letter(fig, b, 'B', 0.245, 0.06); save(fig, 'main_figure_2')


def figure_3():
    E = load('bayes_main')['exploratory (all outcomes)']; S = load('spec_curve'); d = series(); strata = sorted(E['strata'], key=lambda s: E['beta_by_stratum'][s]['mean'])
    col = dict(zip(sorted(E['strata']), (BLUE, ORANGE, AQUA, YELLOW, MAGENTA))); name = lambda s: s.split(' / ')[0]
    fig = plt.figure(figsize=(180 * MM, 68 * MM)); a = fig.add_axes((0.105, 0.17, 0.17, 0.72)); b = fig.add_axes((0.36, 0.17, 0.27, 0.72)); c = fig.add_axes((0.845, 0.17, 0.14, 0.72))
    a.axvspan(*E['mu_beta_ci'], color=BAND, lw=0, zorder=0); a.axvline(E['mu_beta'], color=INK, lw=0.9, zorder=1); a.axvline(0, color=MUTED, lw=0.5, ls=(0, (4, 2)))
    for k, s in enumerate(strata):
        v = E['beta_by_stratum'][s]; a.plot(v['ci'], [k, k], color=col[s], lw=1.2, solid_capstyle='round', zorder=2); a.plot(v['mean'], k, 'o', ms=4.6, mfc=col[s], mec='white', mew=0.6, zorder=3)
        a.text(-0.62, k, name(s), fontsize=6.5, color=INK, ha='right', va='center')
    a.text(E['mu_beta'] + 0.1, len(strata) - 0.42, f"$\\mu_\\beta$ = {E['mu_beta']:+.2f}", fontsize=6.5, ha='left', va='bottom', color=INK)
    a.set_xlim(-0.5, 3.5); a.set_ylim(-0.6, len(strata) + 0.1); a.set_yticks([]); a.spines['left'].set_visible(False); a.set_xticks([0, 1, 2, 3]); a.set_xlabel('Slope per unit of z$_P$')
    m = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]; m = m.assign(stratum=m.study + ' / ' + m.outcome_type); m = m[m.stratum.isin(E['strata'])]
    b.axvspan(1, 2, color=BAND, lw=0, zorder=0)
    for s in E['strata']:
        g = m[m.stratum == s]; z = np.linspace(max(0.8, g.zP.min() - 0.35), min(4.0, g.zP.max() + 0.35), 60)
        b.plot(z, 100 / (1 + np.exp(-(E['alphas'][s] + E['beta_by_stratum'][s]['mean'] * z))), color=col[s], lw=1.2, zorder=2)
        b.scatter(g.zP, g.pct, s=10 + g.n / 5, color=col[s], edgecolor='white', linewidth=0.6, zorder=3)
    b.set_xlim(0.8, 4.0); b.set_ylim(0, 60); b.set_xticks([1, 2, 3, 4]); b.set_xlabel('Pressure load, z$_P$'); b.set_ylabel('Outcome (%)')
    b.legend(handles=[Line2D([], [], marker='o', ls='-', lw=1.2, ms=4, color=col[s], mec='white', label=name(s)) for s in sorted(E['strata'])], loc='upper left', fontsize=6.5, labelspacing=0.3, borderaxespad=0.2, handlelength=1.6)
    S = sorted(S, key=lambda s: s['mu_beta']); c.axvline(0, color=MUTED, lw=0.5, ls=(0, (4, 2)))
    short = {'primary outcome (SFSS or early dysfunction)': 'primary outcome', 'secondary outcome (mortality or graft loss)': 'secondary outcome', 'gradients reported or derived only': 'gradients reported or derived',
             'events reported only (no counts from percentages)': 'events reported only', 'one Kyoto publication (Ogura 2010)': 'Kyoto: Ogura 2010 only', 'one Kyoto publication (Uemura 2016)': 'Kyoto: Uemura 2016 only',
             'one Kyoto publication (Yagi 2005)': 'Kyoto: Yagi 2005 only', 'tau prior HalfNormal(0.25)': '$\\tau$ prior HN(0.25)', 'tau prior HalfNormal(1.0)': '$\\tau$ prior HN(1.0)',
             'one stratum per centre, not per publication': 'one stratum per centre', 'assumed CVP 3 mmHg': 'CVP 3 mmHg', 'assumed CVP 7 mmHg': 'CVP 7 mmHg', 'assumed CVP 9 mmHg': 'CVP 9 mmHg'}
    for k, s in enumerate(S):
        c.plot(s['mu_beta_ci'], [k, k], color=BLUE, lw=1.0, solid_capstyle='round'); c.plot(s['mu_beta'], k, 'o', ms=3.6, mfc=BLUE, mec='white', mew=0.5)
        c.text(-0.72, k, short.get(s['spec'], s['spec']), fontsize=6, color=INK, ha='right', va='center')
    c.set_xlim(-0.6, 3.6); c.set_ylim(-0.7, len(S) - 0.3); c.set_yticks([]); c.spines['left'].set_visible(False); c.set_xticks([0, 1, 2, 3]); c.set_xlabel('$\\mu_\\beta$ across specifications')
    letter(fig, a, 'A', 0.095, 0.07); letter(fig, b, 'B', 0.05, 0.07); letter(fig, c, 'C', 0.20, 0.07); save(fig, 'main_figure_3')


def figure_4():
    P = load('predictions')['P4_thresholds']; ref = P['reference_sensitivity']
    rows = [('Clinically significant portal hypertension', 'HVPG 10 mmHg', 'cirrhosis_HVPG_10'), ('Partial graft', 'PVP 15 mmHg at CVP 5 mmHg', 'graft_PVP15_CVP5'),
            ('Resection limit in cirrhosis', 'HVPG 10 mmHg', 'resection_HVPG_10'), ('Variceal bleeding', 'HVPG 12 mmHg', 'variceal_HVPG_12')]
    fig = plt.figure(figsize=(180 * MM, 64 * MM)); a = fig.add_axes((0.285, 0.18, 0.19, 0.72)); b = fig.add_axes((0.60, 0.18, 0.385, 0.72))
    a.axvspan(1, 2, color=BAND, lw=0, zorder=0)
    for k, (lab, sub, key) in enumerate(rows):
        y = len(rows) - 1 - k; lo, mid, hi = ref['6 mmHg reference'][key], ref['5 mmHg reference'][key], ref['4 mmHg reference'][key]
        a.plot([lo, hi], [y, y], color=BLUE, lw=1.2, solid_capstyle='round', zorder=2); a.plot(mid, y, 'o', ms=5, mfc=BLUE, mec='white', mew=0.6, zorder=3)
        a.text(mid, y + 0.2, f'{mid:.1f}', fontsize=6.5, ha='center', va='bottom', color=INK)
        a.text(-0.12, y + 0.12, lab, fontsize=6.5, color=INK, ha='right', va='center'); a.text(-0.12, y - 0.2, sub, fontsize=6, color=MUTED, ha='right', va='center')
    a.set_xlim(0, 3.5); a.set_ylim(-0.6, len(rows) - 0.3); a.set_yticks([]); a.spines['left'].set_visible(False); a.set_xticks([0, 1, 2, 3]); a.set_xlabel('Pressure load, z$_P$')
    removed = np.linspace(0, 75, 200); g = 1 - removed / 100; b.axhspan(1, 2, color=BAND, lw=0, zorder=0); b.axvline(60, color=MUTED, lw=0.5, ls=(0, (4, 2)))
    for hv, colr in ((5, '#86b6ef'), (8, '#5598e7'), (10, BLUE), (12, '#0d366b')):
        z = hv / indices.NORMAL_GRADIENT / g; b.plot(removed, z, color=colr, lw=1.2, zorder=3); k = np.argmin(abs(z - 4.75)) if z.max() > 4.75 else len(z) - 1
        b.text(removed[k] - 1.0 if z.max() > 4.75 else 75.8, min(z[k], 4.75) if z.max() > 4.75 else z[k], f'{hv}', fontsize=6.5, color=INK, ha='right' if z.max() > 4.75 else 'left', va='bottom' if z.max() > 4.75 else 'center')
    b.text(3, 4.7, 'Baseline HVPG (mmHg)', fontsize=6.5, color=INK, ha='left', va='center'); b.text(59, 0.35, 'right hepatectomy', fontsize=6.5, color=MUTED, ha='right', va='center')
    b.set_xlim(0, 80); b.set_ylim(0, 5); b.set_xticks([0, 20, 40, 60, 80]); b.set_xlabel('Parenchyma removed (%)'); b.set_ylabel('z$_P$ of the remnant')
    letter(fig, a, 'A', 0.275, 0.06); letter(fig, b, 'B', 0.05, 0.06); save(fig, 'main_figure_4')


if __name__ == '__main__':
    figure_1(); figure_2(); figure_3(); figure_4()
