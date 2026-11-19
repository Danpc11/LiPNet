"""Every figure of the article, in one style (src/figstyle.py, docs/FIGURES.md).

    python src/paper_figures.py              all ten figures -> results/Figure_1 ... Figure_4, Figure_S1 ... Figure_S6 (.pdf and .png)
    python src/paper_figures.py 2 S3         only those

Figure 1    A the model: one lobule and the network; B, C the optimum in two and three dimensions; D load per lobule against
            graft size; E the flow-pressure plane with the identity zP = zF zR
Figure 2    A the six groups that report both loads in that plane; B their resistance load with its sampling interval (thick)
            and the interval that also carries the uncertainty of the references (thin)
Figure 3    A slope per stratum; B outcome against pressure load with the fitted curve of each stratum; C specification curve
Figure 4    A clinical thresholds on the pressure-load scale; B pressure load of the remnant after resection
Figure S1   A scaling of the optimum with the number of lobules; B, C allometric closure
Figure S2   inlet and outlet obstruction of a whole graft
Figure S3   slope for a patient against the within-group dispersion of the gradient
Figure S4   each recipient of a hemiportocaval shunt
Figure S5   the model against a human vascular cast
Figure S6   growth with the body and regrowth of a graft

Every number comes from data/, from results/*.json or from the model itself; nothing is typed in here except labels.
Figures 1B and 1C recompute the optimum (about a minute). Figures 2, 3, 4 and S3 need the results of run_all.sh.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs
from figstyle import INK, MUTED, GRID, BAND, NEUTRAL, BLUE, PORTAL, SINUSOID, HEPATIC, EVENT, FS_IN, FS_SMALL, LW, LW_REF, DASH, DOT
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch, FancyBboxPatch, RegularPolygon
from matplotlib.path import Path
import indices

fs.use()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); OUT = f'{ROOT}/results'; DATA = f'{ROOT}/data'
MODEL_ZR = (0.81, 1.0)                       # resistance load the network reaches for a partial graft (docs/LIVER.md, section 6)
ZF, ZP, ZR = '$z_F$', '$z_P$', '$z_R$'

series = lambda: indices.compute(pd.read_csv(f'{DATA}/series.tsv', sep='\t'))
load = lambda name: json.load(open(f'{OUT}/{name}.json'))
save = lambda fig, name: fs.save(fig, name, OUT)


def plane(ax, lim=5.0, model=True, label_at=2.9):
    """The flow-pressure plane: donor diagonal, the window between one and two, the wedge the network reaches."""
    x = np.array([0, lim]); ax.fill_between([1, 2], 1, 2, color=BAND, lw=0, zorder=0)
    if model:
        ax.fill_between(x, MODEL_ZR[0] * x, MODEL_ZR[1] * x, zorder=1, **fs.MODEL_TINT)
        ax.text(0.81 * lim, 0.81 * lim * np.mean(MODEL_ZR) - 0.004 * lim, 'Network model', fontsize=FS_SMALL, rotation=42.5, ha='center', va='center', zorder=3)
    ax.plot(x, x, color=INK, lw=LW_REF, ls=DASH, zorder=2); ax.plot([1], [1], 'o', ms=4, mfc=INK, zorder=5, **DOT)
    ax.text(label_at, label_at + 0.06 * lim, f'{ZR} = 1', fontsize=FS_IN, rotation=45, ha='center', va='center')
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect('equal'); ax.set_xticks(range(0, int(lim) + 1)); ax.set_yticks(range(0, int(lim) + 1))
    ax.set_xlabel(f'Flow load, {ZF}'); ax.set_ylabel(f'Pressure load, {ZP}')


# ------------------------------------------------------------------ Figure 1
def _curve(ax, p0, p1, lw, color, z=2):
    (x0, y0), (x1, y1) = p0, p1; xm = (x0 + x1) / 2
    ax.add_patch(PathPatch(Path([(x0, y0), (xm, y0), (xm, y1), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]), fc='none', ec=color, lw=lw, capstyle='round', zorder=z))


def _schematic(ax):
    """The model in one drawing: a lobule with its six triads and the one path the optimum keeps, and the network of two trees."""
    ax.set_xlim(0, 10); ax.set_ylim(0, 10.4); ax.set_aspect('equal'); ax.axis('off')
    # one lobule
    c = np.array([5.0, 8.1]); R = 1.45; th = np.pi / 6 + np.arange(6) * np.pi / 3; corners = c + R * np.c_[np.cos(th), np.sin(th)]
    ax.add_patch(RegularPolygon(c, 6, radius=R, orientation=0, fc=BAND, ec=NEUTRAL, lw=0.8, zorder=1))
    for k, p in enumerate(corners):
        ax.plot([p[0], c[0]], [p[1], c[1]], color=SINUSOID, lw=2.0 if k == 2 else 0.6, alpha=1 if k == 2 else 0.35, solid_capstyle='round', zorder=2)
    ax.plot(corners[:, 0], corners[:, 1], 'o', ms=4, mfc=PORTAL, zorder=4, **DOT); ax.plot(*c, 'o', ms=6, mfc=HEPATIC, zorder=4, **DOT)
    ax.plot([1.9, corners[2, 0]], [corners[2, 1], corners[2, 1]], color=PORTAL, lw=2.0, solid_capstyle='round', zorder=3)
    ax.annotate('', (2.9, corners[2, 1] + 0.42), (1.9, corners[2, 1] + 0.42), arrowprops=dict(arrowstyle='-|>', lw=0.6, color=INK, mutation_scale=5))
    lead = dict(arrowstyle='-', lw=0.4, color=MUTED, shrinkA=1, shrinkB=2)
    ax.annotate('Portal triad', corners[0], (7.25, 9.75), fontsize=FS_IN, ha='left', va='center', arrowprops=lead)
    ax.annotate('Central vein', c, (7.25, 8.1), fontsize=FS_IN, ha='left', va='center', arrowprops=lead)
    ax.annotate('Sinusoidal path\nkept by the optimum', (c + corners[2]) / 2, (0.25, 7.0), fontsize=FS_IN, ha='left', va='center', arrowprops=lead)
    ax.annotate('Candidate paths', (c + corners[5]) / 2, (7.25, 6.6), fontsize=FS_IN, ha='left', va='center', arrowprops=lead)
    ax.text(0.25, 10.05, 'One lobule', fontsize=FS_IN, ha='left', va='center'); ax.text(1.9, corners[2, 1] + 0.75, 'Portal inflow', fontsize=FS_SMALL, ha='left', va='center', color=MUTED)
    # the network: portal tree, one path per lobule, hepatic venous tree
    y = 0.95 + 0.47 * np.arange(9); yc = y.mean(); xs = (0.2, 1.5, 2.9, 4.2); mirror = lambda x: 10 - x
    for tree, col, f in ((0, PORTAL, lambda x: x), (1, HEPATIC, mirror)):
        ax.plot([f(xs[0]), f(xs[1])], [yc, yc], color=col, lw=4.2, solid_capstyle='butt', zorder=2)
        for g in range(3):
            ym = y[3 * g:3 * g + 3].mean(); _curve(ax, (f(xs[1]), yc), (f(xs[2]), ym), 2.4, col)
            for yy in y[3 * g:3 * g + 3]: _curve(ax, (f(xs[2]), ym), (f(xs[3]), yy), 1.2, col)
    for yy in y:
        ax.plot([xs[3], mirror(xs[3])], [yy, yy], color=SINUSOID, lw=1.2, solid_capstyle='butt', zorder=3)
        ax.plot([xs[3], mirror(xs[3])], [yy, yy], 'o', ms=2.6, mfc=SINUSOID, mec='white', mew=0.4, zorder=4)
    for x0, x1 in ((0.25, 1.25), (8.75, 9.75)): ax.annotate('', (x1, yc + 0.55), (x0, yc + 0.55), arrowprops=dict(arrowstyle='-|>', lw=0.6, color=INK, mutation_scale=5))
    ax.text(0.25, 5.45, 'The network', fontsize=FS_IN, ha='left', va='center')
    ax.text(0.2, yc - 0.45, 'Hilum', fontsize=FS_SMALL, ha='left', va='top', color=MUTED); ax.text(9.8, yc - 0.45, 'Cava', fontsize=FS_SMALL, ha='right', va='top', color=MUTED)
    for x0, x1, col, lab in ((0.2, 4.1, PORTAL, 'Portal tree'), (4.3, 5.7, SINUSOID, 'Lobules'), (5.9, 9.8, HEPATIC, 'Hepatic venous tree')):
        ax.plot([x0, x0, x1, x1], [0.5, 0.35, 0.35, 0.5], color=col, lw=0.8, solid_capstyle='butt'); ax.text((x0 + x1) / 2, 0.12, lab, fontsize=FS_IN, ha='center', va='top')


def _network_2d(ax, R=4.5, b=0.75):
    from model import network as lc
    G = lc.build_liver_graph(R=R); gin, gout = lc.split_liver(G); xy = G['xy']; cen, cor, pairs, Lh = lc.hex_lobules(R, 1.0)
    tree = min((lc.best_sinusoid_tree(gin, b, seed=s) for s in range(3)), key=lambda s: s['D'])
    for (x, y) in cen: th = np.pi / 6 + np.arange(7) * np.pi / 3; ax.plot(x + Lh * np.cos(th), y + Lh * np.sin(th), color=GRID, lw=0.4, zorder=0)
    w, keep = tree['w'], tree['keep']; r = np.zeros(gin.m); r[keep] = (w[keep] * gin.L[keep]) ** 0.25; rm = r[keep].max()
    outl = lc.tree_allocation(gout, np.arange(gout.m), b); ro = (outl['w'] * gout.L) ** 0.25
    for k in range(gout.m): i, j = gout.E[k]; ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]], color=HEPATIC, lw=0.25 + 2.0 * ro[k] / ro.max(), alpha=0.55, zorder=1)
    for k in keep:
        i, j = gin.E[k]; sin = gin.kind[k] == 'sin'
        ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]], color=SINUSOID if sin else PORTAL, lw=0.9 if sin else 0.35 + 2.6 * r[k] / rm, zorder=3 if sin else 2)
    ax.plot(*xy[gin.src], 'o', ms=4.5, mfc=INK, zorder=5, **DOT); ax.set_aspect('equal'); ax.axis('off')
    return len(cen)


def _network_3d(ax, R=3.4, b=0.75):
    from model import network as lc, network3d as l3
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    G3 = l3.build_liver3d(R=R); gin, gout = lc.split_liver(G3); ft = l3.FastTree(gin, b)
    best = min((ft.solve(seed=s) for s in range(3)), key=lambda s: s['D']); sol = lc.tree_allocation(gin, best['keep'], b)
    r = np.zeros(gin.m); r[sol['keep']] = (sol['w'][sol['keep']] * gin.L[sol['keep']]) ** 0.25; X = G3['xy']
    o = lc.tree_allocation(gout, np.arange(gout.m), b); ro = (o['w'] * gout.L) ** 0.25; segs, cols, lws = [], [], []
    for k in range(gout.m): i, j = gout.E[k]; segs.append([X[i], X[j]]); cols.append(HEPATIC); lws.append(0.2 + 1.5 * ro[k] / ro.max())
    for k in sol['keep']:
        i, j = gin.E[k]; sin = gin.kind[k] == 'sin'; segs.append([X[i], X[j]]); cols.append(SINUSOID if sin else PORTAL); lws.append(0.6 if sin else 0.3 + 2.4 * r[k] / r.max())
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, np.pi / 2, 30)); Rh = R + 0.15
    ax.plot_surface(Rh * np.cos(u) * np.sin(v), Rh * np.sin(u) * np.sin(v), Rh * np.cos(v), color=MUTED, alpha=0.05, linewidth=0, shade=False)
    ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=lws, alpha=0.9, capstyle='round')); ax.scatter([0], [0], [0], color=INK, s=14, depthshade=False, edgecolor='white', linewidth=0.5)
    ax.set_xlim(-Rh, Rh); ax.set_ylim(-Rh, Rh); ax.set_zlim(0, Rh); ax.view_init(elev=24, azim=-58); ax.set_box_aspect((1, 1, 0.55), zoom=1.3); ax.set_axis_off()
    return G3['n_lob']


def figure_1():
    fig = fs.figure(124); a = fs.axes(fig, 2, 62, 54, 56); b = fs.axes(fig, 64, 66, 52, 52); c = fs.axes(fig, 119, 63, 60, 56, projection='3d')
    d = fs.axes(fig, 13, 11, 60, 42); e = fs.axes(fig, 112, 11, 42, 42)
    _schematic(a); n2 = _network_2d(b); n3 = _network_3d(c)
    fig.text(68 / 180, 121 / 124, f'Two dimensions, {n2} lobules', fontsize=FS_IN, ha='left', va='center'); fig.text(126 / 180, 121 / 124, f'Three dimensions, {n3} lobules', fontsize=FS_IN, ha='left', va='center')
    fig.legend(handles=[Line2D([], [], color=PORTAL, lw=1.6, label='Portal tree'), Line2D([], [], color=SINUSOID, lw=1.6, label='Sinusoidal paths'), Line2D([], [], color=HEPATIC, lw=1.6, alpha=0.75, label='Hepatic venous tree'),
                        Line2D([], [], marker='o', ls='', ms=4.5, mfc=INK, label='Hilum', **DOT)], loc='center', bbox_to_anchor=(121 / 180, 62 / 124), ncol=4, fontsize=FS_IN, columnspacing=1.4)
    grwr = np.linspace(0.4, 2.4, 200)                                  # a whole liver is 2% of body weight: g = GRWR / 2
    for h, col in zip((1.0, 1.5, 2.0), fs.ramp(3)):
        z = h / (grwr / 2); d.plot(grwr, z, color=col, zorder=3); d.text(2.44, z[-1], f'$h$ = {h:g}', fontsize=FS_IN, ha='left', va='center', clip_on=False)
    d.axhspan(1, 2, color=BAND, lw=0, zorder=0); d.text(0.46, 1.5, 'Window', fontsize=FS_IN, color=MUTED, ha='left', va='center')
    d.set_xlim(0.4, 2.4); d.set_ylim(0, 7); d.set_xticks([0.5, 1.0, 1.5, 2.0]); d.set_yticks(range(0, 8))
    d.set_xlabel('Graft size, GRWR (%)'); d.set_ylabel('Load per lobule, $z$ = $h/g$')
    plane(e); e.text(1.14, 0.7, 'Donor', fontsize=FS_IN)
    e.text(0.25, 4.3, f'{ZR} > 1\nPressure above flow', fontsize=FS_IN, ha='left', va='center'); e.text(4.85, 1.2, f'{ZR} < 1\nPressure below flow', fontsize=FS_IN, ha='right', va='center')
    fs.letter(fig, a, 'A', 1.0, 3); fs.letter(fig, b, 'B', 1.0, 3); fs.letter(fig, c, 'C', -2.0, 2); fs.letter(fig, d, 'D', 10, 4); fs.letter(fig, e, 'E', 10, 4); save(fig, 'Figure_1')


# ------------------------------------------------------------------ Figure 2
def label_of(r):
    return r.study + (', ' + ('splenectomy' if r.group == 'splenectomy' else 'no splenectomy') if r.study == 'Wang 2014' else '')


# the mechanistic class of a group is a symbol, never a colour: (marker, filled, label)
CONFIGURATION = {'outflow reconstruction': ('o', True, 'Outflow reconstruction'), 'combined': ('D', True, 'Outflow and pressure control'),
                 'inflow modulation': ('^', True, 'Inflow modulation (splenectomy)'), 'none': ('o', False, 'No manoeuvre')}


def figure_2():
    d = series(); p = d.dropna(subset=['zF', 'zP']); flow_only = d[d.zF.notna() & d.zP.isna() & d.in_main_analysis]; P = load('predictions')['P2_discordance']
    U = {(g['study'], g['group']): g for g in P['uncertainty']['groups']}; med = float(np.median(p.zR))
    fig = fs.figure(76); a = fs.axes(fig, 12, 12, 56, 56); b = fs.axes(fig, 126, 12, 51, 56)
    plane(a); a.plot([0, 5], [0, 5 * med], color=MUTED, lw=LW_REF, zorder=2); a.text(4.62, 4.62 * med + 0.2, f'Median {med:.2f}', fontsize=FS_SMALL, ha='center', va='center', color=MUTED, rotation=np.degrees(np.arctan(med)))
    a.plot(flow_only.zF, np.zeros(len(flow_only)), '|', ms=7, mew=0.9, color=MUTED, clip_on=False, zorder=4)
    order = {(r.study, r.group): k + 1 for k, r in enumerate(sorted(p.itertuples(), key=lambda r: r.zR))}; left = {'Yagi 2006', 'Wu 2011', 'Ishizaki 2012'}
    for _, r in p.iterrows():
        mk, filled, _ = CONFIGURATION[r.configuration]; ms = (3.2 + r.n ** 0.5 / 3.2) * (0.85 if mk == 'D' else 1); dx = (0.09 + ms / 95) * (-1 if r.study in left else 1)
        a.plot(r.zF, r.zP, mk, ms=ms, mfc=BLUE if filled else 'white', mec=BLUE, mew=1.0, zorder=6)
        a.text(r.zF + dx, r.zP, str(order[(r.study, r.group)]), fontsize=FS_IN, ha='right' if r.study in left else 'left', va='center', zorder=7)
    a.text(0.72, 1.28, 'Donor', fontsize=FS_IN, ha='right'); a.text(2.5, 0.22, 'Flow load only', fontsize=FS_IN, color=MUTED, ha='center', va='bottom')
    a.legend(handles=[Line2D([], [], marker=mk, ls='', ms=4.2 if mk != 'D' else 3.8, mfc=BLUE if filled else 'white', mec=BLUE, mew=1.0, label=lab) for mk, filled, lab in CONFIGURATION.values()],
             loc='upper left', handletextpad=0.2, borderaxespad=0.1)
    rows = sorted(p.itertuples(), key=lambda r: r.zR); top = len(rows) - 0.35
    b.fill_between(MODEL_ZR, -0.6, top, zorder=0, **fs.MODEL_TINT); b.plot([1, 1], [-0.6, top], color=INK, lw=LW_REF, ls=DASH, zorder=1)
    for k, r in enumerate(rows):
        y = len(rows) - 1 - k; g = U[(r.study, r.group)]; mk, filled, _ = CONFIGURATION[r.configuration]
        b.plot(g['ci_total'], [y, y], color=BLUE, lw=0.6, solid_capstyle='butt', zorder=2); b.plot(g['ci_sampling'], [y, y], color=BLUE, lw=2.2, solid_capstyle='butt', zorder=2); b.plot(r.zR, y, mk, ms=5 if mk != 'D' else 4.4, mfc=BLUE if filled else 'white', mec=BLUE, mew=1.0, zorder=3)
        b.text(-0.04, y + 0.13, f'{k + 1}  ' + label_of(r), fontsize=FS_IN, ha='right', va='center'); b.text(r.zR, y + 0.2, f'{r.zR:.2f}', fontsize=FS_IN, ha='center', va='bottom')
        b.text(-0.04, y - 0.2, f'$n$ = {int(r.n)}, {ZF} = {r.zF:.2f}', fontsize=FS_SMALL, color=MUTED, ha='right', va='center')
    for yk, lw, lab in ((0.62, 2.2, 'Sampling'), (0.27, 0.6, 'With references')):          # key in the empty corner, away from the band of the model
        b.plot([0.03, 0.09], [yk, yk], color=BLUE, lw=lw, solid_capstyle='butt', zorder=2); b.text(0.11, yk, lab, fontsize=FS_SMALL, color=MUTED, ha='left', va='center')
    b.text(np.mean(MODEL_ZR), top + 0.07, 'Network\nmodel', fontsize=FS_IN, ha='center', va='bottom', linespacing=0.95)
    b.text(-0.04, top + 0.3, f"{P['cohorts_below_one']} of {P['cohorts']} cohorts below 1\nExact sign test, one-sided $P$ = {P['sign_test_cohorts_one_sided_p']:.3f}", fontsize=FS_IN, ha='right', va='center')
    b.set_xlim(0, 1.2); b.set_ylim(-0.6, len(rows) + 0.45); fs.bare(b); b.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]); b.set_xlabel(f'Resistance load, {ZR} = {ZP} / {ZF}')
    fs.letter(fig, a, 'A', 10, 4); fs.letter(fig, b, 'B', 44, 4); save(fig, 'Figure_2')


# ------------------------------------------------------------------ Figure 3
def figure_3():
    E = load('bayes_main')['exploratory (all outcomes)']; S = load('spec_curve'); d = series(); name = lambda s: s.split(' / ')[0]
    strata = sorted(E['strata'], key=lambda s: E['beta_by_stratum'][s]['mean']); centre = {r.study: r.centre_id for r in d.itertuples()}
    hue = dict(zip(sorted({centre[name(s)] for s in strata}), fs.CATEGORICAL)); shapes = {}                # colour is the centre, the symbol the publication
    for s in sorted(strata): shapes[s] = 'os^Dv'[sum(centre[name(t)] == centre[name(s)] for t in shapes)]
    col = {s: hue[centre[name(s)]] for s in strata}
    fig = fs.figure(68); a = fs.axes(fig, 19, 12, 31, 49); b = fs.axes(fig, 66, 12, 48, 49); c = fs.axes(fig, 152, 12, 25, 49)
    top = len(strata) - 0.45; a.fill_betweenx([-0.55, top], *E['mu_beta_ci'], color=BAND, lw=0, zorder=0); a.plot([E['mu_beta']] * 2, [-0.55, top], color=INK, lw=0.9, zorder=1, solid_capstyle='butt'); a.plot([0, 0], [-1.4, top], color=MUTED, lw=0.5, ls=DASH)
    for k, s in enumerate(strata):
        v = E['beta_by_stratum'][s]; a.plot(v['ci'], [k, k], color=col[s], zorder=2); a.plot(v['mean'], k, shapes[s], ms=4.8, mfc=col[s], zorder=3, **DOT)
        a.text(-0.62, k + 0.14, name(s), fontsize=FS_IN, ha='right', va='center'); a.text(-0.62, k - 0.2, centre[name(s)], fontsize=FS_SMALL, color=MUTED, ha='right', va='center')
    pp = E['prob_mu_positive']; new = E['prediction_interval']; a.plot(new, [-0.95, -0.95], color=INK, lw=0.9, zorder=2, solid_capstyle='butt'); a.plot(new, [-0.95, -0.95], '|', ms=4, mew=0.9, color=INK, ls='')
    a.text(-0.62, -0.95, 'New centre', fontsize=FS_IN, ha='right', va='center')
    a.text(E['mu_beta'], len(strata) - 0.3, f"$\\mu_\\beta$ = {E['mu_beta']:+.2f}\n$P$($\\mu_\\beta$ > 0) " + ('> 0.999' if pp > 0.999 else f'= {pp:.3f}'), fontsize=FS_IN, ha='center', va='bottom', linespacing=1.05)
    a.set_xlim(-0.5, 3.5); a.set_ylim(-1.4, len(strata) + 0.75); fs.bare(a); a.set_xticks([0, 1, 2, 3]); a.set_xlabel(f'Slope per unit of {ZP}')
    m = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]; m = m.assign(stratum=m.study + ' / ' + m.outcome_type); m = m[m.stratum.isin(E['strata'])]
    b.axvspan(1, 2, color=BAND, lw=0, zorder=0)
    for s in E['strata']:
        g = m[m.stratum == s]; z = np.linspace(max(0.8, g.zP.min() - 0.35), min(4.0, g.zP.max() + 0.35), 60)
        b.plot(z, 100 / (1 + np.exp(-(E['alphas'][s] + E['beta_by_stratum'][s]['mean'] * z))), color=col[s], zorder=2)
        b.scatter(g.zP, g.pct, s=10 + g.n / 5, marker=shapes[s], color=col[s], edgecolor='white', linewidth=0.6, zorder=3)
    b.set_xlim(0.8, 4.0); b.set_ylim(0, 60); b.set_xticks([1, 2, 3, 4]); b.set_xlabel(f'Pressure load, {ZP}'); b.set_ylabel('Adverse outcome (%)')
    S = sorted(S, key=lambda s: s['mu_beta']); c.axvline(0, color=MUTED, lw=0.5, ls=DASH)
    short = {'primary outcome (SFSS or early dysfunction)': 'Primary outcome', 'secondary outcome (mortality or graft loss)': 'Secondary outcome', 'all outcomes': 'All outcomes',
             'gradients reported or derived only': 'Gradients reported or derived', 'no imputed CVP': 'No imputed CVP',
             'events reported only (no counts from percentages)': 'Events reported only', 'one Kyoto publication (Ogura 2010)': 'Kyoto: Ogura 2010 only', 'one Kyoto publication (Uemura 2016)': 'Kyoto: Uemura 2016 only',
             'one Kyoto publication (Yagi 2005)': 'Kyoto: Yagi 2005 only', 'tau prior HalfNormal(0.25)': 'Prior on $\\tau_\\beta$: HN(0.25)', 'tau prior HalfNormal(1.0)': 'Prior on $\\tau_\\beta$: HN(1.0)',
             'one stratum per centre, not per publication': 'One stratum per centre', 'assumed CVP 3 mmHg': 'CVP 3 mmHg', 'assumed CVP 7 mmHg': 'CVP 7 mmHg', 'assumed CVP 9 mmHg': 'CVP 9 mmHg'}
    for k, s in enumerate(S):
        c.plot(s['mu_beta_ci'], [k, k], color=BLUE, lw=1.0); c.plot(s['mu_beta'], k, 'o', ms=3.8, mfc=BLUE, **DOT)
        lab = short.get(s['spec'], s['spec']); c.text(-0.72, k, lab[0].upper() + lab[1:], fontsize=FS_SMALL, ha='right', va='center')
    c.set_xlim(-0.6, 3.6); c.set_ylim(-0.7, len(S) - 0.3); fs.bare(c); c.set_xticks([0, 1, 2, 3]); c.set_xlabel('$\\mu_\\beta$ across specifications')
    fs.letter(fig, a, 'A', 17, 4); fs.letter(fig, b, 'B', 10, 4); fs.letter(fig, c, 'C', 36, 4); save(fig, 'Figure_3')


# ------------------------------------------------------------------ Figure 4
def figure_4():
    P = load('predictions')['P4_thresholds']; ref = P['reference_sensitivity']
    rows = [('Clinically significant portal hypertension', 'HVPG 10 mmHg', 'cirrhosis_HVPG_10'), ('Partial graft', 'PVP 15 mmHg at CVP 5 mmHg', 'graft_PVP15_CVP5'),
            ('Resection limit in cirrhosis', 'HVPG 10 mmHg', 'resection_HVPG_10'), ('Variceal bleeding', 'HVPG 12 mmHg', 'variceal_HVPG_12')]
    fig = fs.figure(64); a = fs.axes(fig, 51, 12, 34, 45); b = fs.axes(fig, 108, 12, 69, 45)
    a.axvspan(1, 2, color=BAND, lw=0, zorder=0)
    for k, (lab, sub, key) in enumerate(rows):
        y = len(rows) - 1 - k; lo, mid, hi = ref['6 mmHg reference'][key], ref['5 mmHg reference'][key], ref['4 mmHg reference'][key]
        a.plot([lo, hi], [y, y], color=BLUE, zorder=2); a.plot(mid, y, 'o', ms=5, mfc=BLUE, zorder=3, **DOT)
        a.text(mid, y + 0.2, f'{mid:.1f}', fontsize=FS_IN, ha='center', va='bottom')
        a.text(-0.12, y + 0.12, lab, fontsize=FS_IN, ha='right', va='center'); a.text(-0.12, y - 0.2, sub, fontsize=FS_SMALL, color=MUTED, ha='right', va='center')
    a.set_xlim(0, 3.5); a.set_ylim(-0.6, len(rows) - 0.3); fs.bare(a); a.set_xticks([0, 1, 2, 3]); a.set_xlabel(f'Pressure load, {ZP}')
    removed = np.linspace(0, 75, 200); g = 1 - removed / 100; b.axhspan(1, 2, color=BAND, lw=0, zorder=0); b.axvline(60, color=MUTED, lw=0.5, ls=DASH)
    for hv, colr in zip((5, 8, 10, 12), fs.ramp(4)):
        z = hv / indices.NORMAL_GRADIENT / g; b.plot(removed, z, color=colr, zorder=3); out = z.max() > 4.75; k = np.argmin(abs(z - 4.75)) if out else len(z) - 1
        b.text(removed[k] - 1.0 if out else 75.8, 4.75 if out else z[k], f'{hv}', fontsize=FS_IN, ha='right' if out else 'left', va='bottom' if out else 'center')
    b.text(3, 4.7, 'Baseline HVPG (mmHg)', fontsize=FS_IN, ha='left', va='center'); b.text(59, 0.35, 'Right hepatectomy', fontsize=FS_IN, color=MUTED, ha='right', va='center')
    b.set_xlim(0, 80); b.set_ylim(0, 5); b.set_xticks([0, 20, 40, 60, 80]); b.set_xlabel('Parenchyma removed (%)'); b.set_ylabel(f'Pressure load of the remnant, {ZP}')
    fs.letter(fig, a, 'A', 49, 4); fs.letter(fig, b, 'B', 10, 4); save(fig, 'Figure_4')


# ------------------------------------------------------------------ Figure S1
def closure():
    """Exponents of tau0 and of the vascular mass under an invariant portal pressure (predictions.p1_allometry)."""
    from model import scaling as sc
    gs = np.arange(8, 19); rows = []
    for f, flow in ((0.77, 'portal'), (0.88, 'total')):
        for l, n in ((0.0, 0.89), (0.10, 0.59)):
            for b in (1.0, 0.9, 0.85, 0.8, 0.75, 2 / 3, 0.5):
                R = [sc.canopy_optimum(g, b, 3) for g in gs]; N = np.array([r['N'] for r in R])
                eD = np.polyfit(np.log(N), np.log([r['D'] for r in R]), 1)[0]; eV = np.polyfit(np.log(N), np.log([r['V'] for r in R]), 1)[0]; c = b * (f + eD * n + 3 * l) / 2
                rows.append(dict(flow=flow, l=l, b=b, tau0=(2 * f - (2 + b) * c / b + eD * n + 3 * l) / 2, V=c / b + eV * n))
    return pd.DataFrame(rows)


def figure_s1():
    from model import scaling as sc
    g2 = pd.read_csv(f'{DATA}/graph_scaling_all.tsv', sep='\t'); g3 = pd.read_csv(f'{DATA}/graph3d_all.tsv', sep='\t'); cl = closure()
    fig = fs.figure(136); a1 = fs.axes(fig, 14, 84, 70, 44); a2 = fs.axes(fig, 105, 84, 70, 44); b = fs.axes(fig, 14, 12, 70, 44); c = fs.axes(fig, 105, 12, 70, 44)
    def trend(ax, N, D, col, marker, lab, last=None):
        k = slice(-last, None) if last else slice(None); e, i = np.polyfit(np.log(N[k]), np.log(D[k]), 1); x = np.array([N.min(), N.max()])
        ax.plot(x, np.exp(i) * x ** e, color=col, lw=1.0, zorder=2); ax.plot(N, D, marker, ms=4.6, mfc=col, ls='', zorder=3, label=lab, **DOT)
    bs = sorted(g2.b.unique())
    for bb, col in zip(bs, fs.ramp(len(bs))):
        gr = g2[g2.b == bb]; Dt = (gr.D_in_perF2 + gr.D_out / gr.N ** 2).values; e = np.polyfit(np.log(gr.N), np.log(Dt), 1)[0]; trend(a1, gr.N.values, Dt / Dt[0], col, 'o', f'$b$ = {bb:g}, exponent {e:.2f}')
    bs = sorted(g3.b.unique())
    for bb, col in zip(bs, fs.ramp(len(bs))):
        gr = g3[g3.b == bb]; D = gr.D_tot_perF2.values; eg = np.polyfit(np.log(gr.N.values[-3:]), np.log(D[-3:]), 1)[0]
        ana = [sc.canopy_optimum(g, bb, 3) for g in range(4, 9)]; ec = np.polyfit(np.log([q['N'] for q in ana][-3:]), np.log([q['D'] for q in ana][-3:]), 1)[0]
        trend(a2, gr.N.values, D / D[0], col, 's', f'$b$ = {bb:.3g}, graph {eg:.2f}, closed form {ec:.2f}', last=3)
    for ax, lab in ((a1, 'Two-dimensional disc'), (a2, 'Three-dimensional hemisphere')):
        ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('Number of lobules, $N$'); ax.set_ylabel('Minimum dissipation, $D^{*}/F^{2}$ (normalised)')
        ax.legend(loc='upper left', title=lab, title_fontsize=FS_IN, alignment='left', handletextpad=0.2)
    a1.set_xticks([40, 100, 300]); a1.set_xticklabels(['40', '100', '300']); a1.minorticks_off(); a2.set_xticks([100, 300, 1000]); a2.set_xticklabels(['100', '300', '1,000']); a2.xaxis.set_minor_formatter(plt.NullFormatter())
    hue = {'portal': PORTAL, 'total': SINUSOID}                       # the sinusoids carry the total hepatic flow
    for (flow, l), gr in cl.groupby(['flow', 'l']):
        gr = gr.sort_values('b'); kw = dict(color=hue[flow], ls='-' if l > 0 else (0, (1, 1.6)), marker='o' if l > 0 else 's', ms=3.6, mfc=hue[flow], mec='white', mew=0.5)
        b.plot(gr.b, gr.tau0, **kw); c.plot(gr.b, gr.V, **kw)
    b.axhspan(-0.05, 0.05, color=BAND, lw=0, zorder=0); b.axhline(0, color=INK, lw=0.5); b.text(0.985, 0.058, '±0.05', fontsize=FS_SMALL, color=MUTED, ha='right', va='bottom')
    c.axhspan(0.80, 0.92, color=BAND, lw=0, zorder=0); c.text(0.515, 0.812, 'Observed, hepatic blood volume: 0.86 (0.80 to 0.92)', fontsize=FS_SMALL, color=MUTED, ha='left', va='bottom')
    inside = int(((cl.V >= 0.80) & (cl.V <= 0.92)).sum()); c.text(0.985, 1.088, f'{inside} of {len(cl)} combinations inside', fontsize=FS_SMALL, color=MUTED, ha='right', va='top')
    for ax, lab in ((b, 'Exponent of $\\tau_0$ with body mass'), (c, 'Vascular-mass exponent with body mass')):
        ax.set_xlim(0.47, 1.03); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0]); ax.set_xlabel('Maintenance exponent, $b$'); ax.set_ylabel(lab)
    b.set_ylim(-0.09, 0.09); b.set_yticks([-0.08, -0.04, 0, 0.04, 0.08]); c.set_ylim(0.78, 1.10)
    H = [Line2D([], [], color=hue[f], ls='-' if l else (0, (1, 1.6)), marker='o' if l else 's', ms=3.6, mec='white', mew=0.5, label=f"{'Portal' if f == 'portal' else 'Total hepatic'} flow ~ $M^{{{0.77 if f == 'portal' else 0.88}}}$, lobule {'~ $M^{0.10}$' if l else 'invariant'}")
         for f in ('portal', 'total') for l in (0.10, 0.0)]
    fig.legend(handles=H, loc='center', bbox_to_anchor=(0.5, 67 / 136), ncol=2, columnspacing=2.5, handlelength=2.6)
    fs.letter(fig, a1, 'A', 12, 4); fs.letter(fig, b, 'B', 12, 4); fs.letter(fig, c, 'C', 12, 4); save(fig, 'Figure_S1')


# ------------------------------------------------------------------ Figure S2
def figure_s2(h=1.5):
    from model import whole_graft as wg
    dd = np.linspace(1.0, 0.5, 26); out, inn = wg.sweep('out', d=dd, h=h), wg.sweep('in', d=dd, h=h)
    fig = fs.figure(62); a = fs.axes(fig, 2, 6, 66, 52); b = fs.axes(fig, 82, 12, 44, 44); c = fs.axes(fig, 134, 12, 44, 44)
    a.set_xlim(0, 12); a.set_ylim(0, 9.4); a.axis('off'); tint = lambda col: tuple(0.82 + 0.18 * v for v in plt.matplotlib.colors.to_rgb(col))
    boxes = [(0.25, 'Recipient\nportal trunk', BAND, MUTED), (2.85, 'Inlet\n$R_{in}$ ~ $d_{in}^{\\,-4}$', tint(PORTAL), PORTAL), (5.45, 'Lobules\n$R_{lob}$', tint(SINUSOID), SINUSOID), (8.05, 'Outlet\n$R_{out}$ ~ $d_{out}^{\\,-4}$', tint(HEPATIC), HEPATIC)]
    for x, lab, fc, ec in boxes:
        a.add_patch(FancyBboxPatch((x, 3.4), 2.15, 1.9, boxstyle='round,pad=0.03,rounding_size=0.18', fc=fc, ec=ec, lw=0.8)); a.text(x + 1.075, 4.35, lab, ha='center', va='center', fontsize=FS_SMALL, linespacing=1.25)
    arrow = dict(arrowstyle='-|>', lw=0.6, color=INK, mutation_scale=5, shrinkA=0, shrinkB=0)
    for x in (2.44, 5.04, 7.64): a.annotate('', (x + 0.37, 4.35), (x, 4.35), arrowprops=arrow)
    a.annotate('', (10.9, 4.35), (10.24, 4.35), arrowprops=arrow); a.text(11.0, 4.35, 'Cava', fontsize=FS_SMALL, ha='left', va='center')
    a.plot([1.32, 1.32, 10.55, 10.55], [3.36, 1.9, 1.9, 4.2], color=MUTED, lw=0.7, ls=DASH); a.text(5.9, 1.55, 'Collateral path $R_{coll}$, in parallel', ha='center', va='top', fontsize=FS_SMALL, color=MUTED)
    for x, y0, lab in ((1.32, 5.34, 'Portal pressure,\ntrunk'), (5.22, 4.35, 'Portal pressure,\ngraft')):
        a.plot([x, x], [y0, 6.5], color=MUTED, lw=0.5); a.plot([x], [y0], 'o', ms=2.2, color=MUTED); a.text(x, 6.65, lab, fontsize=FS_SMALL, ha='center', va='bottom', color=MUTED)
    x = [r['diameter_ratio'] for r in out]
    b.plot(x, [r['zP'] for r in out], color=HEPATIC, label=f'Outlet narrowed, {ZP}'); b.plot(x, [r['zF'] for r in out], color=HEPATIC, ls=DASH, zorder=4, label=f'Outlet narrowed, {ZF}')
    b.plot(x, [r['zP'] for r in inn], color=PORTAL, label=f'Inlet narrowed, {ZP} = {ZF}'); b.plot(x, [r['portal_trunk_pressure'] for r in inn], color=PORTAL, lw=1.0, ls=(0, (1, 1.6)), label='Inlet narrowed, trunk pressure')
    b.invert_xaxis(); b.set_xticks([1.0, 0.9, 0.8, 0.7, 0.6, 0.5]); b.set_xlabel('Anastomosis diameter (fraction of normal)'); b.set_ylabel('Load (donor = 1)'); b.set_ylim(0, 4); b.set_yticks([0, 1, 2, 3, 4]); b.legend(loc='upper left', handlelength=2.3)
    plane(c, lim=3.0, model=False, label_at=2.45)
    for rows, col in ((out, HEPATIC), (inn, PORTAL)):
        xf, yp = [r['zF'] for r in rows], [r['zP'] for r in rows]; c.plot(xf, yp, color=col, lw=1.5, zorder=4); c.plot(xf[-1], yp[-1], 'o', ms=4.5, mfc=col, zorder=5, **DOT)
    lead = dict(arrowstyle='-', color=MUTED, lw=0.4, shrinkA=1, shrinkB=2); c.plot(h, h, 'o', ms=4.5, mfc=INK, zorder=6, **DOT)
    c.annotate(f'Normal anastomoses,\ninflow {h:g} × donor', (h, h), (1.5, 0.55), fontsize=FS_SMALL, ha='left', va='center', arrowprops=lead)
    c.annotate('Outlet at half\nits diameter', (out[-1]['zF'], out[-1]['zP']), (0.12, 2.6), fontsize=FS_SMALL, ha='left', va='center', arrowprops=lead)
    c.annotate('Inlet at half\nits diameter', (inn[-1]['zF'], inn[-1]['zP']), (0.12, 1.85), fontsize=FS_SMALL, ha='left', va='center', arrowprops=lead)
    c.text(0.92, 1.12, 'Donor', fontsize=FS_SMALL, ha='right', va='bottom')
    fs.letter(fig, a, 'A', 0, 0.5); fs.letter(fig, b, 'B', 10, 4); fs.letter(fig, c, 'C', 10, 4); save(fig, 'Figure_S2')


# ------------------------------------------------------------------ Figures S3 and S4
def figure_s3():
    E = load('ipd_ecological'); sw = E['sweep']; x = np.array([s['sigma_mmHg'] for s in sw]); y = np.array([s['mu_beta'] for s in sw]); ci = np.array([s['mu_beta_ci'] for s in sw])
    fig = fs.figure(64, fs.SINGLE); ax = fs.axes(fig, 12, 12, 73, 48)
    for sd in E['reported_within_group_sd_mmHg']: ax.axvline(sd, color=NEUTRAL, lw=0.7, zorder=0)
    ax.fill_between(x, ci[:, 0], ci[:, 1], zorder=1, **fs.MODEL_TINT); ax.axhline(0, color=INK, lw=0.5); ax.axhline(y[0], color=INK, lw=LW_REF, ls=DASH, zorder=2)
    ax.plot(x, y, color=BLUE, zorder=3); ax.plot(x, y, 'o', ms=4.5, mfc=BLUE, zorder=4, **DOT)
    ax.text(0.08, y[0] - 0.09, f'Group-mean estimate, {y[0]:+.2f}', fontsize=FS_SMALL, ha='left', va='top'); ax.text(np.mean(E['reported_within_group_sd_mmHg']), 3.82, 'Dispersions reported', fontsize=FS_SMALL, color=MUTED, ha='center', va='top')
    pp = min(s.get('prob_positive', s.get('prob_mu_positive', 1.0)) for s in sw); ax.text(0.08, 0.5, f'$P$($\\mu_\\beta$ > 0) > {np.floor(pp * 100) / 100:.2f} at every dispersion', fontsize=FS_SMALL, ha='left', va='bottom')
    ax.text(0.08, ci[0, 1] - 0.08, '95% credible interval', fontsize=FS_SMALL, color=MUTED, ha='left', va='top')
    ax.set_xlim(-0.2, 6.2); ax.set_ylim(-0.15, 3.9); ax.set_xticks(range(0, 7)); ax.set_yticks([0, 1, 2, 3]); ax.set_xlabel('Within-group SD of the gradient (mmHg)'); ax.set_ylabel('Slope for a patient, $\\mu_\\beta$')
    save(fig, 'Figure_S3')


def figure_s4():
    from ipd import patient_level as pl
    d = pl.load().dropna(subset=['gradient_before', 'gradient_after']); bad = (d.sfss.astype(bool) | d.graft_loss_or_death.astype(bool)).values
    fig = fs.figure(68, fs.SINGLE); ax = fs.axes(fig, 12, 9, 60, 55); ax.plot([-0.55, 1.03], [10, 10], color=INK, lw=LW_REF, ls=DASH, zorder=1)
    for (_, r), ev in zip(d.iterrows(), bad):
        col = EVENT if ev else '#9a9993'; ax.plot([0, 1], [r.gradient_before, r.gradient_after], color=col, lw=1.5 if ev else 0.9, zorder=4 if ev else 2); ax.plot([0, 1], [r.gradient_before, r.gradient_after], 'o', ms=4.5 if ev else 3.8, mfc=col, ls='', zorder=5 if ev else 3, **DOT)
    e = d[bad].iloc[0]; ax.text(1.06, e.gradient_after, 'Small-for-size\nsyndrome,\ngraft loss', fontsize=FS_SMALL, ha='left', va='center'); ax.text(1.06, 10, f'{ZP} = 2', fontsize=FS_SMALL, ha='left', va='center')
    ax.text(-0.5, 3.0, f'{len(d)} recipients; one event,\nin the recipient with\nthe highest gradient\n(exact $P$ = {1 / len(d):.2f})', fontsize=FS_SMALL, ha='left', va='center', color=MUTED)
    med = (d.gradient_before.median(), d.gradient_after.median()); ax.plot([-0.1, 1.1], med, '_', ms=9, mew=1.2, color=INK, ls='', zorder=6); ax.text(-0.17, med[0], f'Median {med[0]:g}', fontsize=FS_SMALL, ha='right', va='center')
    ax.text(1.06, med[1] - 0.9, f'Median {med[1]:g}', fontsize=FS_SMALL, ha='left', va='center')
    ax.set_xlim(-0.55, 1.3); ax.set_ylim(0, 25); ax.set_xticks([0, 1]); ax.set_xticklabels(['Shunt clamped', 'Shunt open']); ax.set_yticks([0, 5, 10, 15, 20, 25]); ax.set_ylabel('Portocaval gradient (mmHg)'); ax.spines['bottom'].set_bounds(0, 1)
    save(fig, 'Figure_S4')


# ------------------------------------------------------------------ Figures S5 and S6
def figure_s5():
    from model import human
    fig = fs.figure(62); a = fs.axes(fig, 11, 11, 41, 44); b = fs.axes(fig, 66, 11, 32, 44); c = fs.axes(fig, 134, 11, 43, 44)
    g, r, l, N = human.cast_arrays('portal vein'); x = np.logspace(0, np.log10(6e4), 80)
    for bb, col in zip((0.5, 0.75, 1.0), fs.ramp(3)): a.plot(x, human.radius_at(human.liver(bb, liver_g=human.CAST_LIVER_G), x), color=col, zorder=2, label=f'Model, $b$ = {bb:g}')
    a.plot(N, r, 'o', ms=4.4, mfc=INK, zorder=3, label='Human cast', **DOT); a.legend(loc='lower left', handlelength=1.4)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(0.6, 8e4); a.set_ylim(0.1, 25); a.set_xticks([1, 10, 100, 1e3, 1e4]); a.set_xticklabels(['1', '10', '10$^2$', '10$^3$', '10$^4$']); a.set_yticks([0.1, 1, 10]); a.set_yticklabels(['0.1', '1', '10'])
    a.set_xlabel('Vessels in the generation'); a.set_ylabel('Portal vein radius (mm)')
    for name, col, mk, lab, xy, ha in (('portal vein', PORTAL, 'o', 'Portal vein', (9.6, 0.26), 'center'), ('hepatic vein', HEPATIC, 's', 'Hepatic veins', (5.9, 1.10), 'right')):
        p = human.cast_poiseuille(name); b.plot(p['generations'], p['drop_mmHg'], color=col, zorder=2); b.plot(p['generations'], p['drop_mmHg'], mk, ms=4.4, mfc=col, ls='', zorder=3, **DOT); b.text(xy[0], xy[1], lab, fontsize=FS_IN, ha=ha, va='bottom')
    b.set_xlim(0.4, 11.6); b.set_ylim(0, 1.5); b.set_xticks([1, 3, 5, 7, 9, 11]); b.set_yticks([0, 0.5, 1.0, 1.5]); b.set_xlabel('Generation'); b.set_ylabel('Pressure drop per generation (mmHg)')
    c1, c0 = human.liver(0.75), human.liver(0.75, lobule_drop=None); pv = 0.5 * sum(human.cast_tree_drop('portal vein', G) for G in human.CAST_GENERATIONS['portal vein']) / 5 * 100
    sh = lambda q: [(q['share_inlet'] * 100, PORTAL), (q['share_sinusoid'] * 100, SINUSOID), (q['share_outlet'] * 100, HEPATIC)]
    rows = [('Model, lobule imposed', sh(c1)), ('Model, lobule at optimum', sh(c0)), ('Rat (Shibayama 1975)', [(42 / 90 * 100, PORTAL), (40 / 90 * 100, SINUSOID), (8 / 90 * 100, HEPATIC)]),
            ('Rat (Bohlen 1991)', [(56, NEUTRAL), (44, HEPATIC)]), ('Dog (Bohlen 1991)', [(69, NEUTRAL), (31, HEPATIC)]), ('Human cast (Poiseuille)', [(pv, PORTAL), (100 - pv, NEUTRAL)])]
    for k, (lab, segs) in enumerate(rows):
        y = len(rows) - 1 - k; left = 0
        for v, col in segs:
            c.barh(y, v, left=left, height=0.62, color=col, edgecolor='white', linewidth=1.0, zorder=2)
            if v >= 9 and col != NEUTRAL: c.text(left + v / 2, y, f'{v:.0f}', ha='center', va='center', fontsize=FS_IN, color='white', zorder=3)
            left += v
        c.text(-2.5, y, lab, ha='right', va='center', fontsize=FS_IN)
    c.set_xlim(0, 100); c.set_ylim(-0.6, len(rows) - 0.4); fs.bare(c); c.set_xticks([0, 25, 50, 75, 100]); c.set_xlabel('Share of the portocaval pressure drop (%)')
    fig.legend(handles=[Patch(color=PORTAL, label='Portal tree'), Patch(color=SINUSOID, label='Sinusoids'), Patch(color=HEPATIC, label='Hepatic veins'), Patch(color=NEUTRAL, label='Not resolved')],
               loc='lower right', bbox_to_anchor=(177.5 / 180, 56.2 / 62), ncol=4, handlelength=0.9, handleheight=0.8, columnspacing=0.9, handletextpad=0.4, borderaxespad=0)
    fs.letter(fig, a, 'A', 10, 4); fs.letter(fig, b, 'B', 10, 4); fs.letter(fig, c, 'C', 32, 4); save(fig, 'Figure_S5')


def figure_s6():
    from model import human
    fig = fs.figure(118); a = fs.axes(fig, 13, 69, 66, 42); b = fs.axes(fig, 108, 69, 66, 42); c = fs.axes(fig, 13, 11, 66, 42); d = fs.axes(fig, 130, 11, 44, 42)
    sizes = np.geomspace(human.NEWBORN_LIVER_G, human.ADULT_LIVER_G, 25); V = sizes / human.DENSITY
    for bb, col in zip((0.6, 0.75, 0.93), fs.ramp(3)): a.plot(V, [r['portal_diameter_mm'] for r in human.growth(bb, 1.0, sizes)['rows']], color=col, label=f'Model, $b$ = {bb:g}', zorder=2)
    for k, (v, (w, dia, vel)) in enumerate(zip(human.measured_growth()['liver_volume_mL'], human.PORTAL_BY_WEIGHT)):
        a.plot([v, v], dia, color=INK, lw=1.6, zorder=3, label='Children, range by body weight' if k == 0 else None)
    a.plot([human.urata_volume(70)], [human.PORTAL_DIAMETER_MM['adults'][0]], 's', ms=4.2, mfc=INK, zorder=4, label='Adults, pooled ultrasound', **DOT)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(120, 1700); a.set_ylim(2.5, 60); a.set_yticks([3, 5, 10, 20]); a.set_yticklabels(['3', '5', '10', '20']); a.set_ylabel('Portal vein diameter (mm)'); a.legend(loc='upper left', handlelength=1.4, borderaxespad=0.1)
    for s, ls, lab, xy in ((0.0, '-', 'Lobules added ($s$ = 0)', (150, 1.04)), (0.5, DASH, '$s$ = 0.5', (170, 0.62)), (1.0, (0, (1, 1.6)), 'Lobules enlarged ($s$ = 1)', (330, 0.14))):
        b.plot(V, [r['lobule_drop_mmHg'] for r in human.growth(0.75, s, sizes)['rows']], color=INK, ls=ls, zorder=2); b.text(xy[0], xy[1], lab, fontsize=FS_IN, ha='left', va='bottom')
    b.set_xscale('log'); b.set_xlim(120, 1700); b.set_ylim(0, 1.25); b.set_yticks([0, 0.5, 1.0]); b.set_ylabel('Pressure drop across the lobule (mmHg)')
    for ax in (a, b): ax.set_xticks([150, 300, 600, 1200]); ax.set_xticklabels(['150', '300', '600', '1,200']); ax.minorticks_off(); ax.set_xlabel('Liver volume (mL)')
    names = list(dict.fromkeys(q[0] for q in human.GRAFT_VOLUME)); size = lambda n: 6 + 0.12 * n
    for nm in names:
        pts = sorted((t, x, n) for q, n, t, x in human.GRAFT_VOLUME if q == nm); t = [p[0] for p in pts]; x = [p[1] for p in pts]
        c.plot(t, x, color=BLUE, lw=0.8, alpha=0.6, zorder=2); c.scatter(t, x, s=[size(p[2]) for p in pts], color=BLUE, edgecolor='white', linewidth=0.6, zorder=3)
    c.axhline(1, color=GRID, lw=0.5, zorder=0); c.set_xscale('log'); c.set_xlim(10, 240); c.set_ylim(0.9, 2.4); c.set_xticks([14, 30, 60, 90, 180]); c.set_xticklabels(['14', '30', '60', '90', '180']); c.minorticks_off(); c.set_yticks([1.0, 1.5, 2.0])
    c.set_xlabel('Days after transplantation'); c.set_ylabel('Graft volume / implanted volume')
    for n, yy in ((15, 2.33), (60, 2.2), (300, 2.03)): c.scatter([150], [yy], s=size(n), color=BLUE, edgecolor='white', linewidth=0.6, zorder=3); c.text(162, yy, f'$n$ = {n}', fontsize=FS_IN, va='center', ha='left')
    w = human.first_week(); rows = [('Model, calibres unchanged', w['unchanged']['pressure_load'], None, True), ('Model, calibres grow with the graft', w['proportional']['pressure_load'], None, False)]
    for lab, key in (('Measured, portal pressure ≥ 20 mmHg', 'high'), ('Measured, portal pressure < 20 mmHg', 'low')): rows.append((lab, w['observed_gradient_if_caval_pressure_8'][key], w['observed_portal_pressure'][key], None))
    for i, (lab, lo, hi, filled) in enumerate(rows):
        y = len(rows) - 1 - i
        if hi is None: d.plot([lo], [y], 'o', ms=5, mfc=BLUE if filled else 'white', mec=BLUE, mew=1.0, zorder=3); d.text(lo - 0.06, y, f'{lo:.2f}', fontsize=FS_IN, ha='right', va='center')
        else: d.plot([lo, hi], [y, y], color=INK, lw=2.2, zorder=3); d.text(lo - 0.06, y, f'{lo:.2f} to {hi:.2f}', fontsize=FS_IN, ha='right', va='center')
        d.text(-0.04, y, lab, fontsize=FS_IN, ha='right', va='center')
    d.axvline(1, color=INK, lw=LW_REF, ls=DASH, zorder=0); d.set_xlim(0, 1.1); d.set_ylim(-0.6, len(rows) - 0.3); fs.bare(d); d.set_xticks([0, 0.5, 1.0]); d.set_xlabel('Portal pressure at one week /\nportal pressure at reperfusion')
    fs.letter(fig, a, 'A', 11, 4); fs.letter(fig, b, 'B', 11, 4); fs.letter(fig, c, 'C', 11, 4); fs.letter(fig, d, 'D', 54, 4); save(fig, 'Figure_S6')


READS = {                                    # what each figure is drawn from; hashed when the figure is drawn (results/figures_provenance.json)
    '1': ('src/model/network.py', 'src/model/network3d.py', 'src/model/scaling.py'),
    '2': ('data/series.tsv', 'results/predictions.json', 'src/indices.py'),
    '3': ('data/series.tsv', 'results/bayes_main.json', 'results/spec_curve.json', 'src/indices.py'),
    '4': ('results/predictions.json', 'src/indices.py'),
    'S1': ('data/graph_scaling_all.tsv', 'data/graph3d_all.tsv', 'src/model/scaling.py'),
    'S2': ('src/model/whole_graft.py',),
    'S3': ('results/ipd_ecological.json',),
    'S4': ('src/ipd/patient_level.py', 'data/patient_level'),
    'S5': ('src/model/human.py', 'src/model/liver.py'),
    'S6': ('src/model/human.py', 'src/model/liver.py')}
FIGURES = {'1': figure_1, '2': figure_2, '3': figure_3, '4': figure_4, 'S1': figure_s1, 'S2': figure_s2, 'S3': figure_s3, 'S4': figure_s4, 'S5': figure_s5, 'S6': figure_s6}

def provenance(k):
    """Record what figure k was drawn from: the hash of every file it reads, of the two scripts that draw it, and of the
    two images written. python src/manifest.py check then says whether a figure on disk is still the one its inputs give."""
    import glob, hashlib
    sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
    def hashes(paths):
        out = {}
        for p in paths:
            for f in ([f'{ROOT}/{p}'] if os.path.isfile(f'{ROOT}/{p}') else sorted(x for x in glob.glob(f'{ROOT}/{p}/**/*', recursive=True) if os.path.isfile(x))):
                out[os.path.relpath(f, ROOT).replace(os.sep, '/')] = sha(f)
        return out
    name = f"Figure_{k}"; f = f'{OUT}/figures_provenance.json'; rec = json.load(open(f)) if os.path.exists(f) else {}
    rec[name] = dict(function=FIGURES[k].__name__, reads=hashes(READS[k]), drawn_by=hashes(('src/paper_figures.py', 'src/figstyle.py')),
                     files={e: sha(f'{OUT}/{name}.{e}') for e in ('pdf', 'png')})
    json.dump(dict(sorted(rec.items())), open(f, 'w'), indent=1)


if __name__ == '__main__':
    for k in (sys.argv[1:] or FIGURES): FIGURES[k.upper()](); provenance(k.upper())
