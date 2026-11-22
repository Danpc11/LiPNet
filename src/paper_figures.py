"""Every figure of the article, in one style (src/figstyle.py).

    python src/paper_figures.py              all eleven figures -> results/Figure_1 ... Figure_4, Figure_S1 ... Figure_S7 (.pdf and .png)
    python src/paper_figures.py 2 S3         only those

Figure 1    A the model: one lobule and the network; B, C the optimum in two and three dimensions; D load per lobule against
            graft size; E the flow-pressure plane; F the portal radii of a human cast against the model
Figure 2    A the three terms of the resistance load (diagram); B the six groups that report both loads in the flow-pressure plane,
            with the two clinical thresholds; C their resistance load against the rigid network with normal blood and with the blood
            of the recipient; D the gradient at the flow threshold against the haematocrit; E the contrast inside one cohort
Figure 3    each term where the other does not act: A to C the viscosity law in perfused livers; D to F published hepatectomies
            against the rigid network and the three relations
Figure 4    A slope per stratum; B outcome against pressure load with the fitted curve of each stratum; C specification curve;
            D the primary slope as the risk of a group against the difference of its gradient from another group of the cohort;
            E a whole graft with a narrowed outlet or inlet (a prediction)
Figure S1   A scaling of the optimum with the number of lobules; B, C allometric closure
Figure S2   inlet and outlet obstruction of a whole graft
Figure S3   slope of the integrated model against the spread of the gradient assumed within groups, and over the spreads imputed
Figure S4   each recipient of a hemiportocaval shunt
Figure S5   the model against a human vascular cast
Figure S6   growth with the body and regrowth of a graft
Figure S7   clinical thresholds on the pressure-load scale (exploratory)

Every number comes from data/, from results/*.json or from the model itself; nothing is typed in here except labels.
Figures 1B and 1C recompute the optimum (about a minute). Figures 2, 3, 4, S3 and S7 need the results of run_all.sh.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs
from figstyle import INK, MUTED, GRID, BAND, NEUTRAL, BLUE, RED, PORTAL, SINUSOID, HEPATIC, EVENT, FS, FS_IN, FS_SMALL, LW_REF, DASH, DOT
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch, FancyBboxPatch, RegularPolygon, Ellipse
from matplotlib.path import Path
import indices

fs.use()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); OUT = f'{ROOT}/results'; DATA = f'{ROOT}/data'
MODEL_ZR = (0.81, 1.0)                       # resistance load the network reaches for a partial graft (model/liver.py, two_inflows)
ZF, ZP, ZR = '$z_F$', '$z_P$', '$z_R$'

series = lambda: indices.compute(pd.read_csv(f'{DATA}/series.tsv', sep='\t'))
load = lambda name: json.load(open(f'{OUT}/{name}.json'))
save = lambda fig, name: fs.save(fig, name, OUT)


def plane(ax, lim=5.0, model=True, label_at=2.9, label='Network model'):
    """The flow-pressure plane: donor diagonal, the window between one and two and the wedge the rigid network reaches."""
    x = np.array([0, lim]); ax.fill_between([1, 2], 1, 2, color=BAND, lw=0, zorder=0)
    if model:
        ax.fill_between(x, MODEL_ZR[0] * x, MODEL_ZR[1] * x, zorder=1, **fs.MODEL_TINT)
        ax.text(0.81 * lim, 0.81 * lim * np.mean(MODEL_ZR) - 0.004 * lim, label, fontsize=FS_SMALL, rotation=42.5, ha='center', va='center', zorder=3)
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
    G = lc.build_liver_graph(R=R); gin, gout = lc.split_liver(G); xy = G['xy']; cen, _, _, Lh = lc.hex_lobules(R, 1.0)
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


def _cast(a, legend='lower left'):
    """Radius of the portal vein per generation: a human corrosion cast against the model with no fitted scale."""
    from model import human
    g, r, l, N = human.cast_arrays('portal vein'); x = np.logspace(0, np.log10(6e4), 80)
    for bb, col in zip((0.5, 0.75, 1.0), fs.ramp(3)): a.plot(x, human.radius_at(human.liver(bb, liver_g=human.CAST_LIVER_G), x), color=col, zorder=2, label=f'Model, $b$ = {bb:g}')
    a.plot(N, r, 'o', ms=4.4, mfc=INK, zorder=3, label='Human cast', **DOT); a.legend(loc=legend, handlelength=1.4)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(0.6, 8e4); a.set_ylim(0.1, 25); a.set_xticks([1, 10, 100, 1e3, 1e4]); a.set_xticklabels(['1', '10', '10$^2$', '10$^3$', '10$^4$']); a.set_yticks([0.1, 1, 10]); a.set_yticklabels(['0.1', '1', '10'])
    a.set_xlabel('Vessels in the generation'); a.set_ylabel('Portal vein radius (mm)')


def figure_1():
    fig = fs.figure(124); a = fs.axes(fig, 2, 62, 54, 56); b = fs.axes(fig, 64, 66, 52, 52); c = fs.axes(fig, 119, 63, 60, 56, projection='3d')
    d = fs.axes(fig, 12, 11, 40, 42); e = fs.axes(fig, 75, 11, 42, 42); f = fs.axes(fig, 135, 11, 42, 42)
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
    _cast(f); fs.letter(fig, a, 'A', 1.0, 3); fs.letter(fig, b, 'B', 1.0, 3); fs.letter(fig, c, 'C', -2.0, 2); fs.letter(fig, d, 'D', 10, 4); fs.letter(fig, e, 'E', 10, 4); fs.letter(fig, f, 'F', 11, 4); save(fig, 'Figure_1')


# ------------------------------------------------------------------ diagrams
def _mm(fig, left, bottom, width, height):
    """An axes without frame whose data units are millimetres."""
    ax = fs.axes(fig, left, bottom, width, height); ax.set_xlim(0, width); ax.set_ylim(0, height); ax.axis('off'); return ax


def _tube(ax, x0, x1, y, h, cells=0, arrows=0, seed=1):
    """A vessel as a rounded tube h mm across, with red cells or flow arrows inside."""
    ax.add_patch(FancyBboxPatch((x0, y - h / 2), x1 - x0, h, boxstyle=f'round,pad=0,rounding_size={0.48 * h}', fc='#fbf1f0' if cells else BAND, ec=INK, lw=0.6, zorder=2))
    rng = np.random.default_rng(seed); pts = []
    while len(pts) < cells:
        p = np.array([rng.uniform(x0 + 1.6, x1 - 1.6), rng.uniform(y - h / 2 + 0.8, y + h / 2 - 0.8)])
        if all(np.hypot(*((p - q) * (1, 1.5))) > 1.5 for q in pts): pts.append(p)
    for p in pts: ax.add_patch(Ellipse(p, 1.35, 0.85, angle=rng.uniform(-25, 25), fc=RED, ec='none', alpha=0.85, zorder=3))
    for k in range(arrows):
        yy = y + (k - (arrows - 1) / 2) * h / (arrows + 0.6)
        ax.annotate('', (x1 - 2.2, yy), (x0 + 2.2, yy), arrowprops=dict(arrowstyle='-|>', lw=0.6, color=INK, mutation_scale=4.5, shrinkA=0, shrinkB=0), zorder=3)


def _terms(ax, V):
    """The three terms of the resistance load of a graft: architecture, blood, distension."""
    S = V['clinical']['span']['normal_blood']['rigid']; f = V['clinical']['factor']; top, cap, val = 26.8, 5.6, 2.2
    for k, s in enumerate(('Pressure load =', 'flow load ×', 'resistance load')): ax.text(0, 21.4 - 3.1 * k, s, fontsize=FS_IN, ha='left', va='center')
    ax.text(0, 9.4, f'{ZP} = {ZF} × {ZR}', fontsize=FS + 1, ha='left', va='center'); ax.text(34.5, 15, f'{ZR} =', fontsize=FS + 1, ha='right', va='center')
    # architecture: two trees joined by one path per lobule
    y = 8.0 + 1.75 * np.arange(9); yc = y.mean(); xs = np.array([37.0, 41.5, 46.5, 51.0]); mirror = lambda x: 111.0 - x
    for col, g in ((PORTAL, lambda x: x), (HEPATIC, mirror)):
        ax.plot([g(xs[0]), g(xs[1])], [yc, yc], color=col, lw=3.0, solid_capstyle='butt', zorder=2)
        for k in range(3):
            ym = y[3 * k:3 * k + 3].mean(); _curve(ax, (g(xs[1]), yc), (g(xs[2]), ym), 1.8, col)
            for yy in y[3 * k:3 * k + 3]: _curve(ax, (g(xs[2]), ym), (g(xs[3]), yy), 0.9, col)
    for yy in y: ax.plot([xs[3], mirror(xs[3])], [yy, yy], color=SINUSOID, lw=0.9, solid_capstyle='butt', zorder=3)
    ax.text(55.5, top, 'Architecture', fontsize=FS, fontweight='bold', ha='center', va='center')
    ax.text(55.5, cap, 'Rigid network, normal blood', fontsize=FS_IN, ha='center', va='center'); ax.text(55.5, val, f'{S[0]:.2f} to {S[1]:.2f}', fontsize=FS_IN, color=MUTED, ha='center', va='center')
    # blood: the same vessel at two haematocrits
    ax.text(80.0, 15, '×', fontsize=FS + 3, ha='center', va='center')
    _tube(ax, 86, 112, 19.2, 5.0, cells=26, seed=3); _tube(ax, 86, 112, 12.0, 5.0, cells=14, seed=5)
    hr, hc = V['clinical']['haematocrit_reference'], V['clinical']['haematocrit_recipient']
    ax.text(113.5, 19.2, f'Normal,\n{100 * hr[0]:.0f} to {100 * hr[1]:.0f}%', fontsize=FS_SMALL, ha='left', va='center', linespacing=1.0)
    ax.text(113.5, 12.0, f'Recipient,\n{100 * hc[0]:.0f} to {100 * hc[1]:.0f}%', fontsize=FS_SMALL, ha='left', va='center', linespacing=1.0)
    ax.text(104.5, top, 'Blood', fontsize=FS, fontweight='bold', ha='center', va='center')
    ax.text(104.5, cap, 'Viscosity at the haematocrit of the recipient', fontsize=FS_IN, ha='center', va='center'); ax.text(104.5, val, f'× {f[0]:.2f} to {f[1]:.2f}', fontsize=FS_IN, color=MUTED, ha='center', va='center')
    # distension: the same vessel at donor flow and at the flow of a graft
    ax.text(131.5, 15, '×', fontsize=FS + 3, ha='center', va='center')
    _tube(ax, 137, 159, 20.2, 2.6, arrows=1); _tube(ax, 137, 159, 12.6, 5.6, arrows=3)
    ax.text(160.5, 20.2, 'Donor flow', fontsize=FS_SMALL, ha='left', va='center'); ax.text(160.5, 12.6, 'Flow of\na graft', fontsize=FS_SMALL, ha='left', va='center', linespacing=1.0)
    ax.text(153.5, top, 'Distension', fontsize=FS, fontweight='bold', ha='center', va='center')
    ax.text(153.5, cap, 'Vessels widen as pressure rises', fontsize=FS_IN, ha='center', va='center'); ax.text(153.5, val, '× less than 1, falling with flow', fontsize=FS_IN, color=MUTED, ha='center', va='center')


def _liver(ax, x, y, w, h, **kw):
    kw = dict(dict(fc='#ead9cf', ec=INK, lw=0.6, zorder=2), **kw)
    ax.add_patch(FancyBboxPatch((x + 2, y + 2), w - 4, h - 4, boxstyle='round,pad=2,rounding_size=3.2', **kw))


def _perfusion(ax):
    """A liver perfused through the portal vein at a fixed flow with blood of four haematocrits."""
    ax.text(0, 46, 'Perfused livers', fontsize=FS, fontweight='bold', ha='left', va='center')
    for k, (hct, alpha) in enumerate(((0, 0.0), (10, 0.3), (20, 0.6), (30, 0.9))):
        x = 3 + 9.5 * k; ax.add_patch(FancyBboxPatch((x + 0.6, 31.6), 4.8, 5.8, boxstyle='round,pad=0.6,rounding_size=1.2', fc='white', ec=INK, lw=0.6, zorder=2))
        ax.add_patch(FancyBboxPatch((x + 0.6, 31.6), 4.8, 3.8, boxstyle='round,pad=0.6,rounding_size=1.2', fc=RED, ec='none', alpha=alpha, zorder=3)); ax.text(x + 3, 28.4, f'{hct}%' if k == 3 else f'{hct}', fontsize=FS_SMALL, ha='center', va='center')
    ax.text(41.5, 34.5, 'Haematocrit of\nthe perfusate', fontsize=FS_SMALL, color=MUTED, ha='left', va='center', linespacing=1.0)
    _liver(ax, 17, 8, 22, 15); ax.text(28, 15.5, 'Liver', fontsize=FS_IN, ha='center', va='center', zorder=4)
    arrow = lambda p, q, col: ax.annotate('', q, p, arrowprops=dict(arrowstyle='-|>', lw=1.6, color=col, mutation_scale=7, shrinkA=0, shrinkB=0), zorder=3)
    ax.plot([6, 6, 16.5], [27, 15.5, 15.5], color=PORTAL, lw=1.6, solid_capstyle='round', solid_joinstyle='round', zorder=1); arrow((12, 15.5), (17, 15.5), PORTAL); arrow((39, 15.5), (48, 15.5), HEPATIC)
    ax.text(1, 11.8, 'Fixed\nflow', fontsize=FS_SMALL, ha='left', va='center', linespacing=1.0); ax.text(43.5, 11.2, 'Pressure\ndrop', fontsize=FS_SMALL, ha='center', va='center', linespacing=1.0)
    ax.text(0, 2.6, 'Same vessels, same flow: only the blood changes', fontsize=FS_SMALL, color=MUTED, ha='left', va='center')


def _resection(ax):
    """The same liver before and after resection: the remnant takes the whole portal inflow."""
    ax.text(0, 46, 'Hepatectomy', fontsize=FS, fontweight='bold', ha='left', va='center')
    _liver(ax, 1, 20, 22, 15); ax.text(12, 27.5, 'Before', fontsize=FS_IN, ha='center', va='center', zorder=4)
    _liver(ax, 29, 20, 22, 15, fc='white', ec=MUTED, ls=(0, (2, 1.5))); _liver(ax, 29, 20, 9, 15); ax.text(44.5, 27.5, 'Removed', fontsize=FS_SMALL, color=MUTED, ha='center', va='center', zorder=4)
    for x in (12, 33.5): ax.annotate('', (x, 20), (x, 12), arrowprops=dict(arrowstyle='-|>', lw=1.6, color=PORTAL, mutation_scale=7, shrinkA=0, shrinkB=0), zorder=3)
    ax.text(22.8, 15.6, 'Same portal inflow', fontsize=FS_SMALL, ha='center', va='center')
    ax.text(12, 39.6, f'Whole liver, {ZF} = 1', fontsize=FS_SMALL, ha='center', va='center'); ax.text(40, 39.6, f'Remnant $g$, {ZF} = 1/$g$', fontsize=FS_SMALL, ha='center', va='center')
    ax.text(0, 7.0, 'Same animal, same blood:', fontsize=FS_SMALL, color=MUTED, ha='left', va='center'); ax.text(0, 3.6, 'only the flow per gram changes', fontsize=FS_SMALL, color=MUTED, ha='left', va='center')


# ------------------------------------------------------------------ Figure 2
def label_of(r):
    return r.study + (', ' + ('splenectomy' if r.group == 'splenectomy' else 'no splenectomy') if r.study == 'Wang 2014' else '')


# the mechanistic class of a group is a symbol, never a colour: (marker, filled, label)
CONFIGURATION = {'outflow reconstruction': ('o', True, 'Outflow reconstruction'), 'combined': ('D', True, 'Outflow and pressure control'),
                 'inflow modulation': ('^', True, 'Inflow modulation (splenectomy)'), 'none': ('o', False, 'No manoeuvre')}


def figure_2():
    """A the three terms of the resistance load; B the six groups in the flow-pressure plane with the two clinical thresholds; C their
    resistance load against the rigid network with normal blood and with the blood of the recipient; D the gradient at the flow
    threshold against the haematocrit (lighter above the haematocrit at which the law was checked in perfused livers); E the contrast inside one cohort."""
    d = series(); p = d.dropna(subset=['zF', 'zP']); flow_only = d[d.zF.notna() & d.zP.isna() & d.in_main_analysis]; P = load('predictions')['P2_discordance']
    U = {(g['study'], g['group']): g for g in P['uncertainty']['groups']}; V = load('model_viscosity'); T = V['thresholds']; S = V['clinical']['span']; n0, n1 = S['normal_blood']['rigid'], S['blood_of_the_recipient']['rigid']
    fig = fs.figure(170); t = _mm(fig, 3, 137, 174, 30); a = fs.axes(fig, 15, 72, 50, 50); b = fs.axes(fig, 123, 72, 53, 50); c = fs.axes(fig, 15, 11, 50, 33); e = fs.axes(fig, 120, 11, 38, 33)
    _terms(t, V)
    # B: the plane, with the scales of the operating room on the far sides
    lim = 5.0; x = np.array([0, lim]); a.fill_between(x, n0[0] * x, n0[1] * x, facecolor='none', edgecolor=BLUE, hatch='//////', lw=0, alpha=0.55, zorder=1); a.fill_between(x, n1[0] * x, n1[1] * x, color=BLUE, alpha=0.34, lw=0, zorder=1)
    a.plot(x, x, color=INK, lw=LW_REF, ls=DASH, zorder=2); a.plot([1], [1], 'o', ms=4, mfc=INK, zorder=5, **DOT); a.text(0.72, 1.28, 'Donor', fontsize=FS_IN, ha='right')
    thr = dict(color=MUTED, lw=0.7, ls=(0, (1, 1.5)), zorder=2); a.plot([T['zF']] * 2, [0, lim], **thr); a.plot([0, lim], [T['zP']] * 2, **thr)
    a.plot(x, T['zR_where_they_agree'] * x, color=INK, lw=0.7, zorder=3)
    a.text(4.42, 4.42 * T['zR_where_they_agree'] + 0.24, f"{ZR} = {T['zR_where_they_agree']:.2f}", fontsize=FS_IN, rotation=35.8, ha='center', va='center', zorder=4); a.text(4.42, 4.66, f'{ZR} = 1', fontsize=FS_IN, rotation=45, ha='center', va='center')
    a.text(T['zF'] + 0.07, 0.34, 'Flow\nthreshold', fontsize=FS_SMALL, color=MUTED, ha='left', va='center', linespacing=1.0); a.text(0.1, T['zP'] + 0.17, 'Pressure threshold', fontsize=FS_SMALL, color=MUTED, ha='left', va='center')
    a.plot(flow_only.zF, np.zeros(len(flow_only)), '|', ms=7, mew=0.9, color=MUTED, clip_on=False, zorder=4)
    order = {(r.study, r.group): k + 1 for k, r in enumerate(sorted(p.itertuples(), key=lambda r: r.zR))}; left = {'Yagi 2006', 'Wu 2011', 'Ishizaki 2012'}
    for _, r in p.iterrows():
        mk, filled, _ = CONFIGURATION[r.configuration]; ms = (3.2 + r.n ** 0.5 / 3.2) * (0.85 if mk == 'D' else 1); dx = (0.09 + ms / 95) * (-1 if r.study in left else 1)
        a.plot(r.zF, r.zP, mk, ms=ms, mfc=INK if filled else 'white', mec=INK, mew=1.0, zorder=6)
        a.text(r.zF + dx, r.zP, str(order[(r.study, r.group)]), fontsize=FS_IN, ha='right' if r.study in left else 'left', va='center', zorder=7)
    a.set_xlim(0, lim); a.set_ylim(0, lim); a.set_aspect('equal'); a.set_xticks(range(6)); a.set_yticks(range(6)); a.set_xlabel(f'Flow load, {ZF}'); a.set_ylabel(f'Pressure load, {ZP}')
    q = indices.DEFAULT_DONOR_REF; g0 = indices.NORMAL_GRADIENT
    sx = a.secondary_xaxis('top', functions=(lambda z: z * q, lambda v: v / q)); sx.set_xticks([0, 90, 180, 250, 360, 450]); sx.set_xlabel('Portal flow (mL/min/100 g)')
    sy = a.secondary_yaxis('right', functions=(lambda z: z * g0, lambda v: v / g0)); sy.set_yticks([0, 5, 10, 15, 20, 25]); sy.set_ylabel('Portocaval gradient (mmHg)')
    key = [Line2D([], [], marker=mk, ls='', ms=4.2 if mk != 'D' else 3.8, mfc=INK if filled else 'white', mec=INK, mew=1.0, label=lab) for mk, filled, lab in CONFIGURATION.values()]
    key += [Patch(facecolor='none', edgecolor=BLUE, hatch='//////', lw=0, alpha=0.55, label='Rigid network, normal blood'), Patch(color=BLUE, alpha=0.34, lw=0, label='Rigid network, blood of the recipient')]
    fig.legend(handles=key, loc='center', ncol=3, bbox_to_anchor=(0.5, 56.5 / 170), handletextpad=0.4, handlelength=1.3, borderaxespad=0.0, labelspacing=0.45, columnspacing=2.2)
    # C: resistance load per group
    rows = sorted(p.itertuples(), key=lambda r: r.zR); top = len(rows) - 0.35
    b.fill_between(n0, -0.6, top, facecolor='none', edgecolor=BLUE, hatch='//////', lw=0, alpha=0.55, zorder=0); b.fill_between(n1, -0.6, top, color=BLUE, alpha=0.34, lw=0, zorder=0)
    b.plot([1, 1], [-0.6, top], color=INK, lw=LW_REF, ls=DASH, zorder=1); b.plot([T['zR_where_they_agree']] * 2, [-0.6, top], color=INK, lw=0.7, zorder=1)
    for k, r in enumerate(rows):
        y = len(rows) - 1 - k; g = U[(r.study, r.group)]; mk, filled, _ = CONFIGURATION[r.configuration]
        b.plot(g['ci_total'], [y, y], color=INK, lw=0.6, solid_capstyle='butt', zorder=2); b.plot(g['ci_sampling'], [y, y], color=INK, lw=2.2, solid_capstyle='butt', zorder=2); b.plot(r.zR, y, mk, ms=5 if mk != 'D' else 4.4, mfc=INK if filled else 'white', mec=INK, mew=1.0, zorder=3)
        b.text(-0.04, y + 0.13, f'{k + 1}  ' + label_of(r), fontsize=FS_IN, ha='right', va='center'); b.text(r.zR, y + 0.2, f'{r.zR:.2f}', fontsize=FS_IN, ha='center', va='bottom', zorder=4, bbox=dict(fc='white', ec='none', pad=0.4, alpha=0.75))
        b.text(-0.04, y - 0.2, f'$n$ = {int(r.n)}, {ZF} = {r.zF:.2f}', fontsize=FS_SMALL, color=MUTED, ha='right', va='center')
    b.text(T['zR_where_they_agree'], top + 0.12, 'Thresholds agree', fontsize=FS_SMALL, color=MUTED, ha='center', va='bottom')
    for x0, lw, lab in ((-0.86, 2.2, 'Sampling'), (-0.50, 0.6, 'With references')):
        b.plot([x0, x0 + 0.06], [-0.86, -0.86], color=INK, lw=lw, solid_capstyle='butt', clip_on=False, zorder=2); b.text(x0 + 0.08, -0.86, lab, fontsize=FS_SMALL, color=MUTED, ha='left', va='center')
    b.text(-0.04, top + 0.62, f"{P['cohorts_below_one']} of {P['cohorts']} cohorts below 1; exact sign test\n$P$ = {P['sign_test_cohorts_one_sided_p']:.3f} one-sided, {P['sign_test_cohorts_two_sided_p'] + 1e-9:.3f} two-sided", fontsize=FS_IN, ha='left', va='bottom', linespacing=1.15)
    b.set_xlim(0, 1.2); b.set_ylim(-0.6, len(rows) + 0.95); fs.bare(b); b.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]); b.set_xlabel(f'Resistance load, {ZR} = {ZP} / {ZF}')
    # D: what the rigid network gives at the flow threshold as the blood changes
    K = T['curve']; h = 100 * np.array(K['haematocrit']); hr, hc = V['clinical']['haematocrit_reference'], V['clinical']['haematocrit_recipient']
    c.axvspan(100 * hc[0], 100 * hc[1], color=BAND, lw=0, zorder=0); c.axvspan(100 * hr[0], 100 * hr[1], color=BAND, lw=0, zorder=0)
    lo_, hi_ = np.array(K['gradient_low']), np.array(K['gradient_high']); hk = 100 * V['published_check']['haematocrit_checked_up_to']        # the law was set against perfused livers up to here
    c.fill_between(h, lo_, hi_, where=h <= hk + 1e-9, color=BLUE, alpha=0.34, lw=0, zorder=1); c.fill_between(h, lo_, hi_, where=h >= hk - 1e-9, color=BLUE, alpha=0.16, lw=0, zorder=1); c.plot([20, 45], [T['gradient_mmHg']] * 2, **thr)
    c.text(20.4, 12.9, f'Darker: law checked in\nperfused livers, to {hk:.0f}%', fontsize=FS_SMALL, color=MUTED, ha='left', va='center', linespacing=1.1)
    c.text(30, 15.3, 'Recipient', fontsize=FS_IN, ha='center', va='center'); c.text(42.5, 15.3, 'Normal', fontsize=FS_IN, ha='center', va='center'); c.text(20.4, T['gradient_mmHg'] + 0.55, 'Pressure threshold', fontsize=FS_SMALL, color=MUTED, ha='left', va='center')
    ri = T['rise_of_the_gradient_from_0p25_to_0p35']; c.text(44.6, 4.9, f"25 to 35%: gradient +{100 * (ri[0] - 1):.0f} to {100 * (ri[1] - 1):.0f}%", fontsize=FS_IN, ha='right', va='center')
    c.set_xlim(20, 45); c.set_ylim(4, 16); c.set_xticks([20, 25, 30, 35, 40, 45]); c.set_yticks([4, 8, 12, 16]); c.set_xlabel('Haematocrit (%)'); c.set_ylabel(f"Gradient at {T['flow_mL_min_100g']:.0f} mL/min/100 g (mmHg)")
    # E: the one contrast inside a cohort, as distributions
    W = load('model_distension')['wang_contrast']; edges = np.array(W['histogram_edges']); mids = 0.5 * (edges[1:] + edges[:-1])
    rows_c = [('Observed', W['observed_histogram'], W['observed_ratio'], INK, 0.55)] + [(lab, W['models'][k]['histogram'], W['models'][k]['ratio'], BLUE, alpha)
              for k, lab, alpha in (('rigid network', 'Rigid network', 0.18), ('cube law', 'Cube law', 0.40), ('compliance', 'Compliance', 0.40), ('pressure-flow', 'Pressure-flow', 0.40))]
    e.plot([1, 1], [-0.3, len(rows_c) - 0.1], color=INK, lw=LW_REF, ls=DASH, zorder=1)
    for k, (lab, hh, qq, col, alpha) in enumerate(rows_c):
        y = len(rows_c) - 1 - k; hh = np.asarray(hh, float); hh = 0.85 * hh / hh.max()
        e.fill_between(mids, y, y + hh, color=col, alpha=alpha, lw=0, zorder=2); e.plot(mids, y + hh, color=col, lw=0.6, zorder=3); e.plot(qq[1], y, '|', ms=5, mew=1.0, color=col, zorder=4)
        e.text(0.585, y + 0.3, lab, fontsize=FS_IN, ha='right', va='center'); e.text(1.105, y + 0.3, f"{qq[1]:.2f} ({qq[0]:.2f} to {qq[2]:.2f})", fontsize=FS_SMALL, color=MUTED, ha='left', va='center')
    e.set_xlim(0.6, 1.1); e.set_ylim(-0.3, len(rows_c) - 0.1); fs.bare(e); e.set_xticks([0.6, 0.7, 0.8, 0.9, 1.0, 1.1]); e.set_xlabel(f'Ratio of {ZR}, splenectomy to none (Wang 2014)')
    fig.text(3 / 180, 168.5 / 170, 'A', fontsize=fs.FS_LETTER, fontweight='bold', va='center', ha='left', color=INK)
    fs.letter(fig, a, 'B', 12, 10); fs.letter(fig, b, 'C', 44, 10); fs.letter(fig, c, 'D', 12, 4); fs.letter(fig, e, 'E', 19, 4); save(fig, 'Figure_2')


# ------------------------------------------------------------------ Figure 3
STUDY_MARK = {'Golriz 2016': 'o', 'Golriz 2022': 's', 'Bucur 2017': '^', 'Mehrabi 2026': 'D'}


def figure_3():
    """Each term where the other does not act. A to C the viscosity law in livers perfused at a fixed flow; D to F the same
    liver before and after resection against the rigid network and the three relations."""
    H = load('model_hepatectomy'); V = load('model_viscosity'); C = V['published_check']; S = C['states']
    fig = fs.figure(122); da = _mm(fig, 3, 66, 52, 50); b = fs.axes(fig, 72, 72, 40, 40); c = fs.axes(fig, 141, 72, 35, 40)
    dd = _mm(fig, 3, 6, 52, 50); e = fs.axes(fig, 72, 12, 46, 42); f = fs.axes(fig, 141, 12, 35, 42)
    _perfusion(da); _resection(dd)
    # B: total resistance against the haematocrit of the perfusate
    K = C['curve']; b.plot([100 * h for h in K['haematocrit']], K['total_at_rat_size'], color=BLUE, zorder=2, label='Viscosity law')
    names = ['haematocrit 0 over 0.30', 'haematocrit 0.1 over 0.30', 'haematocrit 0.2 over 0.30']
    for h, n in zip((0.0, 0.10, 0.20), names):
        o = S[n]['observed_by_flow']['total']; b.plot(100 * h + np.linspace(-1.3, 1.3, len(o)), o, 'o', ms=3.6, mfc=INK, mec='white', mew=0.5, zorder=3, label='Livers, five flows' if h == 0 else None)
    b.plot(30, 1.0, 'o', ms=3.6, mfc=INK, mec='white', mew=0.5, zorder=3)
    b.set_xlim(-4, 34); b.set_ylim(0, 1.1); b.set_xticks([0, 10, 20, 30]); b.set_yticks([0, 0.25, 0.5, 0.75, 1.0]); b.set_xlabel('Haematocrit of the perfusate (%)'); b.set_ylabel('Resistance relative to 30%')
    b.legend(loc='lower right', handletextpad=0.4, handlelength=1.3, borderaxespad=0.1)
    # C: by segment, without red cells against 30%
    z = S['haematocrit 0 over 0.30']; segs = [('total', 'Whole liver'), ('portal', 'Portal'), ('sinusoidal', 'Sinusoidal'), ('hepatic_venous', 'Hepatic venous')]
    for k, (key, lab) in enumerate(segs):
        y = len(segs) - 1 - k; o = z['observed_by_flow'][key]; pr = z['predicted']['network at rat size'][key]
        c.plot(o, y + np.linspace(-0.16, 0.16, len(o)), 'o', ms=3.6, mfc=INK, mec='white', mew=0.5, zorder=3); c.plot([pr, pr], [y - 0.32, y + 0.32], color=BLUE, lw=1.6, solid_capstyle='butt', zorder=2)
        c.text(-0.03, y, lab, fontsize=FS_IN, ha='right', va='center')
    c.legend(handles=[Line2D([], [], marker='|', ls='', ms=7, mew=1.6, color=BLUE, label='Viscosity law'), Line2D([], [], marker='o', ls='', ms=3.6, mfc=INK, mec='white', mew=0.5, label='Livers, five flows')], loc='lower left', handletextpad=0.3, handlelength=1.0, borderaxespad=0.0, labelspacing=0.25, fontsize=FS_SMALL)
    c.set_xlim(0, 0.8); c.set_ylim(-0.95, len(segs) - 0.4); fs.bare(c); c.set_xticks([0, 0.2, 0.4, 0.6, 0.8]); c.set_xlabel('Resistance without red cells\nrelative to 30%')
    # E: resistance load against flow load, the same liver before and after
    B = H['band']; zz = np.array(B['zF']); lo = lambda m: np.array([x[0] for x in B[m]]); hi = lambda m: np.array([x[1] for x in B[m]])
    e.fill_between(zz, lo('rigid network'), hi('rigid network'), zorder=0, **fs.MODEL_TINT)
    e.fill_between(zz, np.minimum(lo('cube law'), lo('compliance')), np.maximum(hi('cube law'), hi('compliance')), color=BLUE, alpha=0.34, lw=0, zorder=0)
    e.fill_between(zz, lo('pressure-flow'), hi('pressure-flow'), facecolor='none', edgecolor=BLUE, hatch='//////', lw=0, alpha=0.55, zorder=0)
    e.plot([0.25, 12], [1, 1], color=INK, lw=LW_REF, ls=DASH, zorder=1)
    mu_ = H['relations_measured_up_to_flow_load']; e.plot([mu_, mu_], [0.1, 0.62], color=MUTED, lw=0.7, ls=(0, (1, 1.5)), zorder=1)       # beyond it the bands of the relations are extrapolations
    e.text(mu_ / 1.045, 0.107, 'Relations measured', fontsize=FS_SMALL, color=MUTED, rotation=90, ha='right', va='bottom')
    for x in H['states']:
        mk = STUDY_MARK[x['study']]; filled = x['shunt'] == 'none'
        e.plot([x['zF'], x['zF']], x['zR_observed_ci'], color=INK, lw=0.7, solid_capstyle='butt', zorder=3); e.plot(x['zF'], x['zR_observed'], mk, ms=4.6 if mk != 'D' else 4.0, mfc=INK if filled else 'white', mec=INK, mew=0.9, zorder=5)
    e.set_xscale('log'); e.set_yscale('log'); e.set_xlim(0.25, 12); e.set_ylim(0.1, 3.4)
    e.set_xticks([0.3, 1, 3, 10]); e.set_xticklabels(['0.3', '1', '3', '10']); e.set_yticks([0.1, 0.3, 1, 3]); e.set_yticklabels(['0.1', '0.3', '1', '3']); e.minorticks_off()
    e.set_xlabel(f'Flow load, {ZF}'); e.set_ylabel(f'Resistance load, {ZR}')
    T = H['trend']; ro, po = T['rho_one_state_per_group'], T['p_two_sided_one_state_per_group']
    for yy, s_, kw in ((2.95, f"Seven resections: $\\rho$ = {T['spearman_rho']:.2f}", dict(fontsize=FS_IN)), (2.36, 'One per group of animals:', dict(fontsize=FS_SMALL, color=MUTED)),
                       (1.93, f"$\\rho$ = {ro[0]:.2f} to {ro[1]:.2f}", dict(fontsize=FS_SMALL, color=MUTED)), (1.58, f"$P$ = {po[0]:.3f} to {po[1]:.3f}", dict(fontsize=FS_SMALL, color=MUTED))):
        e.text(11.3, yy, s_, ha='right', va='center', **kw)
    key = [Line2D([], [], marker=mk, ls='', ms=4.0 if mk != 'D' else 3.6, mfc=INK, mec=INK, mew=0.9, label=s) for s, mk in STUDY_MARK.items()] + [Line2D([], [], marker='o', ls='', ms=4.0, mfc='white', mec=INK, mew=0.9, label='With shunt')]
    key += [Patch(color=BLUE, alpha=0.16, lw=0, label='Rigid network'), Patch(facecolor='none', edgecolor=BLUE, hatch='//////', lw=0, alpha=0.55, label='Stops distending'), Patch(color=BLUE, alpha=0.34, lw=0, label='Keeps distending')]
    e.legend(handles=key, loc='lower left', handletextpad=0.35, handlelength=1.1, borderaxespad=0.0, labelspacing=0.18, fontsize=FS_SMALL)
    # F: observed over the range of each model, one bar per resection: no central prediction is drawn
    models = [('rigid network', 'Rigid network'), ('pressure-flow', 'Stops distending'), ('compliance', 'Compliance'), ('cube law', 'Cube law')]; res = [x for x in H['states'] if x['shunt'] == 'none']; ER = H['errors']['by_configuration']
    f.plot([1, 1], [-0.6, len(models) - 0.3], color=INK, lw=LW_REF, ls=DASH, zorder=1)
    for k, (m, lab) in enumerate(models):
        y = len(models) - 1 - k; f.plot([1, 1], [y + 0.14, y + 0.56], color='white', lw=3, zorder=2)
        for j, x in enumerate(sorted(res, key=lambda x: x['zF'])):
            yy = y + 0.04 - 0.062 * j; f.plot(x['observed_over_model'][m], [yy, yy], color=INK, lw=1.5, solid_capstyle='round', zorder=3)
        f.text(0.165, y + 0.47, lab, fontsize=FS_IN, ha='left', va='center'); f.text(0.165, y + 0.25, f"r.m.s. log error {ER[m][0]:.2f} to {ER[m][1]:.2f}", fontsize=FS_SMALL, color=MUTED, ha='left', va='center')
    f.set_xscale('log'); f.set_xlim(0.16, 2.2); f.set_ylim(-0.6, len(models) - 0.3); fs.bare(f); f.set_xticks([0.2, 0.5, 1, 2]); f.set_xticklabels(['0.2', '0.5', '1', '2']); f.minorticks_off(); f.set_xlabel(f'Observed {ZR} over model range')
    for ax, s, y in ((da, 'A', 120.5), (dd, 'D', 60.5)): fig.text(3 / 180, y / 122, s, fontsize=fs.FS_LETTER, fontweight='bold', va='center', ha='left', color=INK)
    fs.letter(fig, b, 'B', 12, 4); fs.letter(fig, c, 'C', 22, 4); fs.letter(fig, e, 'E', 12, 4); fs.letter(fig, f, 'F', 4, 4); save(fig, 'Figure_3')


# ------------------------------------------------------------------ Figure 4
def _trajectories(c, out, inn, h, small=FS_SMALL):
    """Narrowing of the outlet and of the inlet of a whole graft in the flow-pressure plane."""
    plane(c, lim=3.0, model=False, label_at=2.45)
    for rows, col in ((out, HEPATIC), (inn, PORTAL)):
        xf, yp = [r['zF'] for r in rows], [r['zP'] for r in rows]; c.plot(xf, yp, color=col, lw=1.5, zorder=4); c.plot(xf[-1], yp[-1], 'o', ms=4.5, mfc=col, zorder=5, **DOT)
    lead = dict(arrowstyle='-', color=MUTED, lw=0.4, shrinkA=1, shrinkB=2); c.plot(h, h, 'o', ms=4.5, mfc=INK, zorder=6, **DOT)
    c.annotate(f'Normal anastomoses,\ninflow {h:g} × donor', (h, h), (1.5, 0.55), fontsize=small, ha='left', va='center', arrowprops=lead)
    c.annotate('Outlet at half\nits diameter', (out[-1]['zF'], out[-1]['zP']), (0.12, 2.6), fontsize=small, ha='left', va='center', arrowprops=lead)
    c.annotate('Inlet at half\nits diameter', (inn[-1]['zF'], inn[-1]['zP']), (0.12, 1.85), fontsize=small, ha='left', va='center', arrowprops=lead)
    c.text(0.92, 1.12, 'Donor', fontsize=small, ha='right', va='bottom')


def largest_contrast(d, strata):
    """Largest difference of the portocaval gradient (mmHg) between two groups of one stratum of a fit."""
    m = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]; m = m.assign(stratum=m.study + ' / ' + m.outcome_type); m = m[m.stratum.isin(strata)]
    return float(indices.NORMAL_GRADIENT * (m.groupby('stratum').zP.max() - m.groupby('stratum').zP.min()).max())


def figure_4(h=1.5):
    """A slope per stratum; B outcome against pressure load, each curve dashed beyond the groups of its stratum; C specification curve; D the primary slope as the risk of a group against
    the difference of its gradient from another group of the cohort, dashed beyond the largest difference the cohorts hold; E a whole
    graft with a narrowed outlet or inlet (a prediction no data here test)."""
    from model import whole_graft as wg
    BM = load('bayes_main'); E = BM['exploratory (all outcomes)']; PR = BM['primary (SFSS or early dysfunction)']; S = load('spec_curve'); d = series(); name = lambda s: s.split(' / ')[0]
    strata = sorted(E['strata'], key=lambda s: E['beta_by_stratum'][s]['mean']); centre = {r.study: r.centre_id for r in d.itertuples()}
    hue = dict(zip(sorted({centre[name(s)] for s in strata}), fs.CATEGORICAL)); shapes = {}                # colour is the centre, the symbol the publication
    for s in sorted(strata): shapes[s] = 'os^Dv'[sum(centre[name(t)] == centre[name(s)] for t in shapes)]
    col = {s: hue[centre[name(s)]] for s in strata}
    fig = fs.figure(132); a = fs.axes(fig, 19, 76, 31, 49); b = fs.axes(fig, 66, 76, 48, 49); c = fs.axes(fig, 152, 76, 25, 49); dx = fs.axes(fig, 19, 12, 60, 45); ex = fs.axes(fig, 124, 12, 45, 45)
    top = len(strata) - 0.45; a.fill_betweenx([-0.55, top], *E['mu_beta_ci'], color=BAND, lw=0, zorder=0); a.plot([E['mu_beta']] * 2, [-0.55, top], color=INK, lw=0.9, zorder=1, solid_capstyle='butt')
    for k, s in enumerate(strata):
        v = E['beta_by_stratum'][s]; a.plot(v['ci'], [k, k], color=col[s], zorder=2); a.plot(v['mean'], k, shapes[s], ms=4.8, mfc=col[s], zorder=3, **DOT)
        a.text(-0.62, k + 0.14, name(s), fontsize=FS_IN, ha='right', va='center'); a.text(-0.62, k - 0.2, centre[name(s)], fontsize=FS_SMALL, color=MUTED, ha='right', va='center')
    pp = E['prob_mu_positive']; new = E['prediction_interval']; a.plot(new, [-0.95, -0.95], color=INK, lw=0.9, zorder=2, solid_capstyle='butt'); a.plot(new, [-0.95, -0.95], '|', ms=4, mew=0.9, color=INK, zorder=2)
    a.text(-0.62, -0.95, 'New centre', fontsize=FS_IN, ha='right', va='center')
    a.text(E['mu_beta'], len(strata) - 0.3, f"$\\mu_\\beta$ = {E['mu_beta']:+.2f}\n$P$($\\mu_\\beta$ > 0) " + ('> 0.999' if pp > 0.999 else f'= {pp:.3f}'), fontsize=FS_IN, ha='center', va='bottom', linespacing=1.15)
    a.set_xlim(-0.5, 3.5); a.set_ylim(-1.4, len(strata) + 0.75); fs.bare(a); a.set_xticks([0, 1, 2, 3]); a.set_xlabel(f'Slope per unit of {ZP}')
    m = d[d.in_main_analysis & d.zP.notna() & d.events.notna()]; m = m.assign(stratum=m.study + ' / ' + m.outcome_type); m = m[m.stratum.isin(E['strata'])]
    b.axvspan(1, 2, color=BAND, lw=0, zorder=0)
    for s in E['strata']:
        g = m[m.stratum == s]; z = np.linspace(max(0.8, g.zP.min() - 0.35), min(4.0, g.zP.max() + 0.35), 60); zin = np.linspace(g.zP.min(), g.zP.max(), 40)
        curve = lambda v: 100 / (1 + np.exp(-(E['alphas'][s] + E['beta_by_stratum'][s]['mean'] * v)))
        b.plot(z, curve(z), color=col[s], lw=0.9, ls=(0, (2.2, 1.6)), zorder=2); b.plot(zin, curve(zin), color=col[s], zorder=2)      # solid between the groups of the stratum, dashed beyond them
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
    # D: the primary slope read as a risk inside one cohort
    g0 = indices.NORMAL_GRADIENT; dG = np.linspace(-8, 8, 161); logit = lambda q: np.log(q / (1 - q)); risk = lambda q, beta: 100 / (1 + np.exp(-(logit(q) + beta * dG / g0)))
    span = largest_contrast(d, PR['strata']); inside = np.abs(dG) <= span + 1e-9                           # the largest difference of gradient between two groups of one cohort
    dx.axvline(0, color=MUTED, lw=0.5, ls=DASH, zorder=1); dx.fill_between(dG, risk(0.10, PR['mu_beta_ci'][0]), risk(0.10, PR['mu_beta_ci'][1]), color=BAND, lw=0, zorder=0)
    for q, colr in zip((0.05, 0.10, 0.20), fs.ramp(3)):
        y = risk(q, PR['mu_beta']); dx.plot(dG, y, color=colr, lw=0.9, ls=(0, (2.2, 1.6)), zorder=3); dx.plot(dG[inside], y[inside], color=colr, zorder=4); dx.text(8.25, y[-1], f'{100 * q:.0f}%', fontsize=FS_IN, ha='left', va='center')
    dx.text(-7.6, 73, 'Risk of the reference\ngroup: 5, 10 or 20%', fontsize=FS_IN, ha='left', va='center', linespacing=1.15); dx.text(-7.6, 57, 'Shading, 95% credible\ninterval at 10%', fontsize=FS_SMALL, color=MUTED, ha='left', va='center', linespacing=1.15)
    dx.text(-7.6, 41.5, f'Dashed, beyond the\nlargest difference\nobserved, {span:.1f} mmHg', fontsize=FS_SMALL, color=MUTED, ha='left', va='center', linespacing=1.15)
    dx.set_xlim(-8, 8); dx.set_ylim(0, 80); dx.set_yticks([0, 20, 40, 60, 80]); dx.set_xticks([-8, -4, 0, 4, 8]); dx.set_xlabel('Difference of gradient between groups of a cohort (mmHg)'); dx.set_ylabel('Small-for-size syndrome\nor early dysfunction (%)')
    # E: a whole graft
    dd = np.linspace(1.0, 0.5, 26); _trajectories(ex, wg.sweep('out', d=dd, h=h), wg.sweep('in', d=dd, h=h), h, small=FS_IN)
    fs.letter(fig, a, 'A', 17, 4); fs.letter(fig, b, 'B', 10, 4); fs.letter(fig, c, 'C', 36, 4); fs.letter(fig, dx, 'D', 17, 4); fs.letter(fig, ex, 'E', 11, 4); save(fig, 'Figure_4')


def figure_s7():
    """Clinical thresholds on the pressure-load scale, with the reference gradient at 4, 5 and 6 mmHg. Exploratory."""
    P = load('predictions')['P4_thresholds']; ref = P['reference_sensitivity']
    rows = [('Clinically significant portal hypertension', 'HVPG 10 mmHg', 'cirrhosis_HVPG_10'), ('Partial graft', 'PVP 15 mmHg at CVP 5 mmHg', 'graft_PVP15_CVP5'),
            ('Resection limit in cirrhosis', 'HVPG 10 mmHg', 'resection_HVPG_10'), ('Variceal bleeding', 'HVPG 12 mmHg', 'variceal_HVPG_12')]
    fig = fs.figure(58, fs.SINGLE); a = fs.axes(fig, 51, 12, 34, 40)
    a.axvspan(1, 2, color=BAND, lw=0, zorder=0)
    for k, (lab, sub_, key) in enumerate(rows):
        y = len(rows) - 1 - k; lo, mid, hi = ref['6 mmHg reference'][key], ref['5 mmHg reference'][key], ref['4 mmHg reference'][key]
        a.plot([lo, hi], [y, y], color=BLUE, zorder=2); a.plot(mid, y, 'o', ms=5, mfc=BLUE, zorder=3, **DOT)
        a.text(mid, y + 0.2, f'{mid:.1f}', fontsize=FS_IN, ha='center', va='bottom')
        a.text(-0.12, y + 0.12, lab, fontsize=FS_IN, ha='right', va='center'); a.text(-0.12, y - 0.2, sub_, fontsize=FS_SMALL, color=MUTED, ha='right', va='center')
    a.set_xlim(0, 3.5); a.set_ylim(-0.6, len(rows) - 0.3); fs.bare(a); a.set_xticks([0, 1, 2, 3]); a.set_xlabel(f'Pressure load, {ZP}'); save(fig, 'Figure_S7')


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
    fig = fs.figure(62); a = fs.axes(fig, 2, 6, 78, 52); b = fs.axes(fig, 106, 12, 70, 44)
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
    fs.letter(fig, a, 'A', 0, 0.5); fs.letter(fig, b, 'B', 10, 4); save(fig, 'Figure_S2')


# ------------------------------------------------------------------ Figures S3 and S4
def figure_s3():
    """A the slope of the integrated model against one spread of the gradient assumed for every group; B the same model with the
    reported spreads and, for the series that report none, each value of the grid."""
    E = load('ipd_ecological'); sw = E['sweep']; x = np.array([s['sigma_mmHg'] for s in sw]); y = np.array([s['mu_beta'] for s in sw]); ci = np.array([s['mu_beta_ci'] for s in sw])
    fig = fs.figure(64); ax = fs.axes(fig, 14, 12, 70, 46); bx = fs.axes(fig, 108, 12, 66, 46)
    for sd in E['reported_within_group_sd_mmHg']: ax.axvline(sd, color=NEUTRAL, lw=0.7, zorder=0)
    ax.fill_between(x, ci[:, 0], ci[:, 1], zorder=1, **fs.MODEL_TINT); ax.axhline(0, color=INK, lw=0.5); ax.axhline(y[0], color=INK, lw=LW_REF, ls=DASH, zorder=2)
    ax.plot(x, y, color=BLUE, zorder=3); ax.plot(x, y, 'o', ms=4.5, mfc=BLUE, zorder=4, **DOT)
    ax.text(0.08, y[0] - 0.09, f'Group-mean estimate, {y[0]:+.2f}', fontsize=FS_SMALL, ha='left', va='top'); ax.text(np.mean(E['reported_within_group_sd_mmHg']), 4.5, 'Dispersions reported', fontsize=FS_SMALL, color=MUTED, ha='center', va='top')
    pp = min(s.get('prob_positive', s.get('prob_mu_positive', 1.0)) for s in sw); ax.text(0.08, 0.5, f'$P$($\\mu_\\beta$ > 0) > {np.floor(pp * 100) / 100:.2f} at every dispersion', fontsize=FS_SMALL, ha='left', va='bottom')
    ax.text(0.08, ci[0, 1] - 0.08, '95% credible interval', fontsize=FS_SMALL, color=MUTED, ha='left', va='top')
    ax.set_xlim(-0.2, 6.2); ax.set_ylim(-0.15, 4.6); ax.set_xticks(range(0, 7)); ax.set_yticks([0, 1, 2, 3, 4]); ax.set_xlabel('SD of the gradient assumed within every group (mmHg)'); ax.set_ylabel('Slope of the integrated model, $\\mu_\\beta$')
    # B: the values given to the series that report no dispersion, each varied on its own
    I = E['imputation']; first, second = I['series_without_dispersion']; vals = I['values_mmHg']; cols = fs.ramp(len(vals))
    bx.axhline(0, color=INK, lw=0.5); bx.axhline(y[0], color=INK, lw=LW_REF, ls=DASH, zorder=2)
    for f in I['fits']:
        i, j = vals.index(f['sd_mmHg'][second]), vals.index(f['sd_mmHg'][first]); xx = i + 0.17 * (j - (len(vals) - 1) / 2)
        bx.plot([xx, xx], f['mu_beta_ci'], color=cols[j], lw=1.0, solid_capstyle='butt', zorder=3); bx.plot(xx, f['mu_beta'], 'o', ms=4.2, mfc=cols[j], zorder=4, **DOT)
    bx.legend(handles=[Line2D([], [], marker='o', ls='', ms=4.2, mfc=cols[j], label=f'{v:g}', **DOT) for j, v in enumerate(vals)], title=f'SD given to {first} (mmHg)', title_fontsize=FS_SMALL, fontsize=FS_SMALL,
              loc='upper center', ncol=len(vals), columnspacing=1.0, handletextpad=0.1, borderaxespad=0.0)
    bx.text(-0.42, 3.5, f"{I['groups_imputed']} of {I['groups']} groups report no dispersion; $P$($\\mu_\\beta$ > 0) > {np.floor(I['prob_positive_min'] * 100) / 100:.2f} in every fit", fontsize=FS_SMALL, ha='left', va='bottom')
    bx.set_xlim(-0.5, len(vals) - 0.5); bx.set_ylim(-0.15, 4.6); bx.set_xticks(range(len(vals))); bx.set_xticklabels([f'{v:g}' for v in vals]); bx.set_yticks([0, 1, 2, 3, 4]); bx.set_xlabel(f'SD of the gradient given to {second} (mmHg)'); bx.set_ylabel('Slope of the integrated model, $\\mu_\\beta$')
    fs.letter(fig, ax, 'A', 11, 3); fs.letter(fig, bx, 'B', 11, 3); save(fig, 'Figure_S3')


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
    _cast(a)
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
    for k, (v, (_, dia, _)) in enumerate(zip(human.measured_growth()['liver_volume_mL'], human.PORTAL_BY_WEIGHT)):
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
    '1': ('src/model/network.py', 'src/model/network3d.py', 'src/model/scaling.py', 'src/model/human.py', 'src/model/liver.py'),
    '2': ('data/series.tsv', 'results/predictions.json', 'results/model_distension.json', 'results/model_viscosity.json', 'src/indices.py'),
    '3': ('data/hepatectomy.tsv', 'data/kamikado2004.tsv', 'results/model_hepatectomy.json', 'results/model_viscosity.json'),
    '4': ('data/series.tsv', 'results/bayes_main.json', 'results/spec_curve.json', 'src/model/whole_graft.py', 'src/indices.py'),
    'S1': ('data/graph_scaling_all.tsv', 'data/graph3d_all.tsv', 'src/model/scaling.py'),
    'S2': ('src/model/whole_graft.py',),
    'S3': ('results/ipd_ecological.json',),
    'S4': ('src/ipd/patient_level.py', 'data/patient_level'),
    'S5': ('src/model/human.py', 'src/model/liver.py'),
    'S6': ('src/model/human.py', 'src/model/liver.py'),
    'S7': ('results/predictions.json', 'src/indices.py')}
FIGURES = {'1': figure_1, '2': figure_2, '3': figure_3, '4': figure_4, 'S1': figure_s1, 'S2': figure_s2, 'S3': figure_s3, 'S4': figure_s4, 'S5': figure_s5, 'S6': figure_s6, 'S7': figure_s7}

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
