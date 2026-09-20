import numpy as np, pandas as pd, json, pickle, warnings; warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrow, Rectangle, FancyBboxPatch
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.interpolate import splprep, splev
from scipy.stats import spearmanr
from adjustText import adjust_text
import liver_canopy as lc, liver3d as l3
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
C_IN, C_SIN, C_OUT = '#1f5fbf', '#c62828', '#2e7d32'
def panel(ax, letter): ax.set_title(letter, loc='left', fontweight='bold', fontsize=14, pad=4)
OUT = 'results'
df = pd.read_csv(f'{OUT}/series_pdfs_z.tsv', sep='\t'); L = json.load(open(f'{OUT}/logit_zP.json'))

# ============================== FIG 1: model (hybrid) ==============================
from matplotlib.patches import FancyBboxPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], hspace=0.22, wspace=0.14, left=0.04, right=0.98, top=0.96, bottom=0.05)
def frame(ax, letter, title):
    ax.text(0.0, 1.02, letter, transform=ax.transAxes, fontweight='bold', fontsize=14, va='bottom')
    ax.text(0.5, 1.03, title, transform=ax.transAxes, ha='center', va='bottom', fontsize=7.5, fontweight='bold')
    fig.add_artist(FancyBboxPatch((ax.get_position().x0 - 0.012, ax.get_position().y0 - 0.02), ax.get_position().width + 0.024, ax.get_position().height + 0.075,
                                  boxstyle='round,pad=0.005,rounding_size=0.01', transform=fig.transFigure, fill=False, ec='0.75', lw=0.7, zorder=0))
# --- A
ax = fig.add_subplot(gs[0, 0]); frame(ax, 'A', 'Optimal canopy-to-canopy network, 2D')
G = lc.build_liver_graph(R=4.5); gin, gout = lc.split_liver(G); st = pickle.load(open(f'{OUT}/partB_states.pkl', 'rb'))['0.7500_0.0']
xy = G['xy']; cen, cor, pairs, Lh = lc.hex_lobules(4.5, 1.0)
for (x, y) in cen:
    th = np.pi / 6 + np.arange(7) * np.pi / 3; ax.plot(x + Lh * np.cos(th), y + Lh * np.sin(th), color='0.82', lw=0.4, zorder=0)
w, keep = st['tree_w'], st['tree_keep']; r = np.zeros(gin.m); r[keep] = (w[keep] * gin.L[keep]) ** 0.25; rm = r[keep].max()
outl = lc.tree_allocation(gout, np.arange(gout.m), 0.75); ro = (outl['w'] * gout.L) ** 0.25; rom = ro.max()
for k in range(gout.m):
    i, j = gout.E[k]; ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]], color=C_OUT, lw=0.2 + 2.2 * ro[k] / rom, alpha=0.5, solid_capstyle='round', zorder=1)
for k in keep:
    i, j = gin.E[k]; c = C_SIN if gin.kind[k] == 'sin' else C_IN
    ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]], color=c, lw=0.3 + 3.2 * r[k] / rm, solid_capstyle='round', zorder=2)
ax.plot(*xy[gin.src], 'ko', ms=5, zorder=5); ax.set_aspect('equal'); ax.axis('off')
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=C_IN, lw=2, label='portal tree'), Line2D([], [], color=C_SIN, lw=2, label='sinusoids (one per lobule)'),
                   Line2D([], [], color=C_OUT, lw=2, alpha=0.6, label='hepatic venous tree'), Line2D([], [], marker='o', color='k', lw=0, label='hilum')],
          loc='lower center', fontsize=5.5, ncol=4, bbox_to_anchor=(0.5, -0.07), columnspacing=0.8, handlelength=1.4)
# --- B: 3D with translucent hemisphere and axis triad
ax = fig.add_subplot(gs[0, 1], projection='3d'); ax.text2D(0.0, 1.02, 'B', transform=ax.transAxes, fontweight='bold', fontsize=14, va='bottom')
G3 = l3.build_liver3d(R=3.4); g3in, g3out = lc.split_liver(G3)
ax.text2D(0.5, 1.03, f'Hemispheric liver, 3D ({G3["n_lob"]} lobules)', transform=ax.transAxes, ha='center', va='bottom', fontsize=7.5, fontweight='bold')
fig.add_artist(FancyBboxPatch((ax.get_position().x0 - 0.012, ax.get_position().y0 - 0.02), ax.get_position().width + 0.024, ax.get_position().height + 0.075,
                              boxstyle='round,pad=0.005,rounding_size=0.01', transform=fig.transFigure, fill=False, ec='0.75', lw=0.7, zorder=0))
ft = l3.FastTree(g3in, 0.75); best = min((ft.solve(seed=s_) for s_ in range(3)), key=lambda q: q['D'])
sol = lc.tree_allocation(g3in, best['keep'], 0.75); r3 = np.zeros(g3in.m); r3[sol['keep']] = (sol['w'][sol['keep']] * g3in.L[sol['keep']]) ** 0.25; r3m = r3.max()
X = G3['xy']; segs, cols, lws = [], [], []
o3 = lc.tree_allocation(g3out, np.arange(g3out.m), 0.75); ro3 = (o3['w'] * g3out.L) ** 0.25
for k in range(g3out.m):
    i, j = g3out.E[k]; segs.append([X[i], X[j]]); cols.append('#7fbf7f'); lws.append(0.15 + 1.4 * ro3[k] / ro3.max())
for k in sol['keep']:
    i, j = g3in.E[k]; segs.append([X[i], X[j]]); cols.append(C_SIN if g3in.kind[k] == 'sin' else C_IN); lws.append(0.3 if g3in.kind[k] == 'sin' else 0.3 + 2.8 * r3[k] / r3m)
ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=lws, alpha=0.9))
Rh = 3.55; u, v = np.mgrid[0:2 * np.pi:60j, 0:np.pi / 2:30j]
ax.plot_surface(Rh * np.cos(u) * np.sin(v), Rh * np.sin(u) * np.sin(v), Rh * np.cos(v), color='0.6', alpha=0.10, linewidth=0, antialiased=True, shade=True, zorder=0)
ax.plot_wireframe(Rh * np.cos(u[::6, ::5]) * np.sin(v[::6, ::5]), Rh * np.sin(u[::6, ::5]) * np.sin(v[::6, ::5]), Rh * np.cos(v[::6, ::5]), color='0.75', lw=0.25, alpha=0.6)
th = np.linspace(0, 2 * np.pi, 100); ax.plot(Rh * np.cos(th), Rh * np.sin(th), 0, color='0.6', lw=0.6)
ax.scatter([0], [0], [0], color='k', s=16, depthshade=False)
o = np.array([-3.4, -2.4, 0.2]); L_ = 1.0
for d, lab in ((np.array([1, 0, 0]), 'x'), (np.array([0, 1, 0]), 'y'), (np.array([0, 0, 1]), 'z')):
    ax.quiver(*o, *(d * L_), color='k', arrow_length_ratio=0.25, lw=0.7); ax.text(*(o + d * L_ * 1.35), lab, fontsize=6)
ax.set_xlim(-3.6, 3.6); ax.set_ylim(-3.6, 3.6); ax.set_zlim(0, 3.6); ax.view_init(elev=26, azim=-58); ax.set_box_aspect((1, 1, 0.5), zoom=1.3); ax.set_axis_off()
ax.legend(handles=[Line2D([], [], color=C_IN, lw=2, label='portal tree'), Line2D([], [], color=C_SIN, lw=2, label='sinusoids'), Line2D([], [], color='#7fbf7f', lw=2, label='hepatic venous tree'),
                   Line2D([], [], marker='o', color='k', lw=0, label='hilum')], loc='lower right', fontsize=5.5, bbox_to_anchor=(1.02, -0.02))
# --- C: four conditions, cards + generated liver icons + vector encodings
ax = fig.add_subplot(gs[1, 0]); frame(ax, 'C', 'Representative conditions in the two indices'); ax.set_xlim(0, 4.4); ax.set_ylim(-0.75, 1.55); ax.axis('off')
conds = [('Healthy\ndonor', 1.0, 1.0, '#e6f4e6', 'donor'), ('Small graft,\nstandard\noutflow', 3.0, 3.0, '#fbe3e3', 'small'), ('Graft with\nreconstructed\noutflow', 3.5, 1.3, '#fdf0dc', 'recon'), ('Cirrhosis', 0.7, 2.5, '#ece3f6', 'cirr')]
for k, (name, zf, zp, bg, key) in enumerate(conds):
    x0 = 0.05 + k * 1.085; cx = x0 + 0.5
    ax.add_patch(FancyBboxPatch((x0, -0.5), 1.0, 2.0, boxstyle='round,pad=0.0,rounding_size=0.06', fc=bg, ec='none', zorder=0))
    ax.text(cx, 1.42, f'z$_P$ = {zp:g}', fontsize=6.5, ha='center', va='center', color=C_SIN, fontweight='bold')
    gx = cx - 0.055; ax.add_patch(Rectangle((gx, 0.88), 0.11, 0.42, fc='white', ec='0.35', lw=0.7)); ax.add_patch(Rectangle((gx, 0.88), 0.11, 0.42 * min(zp, 3.5) / 3.5, color=C_SIN, lw=0))
    ax.plot([gx - 0.04, gx + 0.15], [0.88 + 0.42 * 2 / 3.5] * 2, color='k', lw=0.6, ls='--'); ax.text(gx + 0.17, 0.88 + 0.42 * 2 / 3.5, '2', fontsize=5.5, va='center')
    img_ = mpimg.imread(f'{OUT}/liver_{key}.png'); zoom = 0.40 if key == 'donor' else (0.31 if key in ('small', 'recon') else 0.38)
    ax.add_artist(AnnotationBbox(OffsetImage(img_, zoom=zoom), (cx + 0.1, 0.45), frameon=False, zorder=3))
    ax.add_patch(FancyArrow(x0 + 0.03, 0.45, 0.14, 0, width=0.04 * zf, head_width=0.04 * zf + 0.05, head_length=0.07, color=C_IN, lw=0, zorder=4))
    ax.text(cx, 0.02, f'z$_F$ = {zf:g}', fontsize=6.5, ha='center', va='center', color=C_IN, fontweight='bold')
    ax.text(cx, -0.16, name, fontsize=6, ha='center', va='top', linespacing=1.1, fontweight='bold')
ax.text(2.2, -0.54, 'Blue arrow width: portal flow per lobule (z$_F$)\nRed gauge: sinusoidal pressure (z$_P$); dashed mark = 2', fontsize=5.3, ha='center', va='top')
# --- D
ax = fig.add_subplot(gs[1, 1]); frame(ax, 'D', 'Load per lobule on a partial graft'); gfrac = np.linspace(0.25, 1.2, 100)
for h, c in ((1.0, C_IN), (1.5, '#f28c28'), (2.0, C_SIN)): ax.plot(gfrac / 0.5, h / gfrac, color=c, lw=1.3, label=f'recipient inflow {h:g}× donor')
ax.axhspan(1, 2, color='0.9', lw=0); ax.axhline(2, color='k', ls='--', lw=0.6); ax.axhline(1, color='k', lw=0.6)
ax.text(0.53, 1.5, 'target\nwindow', fontsize=6, va='center', ha='left')
ax.set_xlabel('Graft size (≈ GRWR, %)'); ax.set_ylabel('Load per lobule, z = h/g\n(at donor outflow resistance)'); ax.set_ylim(0, 7); ax.set_xlim(0.5, 2.4); ax.legend(fontsize=6, loc='upper right')
plt.savefig(f'{OUT}/JHEP_Fig1.pdf'); plt.savefig(f'{OUT}/JHEP_Fig1.png', dpi=300); plt.close()

# ============================== FIG 2: index, risk, interventions ==============================
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1, 1, 1.15], wspace=0.5))
ax = axs[0]; panel(ax, 'A'); pvp = np.linspace(6, 26, 100)
for cvp, c in ((2, '#1f5fbf'), (5, '#2e7d32'), (8, '#f28c28'), (11, '#c62828')): ax.plot(pvp, (pvp - cvp) / 5, color=c, lw=1.2, label=f'CVP {cvp} mmHg')
ax.axhspan(1, 2, color='0.9', lw=0); ax.axhline(1, color='k', lw=0.6); ax.axhline(2, color='k', lw=0.6, ls='--')
ax.text(6.4, 1.5, 'target window', fontsize=6, va='center'); ax.text(6.4, 2.15, 'z$_P$ = 2', fontsize=6); ax.text(6.4, 0.55, 'steal risk', fontsize=6)
ax.set_xlabel('Final portal venous pressure (mmHg)'); ax.set_ylabel('z$_P$ = (PVP − CVP)/5'); ax.set_ylim(0, 4.5); ax.legend(loc='upper left', fontsize=5.5)
ax = axs[1]; panel(ax, 'B'); z = np.array(L['z'])
ax.fill_between(z, np.array(L['lo']) * 100, np.array(L['hi']) * 100, color=C_SIN, alpha=0.15, lw=0, label='95% bootstrap CI')
ax.plot(z, np.array(L['p']) * 100, color=C_SIN, lw=1.3, label='pooled logistic fit')
s = df[df.desenlace.str.contains('SFSS|disfunción') & df.z_P.notna()]
ax.scatter(s.z_P, s.pct, s=12 + s.n / 3, color='w', edgecolor=C_SIN, lw=0.8, zorder=3, label='published groups')
texts = [ax.text(r.z_P, r.pct, ' '.join(r.estudio.split()[:2]), fontsize=5) for _, r in s.iterrows()]
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.3, color='0.5'), expand=(1.4, 1.8))
ax.axvspan(1, 2, color='0.9', lw=0); ax.set_xlabel('z$_P$'); ax.set_ylabel('Small-for-size syndrome (%)'); ax.set_xlim(0.8, 4); ax.set_ylim(0, 60); ax.legend(loc='upper left', fontsize=5.5)
ax = axs[2]; panel(ax, 'C')
items = [('Hemiportocaval shunt\n(Yamada 2008, shunt open vs clamped, n=10)', -8.8, 'PVF −65% (Troisi 2005)'),
         ('Splenectomy\n(Osman 2017, n=5, median)', -6.0, ''), ('Splenectomy\n(Wang 2014, n=154)', -4.8, ''),
         ('Intentional pressure control,\nmainly splenectomy (Ogura 2010, n=77)', -4.4, ''), ('Splenic artery ligation\n(Troisi 2003, n=13)', np.nan, 'PVF −33%')]
y = np.arange(len(items))[::-1]
for yi, (n, dp, note) in zip(y, items):
    if not np.isnan(dp): ax.barh(yi, dp, 0.55, color=C_SIN, alpha=0.85)
    ax.text(0.4, yi + 0.02, n, fontsize=5.3, ha='left', va='center')
    if note: ax.text(dp - 0.3 if not np.isnan(dp) else -0.4, yi, note, fontsize=5.3, ha='right', va='center', color=C_IN)
ax.set_yticks([]); ax.axvline(0, color='k', lw=0.6); ax.spines['left'].set_visible(False)
ax.set_xlim(-11, 11); ax.set_ylim(-0.7, len(items) - 0.3); ax.set_xticks([-10, -5, 0]); ax.set_xlabel('Change in final PVP (mmHg)')
plt.savefig(f'{OUT}/JHEP_Fig2.pdf', bbox_inches='tight'); plt.savefig(f'{OUT}/JHEP_Fig2.png', dpi=300, bbox_inches='tight'); plt.close()

# ============================== FIG 3: series ==============================
fig, axs = plt.subplots(2, 1, figsize=(5.0, 8.0), gridspec_kw=dict(hspace=0.3))
sf = df.desenlace.str.contains('SFSS|disfunción')
for ax, col, lab, let in ((axs[0], 'z_P', 'z$_P$ = (PVP − CVP)/5', 'A'), (axs[1], 'z_F', 'z$_F$ = graft PVF per 100 g / donor PVF per 100 g', 'B')):
    panel(ax, let); d = df.dropna(subset=[col]); a = d[sf[d.index]]; m = d[~sf[d.index]]
    for st_, g in d.groupby('estudio'):
        if len(g) == 2 and g.desenlace.nunique() == 1: ax.plot(g[col], g.pct, '-', color='0.65', lw=0.7, zorder=1)
    ax.scatter(a[col], a.pct, s=12 + a.n.clip(upper=160) / 2.2, color=C_SIN, edgecolor='k', lw=0.4, zorder=3, label='SFSS / graft dysfunction')
    ax.scatter(m[col], m.pct, s=12 + m.n.clip(upper=160) / 2.2, color=C_IN, marker='s', edgecolor='k', lw=0.4, zorder=3, label='mortality / graft loss')
    if col == 'z_P':
        MAN = {('Ogura 2010','<15'):(1.03,24,'Ogura 2010 (<15)'),('Chan 2011',''):(1.03,18,'Chan 2011'),('Osman 2017','A'):(1.03,12,'Osman 2017 (A)'),('Yamada 2008',''):(1.45,-2.5,'Yamada 2008'),
               ('Yagi 2006',''):(1.55,14,'Yagi 2006'),('Ogura 2010','>=15'):(2.04,36,'Ogura 2010 (≥15)'),('Yagi 2005','L'):(2.04,-2.5,'Yagi 2005 (L)'),('Wang 2014','sin esplenectomía'):(2.72,24,'Wang 2014 (no splenectomy)'),('Wang 2014','esplenectomía'):(2.27,13.5,'Wang 2014 (splenectomy)'),
               ('Wang 2014','<20'):(2.27,3,'Wang 2014 (PVP <20)'),('Osman 2017','B'):(2.3,27,'Osman 2017 (B)'),('Wang 2014','>=20'):(3.23,20,'Wang 2014 (PVP ≥20)'),('Yagi 2005','H'):(3.45,45,'Yagi 2005 (H)')}
        for _, r in d.iterrows():
            for (st_, key), (tx, ty, lab_) in MAN.items():
                if r.estudio == st_ and key in r.grupo:
                    ax.annotate(lab_, (r[col], r.pct), xytext=(tx, ty), fontsize=5, ha='left', va='center', arrowprops=dict(arrowstyle='-', lw=0.3, color='0.5')); break
    else:
        texts = [ax.text(r[col], r.pct, r.estudio, fontsize=5) for _, r in d.iterrows()]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.3, color='0.5'), expand=(1.6, 2.2), force_text=(0.8, 1.4), lim=2000)
    rho, p = spearmanr(d[col], d.pct); ax.text(0.98, 0.96, f'Spearman ρ = {rho:.2f}, p = {p:.3f}\n{len(d)} groups', transform=ax.transAxes, ha='right', va='top', fontsize=6)
    ax.axvline(2, color='k', ls='--', lw=0.6); ax.set_xlabel(lab); ax.set_ylabel('Outcome (%)'); ax.set_ylim(-3, 108); ax.legend(loc='upper left', fontsize=5.5)
axs[0].axvspan(1, 2, color='0.9', lw=0); axs[0].text(2.04, 60, 'PVP 15 mmHg\nat CVP 5', fontsize=5.5); axs[1].text(2.04, 60, '250 mL/min/100 g\n(donor ref. 120)', fontsize=5.5)
plt.savefig(f'{OUT}/JHEP_Fig3.pdf', bbox_inches='tight'); plt.savefig(f'{OUT}/JHEP_Fig3.png', dpi=300, bbox_inches='tight'); plt.close()

# ============================== FIG 4: plane + resection ==============================
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw=dict(wspace=0.32))
ax = axs[0]; panel(ax, 'A')
ax.fill_between([0, 5], [0, 0], [1, 1], color='0.92', lw=0); ax.fill_betweenx([0, 5], 0, 1, color='0.92', lw=0)
ax.add_patch(Rectangle((1, 1), 1, 1, color='#2e7d32', alpha=0.18, lw=0)); ax.text(1.5, 1.22, 'target\nwindow', ha='center', va='center', fontsize=5.5)
ax.plot([0, 5], [0, 5], 'k--', lw=0.7); ax.text(4.35, 4.7, 'donor outflow resistance', fontsize=5.5, rotation=45, ha='center', va='center')
ax.plot(1, 1, 'ko', ms=5, zorder=5); ax.text(1.12, 0.62, 'donor liver', fontsize=6)
def colr(p): return '#c62828' if p > 20 else '#f28c28' if p > 5 else '#2e7d32'
for _, r in df[df.z_F.notna() & df.z_P.isna() & df.desenlace.str.contains('SFSS')].iterrows(): ax.plot(r.z_F, r.z_F, 'o', color=colr(r.pct), ms=3.5 + r.n ** 0.5 / 2, mec='k', mew=0.3, alpha=0.9, zorder=4)
for _, r in df.dropna(subset=['z_F', 'z_P']).iterrows(): ax.plot(r.z_F, r.z_P, 's', color=colr(r.pct), ms=3.5 + r.n ** 0.5 / 2, mec='k', mew=0.3, alpha=0.9, zorder=4)
ax.annotate('grafts with reconstructed\noutflow (Chan 2011,\nYagi 2006, Wang 2014)', xy=(3.8, 1.5), xytext=(2.35, 0.25), fontsize=5.5, arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4'))
ax.annotate('small grafts, standard\noutflow (Troisi, Ou;\nflow-only series on the diagonal)', xy=(2.95, 3.05), xytext=(1.1, 4.1), fontsize=5.5, arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4'))
ax.add_patch(Rectangle((0.4, 2), 0.6, 2, color='#7b3fa0', alpha=0.18, lw=0)); ax.text(0.7, 3.0, 'cirrhosis\nHVPG 10–20\nreduced\nhepatopetal\nflow', ha='center', va='center', fontsize=5.3)
ax.annotate('', xy=(0.72, 1.97), xytext=(0.98, 1.05), arrowprops=dict(arrowstyle='->', color='#7b3fa0', lw=1)); ax.text(0.03, 1.55, 'fibrosis:\nR rises', fontsize=5.3, color='#7b3fa0')
ax.annotate('', xy=(2.45, 2.45), xytext=(1.07, 1.07), arrowprops=dict(arrowstyle='->', color='0.3', lw=1)); ax.text(1.55, 2.25, 'partial graft or\nresection: × h/g', fontsize=5.5, rotation=45, ha='left', va='bottom')
ax.text(2.35, 0.06, 'grey: below normal (vessel regression, steal)', fontsize=5.3)
ax.set_xlim(0, 5); ax.set_ylim(0, 5); ax.set_xlabel('z$_F$ (portal flow per lobule / donor)'); ax.set_ylabel('z$_P$ (sinusoidal pressure / normal)'); ax.set_aspect('equal')
ax = axs[1]; panel(ax, 'B'); g = np.linspace(0.3, 1, 100)
for hv, lab, c in ((5, 'HVPG 5 mmHg (normal)', '#2e7d32'), (8, 'HVPG 8 mmHg', '#8a9a1e'), (10, 'HVPG 10 mmHg (CSPH)', '#f28c28'), (12, 'HVPG 12 mmHg', '#c62828')): ax.plot(100 * (1 - g), hv / 5 / g, color=c, lw=1.3, label=lab)
ax.axhspan(1, 2, color='0.9', lw=0); ax.axhline(2, color='k', ls='--', lw=0.7); ax.text(1, 2.1, 'z$_P$ = 2', fontsize=6)
ax.axvline(60, color='0.5', lw=0.6); ax.text(59, 4.55, 'right hepatectomy\n(≈ 60%)', fontsize=5.5, ha='right')
ax.set_xlabel('Parenchyma removed (%)'); ax.set_ylabel('z$_P$ of the remnant\nimmediately after resection'); ax.set_ylim(0, 5); ax.set_xlim(0, 70); ax.legend(fontsize=5.5, loc='upper left')
plt.savefig(f'{OUT}/JHEP_Fig4.pdf', bbox_inches='tight'); plt.savefig(f'{OUT}/JHEP_Fig4.png', dpi=300, bbox_inches='tight'); plt.close()
print('figures done')
