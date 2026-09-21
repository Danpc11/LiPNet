"""Each plot is a separate function and a separate file.  Run from the repository root:

    python src/plots.py all
    python src/plots.py nomogram risk_curve series_pressure series_flow

Available plots (results/<name>.pdf and .png):
    nomogram           zP against final PVP for several CVP values
    within_study       outcome vs zP by study: common slope, one intercept per stratum
    forest             per-stratum slope with CI and the random-effects pooled estimate
    centred            every group on one scale: odds and gradient relative to its own cohort mean
    risk_change        odds ratio and absolute risk implied by a change in gradient, by baseline rate
    series_pressure    outcome (%) against zP, all groups with a pressure index
    series_flow        outcome (%) against zF, all groups with a flow index
    plane              flow-pressure plane: donor diagonal, window, published groups, cirrhosis region
    resection          zP of the remnant after resection for several baseline HVPG values
    load_curve         per-lobule load z = h/g against graft size
    interventions      change in final PVP produced by each inflow-modulation manoeuvre in the series
    network_2d         dissipation-optimal perfusion network on a 2D lobule field (recomputed, ~20 s)
    network_3d         the same on a hemispheric liver (recomputed, ~30 s)
    scaling            exponent of D*/F2 with the number of lobules, 2D and 3D graphs vs closed form (cached campaign data)
    allometry          allometric closure across species: set-point and vascular-mass exponents vs b
    shear_profile      wall shear along the three human hepatic trees from measured diameter ratios
    certificate        heuristic vs exhaustive optimum on a 7-lobule domain (cached)
"""
import json, os, sys, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))

ROOT = os.path.join(os.path.dirname(__file__), '..'); DATA = f'{ROOT}/data'; OUT = f'{ROOT}/results'
plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
C_IN, C_SIN, C_OUT = '#1f5fbf', '#c62828', '#2e7d32'
import indices
series = lambda: indices.compute(pd.read_csv(f'{DATA}/series.tsv', sep='\t'))   # indices always recomputed from the raw columns
def save(fig, name):
    os.makedirs(OUT, exist_ok=True); fig.savefig(f'{OUT}/{name}.pdf', bbox_inches='tight'); fig.savefig(f'{OUT}/{name}.png', dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', name)
import indices
def fit():
    f = f'{OUT}/within_study_fit.json'
    if not os.path.exists(f): indices.main()
    return json.load(open(f))

# ------------------------------------------------------------------ clinical plots
def nomogram():
    fig, ax = plt.subplots(figsize=(3.6, 3.2)); pvp = np.linspace(6, 26, 100)
    for cvp, c in ((2, '#1f5fbf'), (5, '#2e7d32'), (8, '#f28c28'), (11, '#c62828')): ax.plot(pvp, (pvp - cvp) / 5, color=c, lw=1.3, label=f'CVP {cvp} mmHg')
    ax.axhspan(1, 2, color='0.9', lw=0); ax.axhline(1, color='k', lw=0.6); ax.axhline(2, color='k', lw=0.6, ls='--')
    ax.text(25.8, 1.5, 'target window', fontsize=7, va='center', ha='right'); ax.text(25.8, 2.15, 'z$_P$ = 2', fontsize=7, ha='right'); ax.text(6.4, 0.55, 'steal risk', fontsize=7)
    ax.set_xlabel('Final portal venous pressure (mmHg)'); ax.set_ylabel('z$_P$ = (PVP − CVP)/5'); ax.set_ylim(0, 4.5); ax.legend(loc='upper left', fontsize=7); save(fig, 'nomogram')

def within_study():
    """Observed outcome against zP, one curve per stratum: same slope, different level."""
    L = fit(); fig, ax = plt.subplots(figsize=(4.6, 3.6))
    cols = dict(zip(L['studies'], ['#1f5fbf', '#c62828', '#2e7d32', '#f28c28', '#7b3fa0']))
    lab = lambda st: st.replace(' / SFSS_or_dysfunction', '').replace(' / mortality_or_graft_loss', '')
    zz = np.linspace(0.8, 4, 60)
    for st, a in L['alphas'].items():
        ax.plot(zz, 100 / (1 + np.exp(-(a + L['beta'] * zz))), color=cols[st], lw=1.4, alpha=0.9)
    for g in L['study_groups']:
        ax.scatter(g['zP'], g['pct'], s=14 + g['n'] / 4, color=cols[g['stratum']], edgecolor='k', lw=0.4, zorder=3)
    for st in L['studies']:                       # label each curve at its right end instead of a legend box
        gs = [g['zP'] for g in L['study_groups'] if g['stratum'] == st]
        ze = min(3.95, max(gs) + 0.35); pe = 100 / (1 + np.exp(-(L['alphas'][st] + L['beta'] * ze)))
        if pe > 57: ze = (np.log(0.57 / 0.43) - L['alphas'][st]) / L['beta']; pe = 57
        ax.annotate(lab(st), (ze, pe), xytext=(4, 0), textcoords='offset points', color=cols[st], fontsize=7, va='center')
    ax.axvspan(1, 2, color='0.92', lw=0)
    ax.set_xlabel('z$_P$'); ax.set_ylabel('Outcome (%)'); ax.set_xlim(0.8, 4); ax.set_ylim(0, 60)
    ax.set_title('Same slope, different level', fontsize=8, loc='left')
    save(fig, 'within_study')

def forest():
    """Per-stratum slope with its confidence interval and the random-effects pooled estimate."""
    import json as _json
    f = f'{OUT}/meta_slope.json'
    if not os.path.exists(f): indices.main()
    M = _json.load(open(f))['meta']
    rows = sorted(M['strata'], key=lambda r: r['beta'])
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    lab = lambda s: s.replace(' / SFSS_or_dysfunction', ' (SFSS)').replace(' / mortality_or_graft_loss', ' (mortality)')
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        lo, hi = r['beta'] - 1.96 * r['se'], r['beta'] + 1.96 * r['se']
        ax.plot([np.exp(lo), np.exp(hi)], [yi, yi], color='0.35', lw=1.2)
        ax.plot(np.exp(r['beta']), yi, 's', color=C_IN, ms=4 + min(10, r['events'] ** 0.5), mec='k', mew=0.4)
        ax.text(0.9, yi + 0.33, f"{lab(r['stratum'])}  ({r['events']}/{r['patients']}, z$_P$ span {r['zP_span']:.2f})", fontsize=6.5, va='bottom', ha='left')
    ax.axvline(np.exp(M['beta']), color=C_SIN, lw=1.2)
    ax.axvspan(np.exp(M['beta_ci'][0]), np.exp(M['beta_ci'][1]), color=C_SIN, alpha=0.12, lw=0)
    ax.axvline(1, color='k', lw=0.8, ls='--')
    ax.text(0.98, 0.02, f"random effects: OR {M['OR_per_zP']:.2f} ({M['OR_per_zP_ci'][0]:.2f}–{M['OR_per_zP_ci'][1]:.1f})\n"
                        f"tau$^2$ = {M['tau2']:.2f}, I$^2$ = {M['I2']:.0f}%, Q = {M['Q']:.2f} (p = {M['Q_p']:.2f})",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5)
    ax.set_xscale('log'); ax.set_xlim(0.7, 1e4); ax.set_yticks([]); ax.set_ylim(-1.1, len(rows) - 0.1)
    ax.spines['left'].set_visible(False); ax.set_xlabel('Odds ratio per unit of z$_P$ (5 mmHg of gradient)')
    save(fig, 'forest')

def centred():
    """Every group on one scale: outcome odds and gradient relative to that cohort's own average."""
    import json as _json
    f = f'{OUT}/meta_slope.json'
    if not os.path.exists(f): indices.main()
    C = _json.load(open(f))['centred']; P = pd.DataFrame(C['points'])
    lab = lambda s: s.replace(' / SFSS_or_dysfunction', ' (SFSS)').replace(' / mortality_or_graft_loss', ' (mortality)')
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    strata_ = sorted(P.stratum.unique()); cols = dict(zip(strata_, ['#1f5fbf', '#c62828', '#2e7d32', '#f28c28', '#7b3fa0']))
    x = np.linspace(P.zP_c.min() * 1.15, P.zP_c.max() * 1.15, 20)
    ax.fill_between(x, np.exp(C['beta_ci'][0] * x), np.exp(C['beta_ci'][1] * x), color='0.85', lw=0, zorder=0)
    ax.plot(x, np.exp(C['beta'] * x), color='k', lw=1.4, zorder=1)
    for st in strata_:
        g = P[P.stratum == st]
        ax.scatter(g.zP_c, np.exp(g.log_odds_c), s=20 + g.n * 0.8, color=cols[st], edgecolor='k', lw=0.4, label=lab(st), zorder=3)
        if len(g) == 2: ax.plot(g.zP_c, np.exp(g.log_odds_c), '-', color=cols[st], lw=0.8, alpha=0.6, zorder=2)
    ax.axhline(1, color='k', lw=0.6, ls=':'); ax.axvline(0, color='k', lw=0.6, ls=':')
    ax.set_yscale('log')
    ax.set_xlabel('z$_P$ relative to the cohort mean'); ax.set_ylabel("Outcome odds relative to the cohort mean")
    ax.text(0.03, 0.97, f"OR {C['OR_per_zP']:.2f} per unit z$_P$ ({C['OR_per_zP_ci'][0]:.2f}–{C['OR_per_zP_ci'][1]:.2f})\n"
                        f"{C['OR_per_mmHg']:.2f} per mmHg · weighted R$^2$ = {C['r2_weighted']:.2f}",
            transform=ax.transAxes, va='top', fontsize=6.5)
    ax.legend(loc='lower right', fontsize=6)
    save(fig, 'centred')

def risk_change():
    """Odds ratio and absolute risk implied by lowering the gradient, for several baseline rates.
    Uses the intercept-free centred slope, the same estimate the calculator applies."""
    import json as _json
    f = f'{OUT}/meta_slope.json'
    if not os.path.exists(f): indices.main()
    C = _json.load(open(f))['centred']
    L = dict(beta=C['beta'], beta_ci=C['beta_ci'], normal_gradient_mmHg=indices.NORMAL_GRADIENT); fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2))
    ax = axs[0]; dg = np.linspace(-12, 4, 100)
    OR = np.exp(L['beta'] * dg / L['normal_gradient_mmHg'])
    lo = np.exp(L['beta_ci'][0] * dg / L['normal_gradient_mmHg']); hi = np.exp(L['beta_ci'][1] * dg / L['normal_gradient_mmHg'])
    ax.fill_between(dg, np.minimum(lo, hi), np.maximum(lo, hi), color=C_SIN, alpha=0.15, lw=0, label='95% bootstrap band')
    ax.plot(dg, OR, color=C_SIN, lw=1.4); ax.axhline(1, color='k', lw=0.6); ax.axvline(0, color='k', lw=0.6)
    ax.set_yscale('log'); ax.set_xlabel('Change in portocaval gradient (mmHg)'); ax.set_ylabel('Odds ratio for the outcome')
    ax.set_ylim(2e-2, 5)
    for g, t in ((-5, '−5 mmHg'), (-10, '−10 mmHg')):
        o = np.exp(L['beta'] * g / L['normal_gradient_mmHg']); ax.plot(g, o, 'o', color=C_SIN); ax.annotate(f'{t}: OR {o:.2f}', (g, o), xytext=(6, -2), textcoords='offset points', fontsize=6.5)
    ax.legend(fontsize=6.5, loc='upper left')
    ax = axs[1]
    for p0, c in ((0.05, '#2e7d32'), (0.10, '#8a9a1e'), (0.20, '#f28c28'), (0.35, '#c62828')):
        p = [100 * indices.risk_after(p0, g / L['normal_gradient_mmHg'], L['beta']) for g in dg]
        ax.plot(dg, p, color=c, lw=1.4, label=f'{int(p0*100)}%')
    ax.axvline(0, color='k', lw=0.6); ax.set_xlabel('Change in portocaval gradient (mmHg)')
    ax.set_ylabel('Outcome risk (%)'); ax.set_ylim(0, 80)
    ax.legend(fontsize=6.5, title="baseline rate at the\nstarting gradient", title_fontsize=6.5, loc='upper left')
    fig.tight_layout(); save(fig, 'risk_change')

def _series(col, name, xlabel, marker_text):
    d0 = series().dropna(subset=[col]); d = d0[d0.in_main_analysis]; x = d0[~d0.in_main_analysis]
    a = d[d.outcome_type == 'SFSS_or_dysfunction']; m = d[d.outcome_type != 'SFSS_or_dysfunction']
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    for st, g in d.groupby('study'):
        if len(g) == 2 and g.outcome.nunique() == 1: ax.plot(g[col], g.pct, '-', color='0.65', lw=0.7, zorder=1)
    ax.scatter(a[col], a.pct, s=12 + a.n.clip(upper=160) / 2.2, color=C_SIN, edgecolor='k', lw=0.4, zorder=3, label='SFSS / graft dysfunction')
    ax.scatter(m[col], m.pct, s=12 + m.n.clip(upper=160) / 2.2, color=C_IN, marker='s', edgecolor='k', lw=0.4, zorder=3, label='mortality / graft loss')
    try:
        from adjustText import adjust_text
        texts = [ax.text(r[col], r.pct, f'{r.study}: {r.group}'[:38], fontsize=5) for _, r in d0.iterrows()]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.3, color='0.5'), expand=(1.6, 2.2), lim=2000)
    except ImportError:
        for _, r in d0.iterrows(): ax.annotate(r.study, (r[col], r.pct), xytext=(4, 2), textcoords='offset points', fontsize=5)
    rho, p = spearmanr(d[col], d.pct); txt = f'Spearman ρ = {rho:.2f}, p = {p:.3f}\n{len(d)} groups'
    if len(x):
        ax.scatter(x[col], x.pct, s=12 + x.n.clip(upper=160) / 2.2, facecolor='none', edgecolor='0.4', lw=0.8, zorder=3, label='defined by a cut-off (not in ρ)')
        rho2, p2 = spearmanr(d0[col], d0.pct); txt += f'\nwith cut-off groups: ρ = {rho2:.2f}, p = {p2:.3f}'
    ax.text(0.98, 0.80, txt, transform=ax.transAxes, ha='right', va='top', fontsize=7)
    ax.axvline(2, color='k', ls='--', lw=0.6); ax.text(2.04, 60, marker_text, fontsize=6)
    if col == 'zP': ax.axvspan(1, 2, color='0.9', lw=0)
    ax.set_xlabel(xlabel); ax.set_ylabel('Outcome (%)'); ax.set_ylim(-3, 108); ax.legend(loc='upper left', fontsize=6.5); save(fig, name)
def series_pressure(): _series('zP', 'series_pressure', 'z$_P$ = (PVP − CVP)/5', 'PVP 15 mmHg\nat CVP 5')
def series_flow(): _series('zF', 'series_flow', 'z$_F$ = graft PVF per 100 g / donor PVF per 100 g', '250 mL/min/100 g\n(donor ref. 120)')

def plane():
    d = series(); fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.fill_between([0, 5], [0, 0], [1, 1], color='0.92', lw=0); ax.fill_betweenx([0, 5], 0, 1, color='0.92', lw=0)
    ax.add_patch(Rectangle((1, 1), 1, 1, color=C_OUT, alpha=0.18, lw=0)); ax.text(1.5, 1.22, 'target\nwindow', ha='center', va='center', fontsize=6.5)
    ax.plot([0, 5], [0, 5], 'k--', lw=0.7); ax.text(4.35, 4.7, 'donor outflow resistance', fontsize=6.5, rotation=45, ha='center', va='center')
    ax.plot(1, 1, 'ko', ms=5, zorder=5); ax.text(1.12, 0.62, 'donor liver', fontsize=7)
    colr = lambda p: '#c62828' if p > 20 else '#f28c28' if p > 5 else '#2e7d32'
    for _, r in d[d.zF.notna() & d.zP.isna() & (d.outcome_type == 'SFSS_or_dysfunction') & d.in_main_analysis].iterrows(): ax.plot(r.zF, r.zF, 'o', color=colr(r.pct), ms=3.5 + r.n ** 0.5 / 2, mec='k', mew=0.3, alpha=0.9, zorder=4)
    for _, r in d.dropna(subset=['zF', 'zP']).iterrows(): ax.plot(r.zF, r.zP, 's', color=colr(r.pct), ms=3.5 + r.n ** 0.5 / 2, mec='k', mew=0.3, alpha=0.9, zorder=4)
    ax.annotate('grafts with reconstructed\noutflow', xy=(3.8, 1.5), xytext=(2.4, 0.3), fontsize=6.5, arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4'))
    ax.annotate('small grafts, standard outflow\n(flow-only series on the diagonal)', xy=(2.95, 3.05), xytext=(1.1, 4.1), fontsize=6.5, arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4'))
    ax.add_patch(Rectangle((0.4, 2), 0.6, 2, color='#7b3fa0', alpha=0.18, lw=0)); ax.text(0.7, 3.0, 'cirrhosis\nHVPG 10–20', ha='center', va='center', fontsize=6.5)
    ax.annotate('', xy=(2.45, 2.45), xytext=(1.07, 1.07), arrowprops=dict(arrowstyle='->', color='0.3', lw=1)); ax.text(1.55, 2.25, 'partial graft or\nresection: × h/g', fontsize=6.5, rotation=45)
    ax.set_xlim(0, 5); ax.set_ylim(0, 5); ax.set_xlabel('z$_F$ (portal flow per lobule / donor)'); ax.set_ylabel('z$_P$ (sinusoidal pressure / normal)'); ax.set_aspect('equal'); save(fig, 'plane')

def resection():
    fig, ax = plt.subplots(figsize=(3.8, 3.2)); g = np.linspace(0.3, 1, 100)
    for hv, lab, c in ((5, 'HVPG 5 mmHg (normal)', '#2e7d32'), (8, 'HVPG 8 mmHg', '#8a9a1e'), (10, 'HVPG 10 mmHg (CSPH)', '#f28c28'), (12, 'HVPG 12 mmHg', '#c62828')): ax.plot(100 * (1 - g), hv / 5 / g, color=c, lw=1.3, label=lab)
    ax.axhspan(1, 2, color='0.9', lw=0); ax.axhline(2, color='k', ls='--', lw=0.7); ax.axvline(60, color='0.5', lw=0.6); ax.text(59, 4.55, 'right hepatectomy\n(≈ 60%)', fontsize=6.5, ha='right')
    ax.set_xlabel('Parenchyma removed (%)'); ax.set_ylabel('z$_P$ of the remnant after resection'); ax.set_ylim(0, 5); ax.set_xlim(0, 70); ax.legend(fontsize=6.5, loc='upper left'); save(fig, 'resection')

def load_curve():
    fig, ax = plt.subplots(figsize=(3.8, 3.2)); gf = np.linspace(0.25, 1.2, 100)
    for h, c in ((1.0, C_IN), (1.5, '#f28c28'), (2.0, C_SIN)): ax.plot(gf / 0.5, h / gf, color=c, lw=1.3, label=f'recipient inflow {h:g}× donor')
    ax.axhspan(1, 2, color='0.9', lw=0); ax.axhline(2, color='k', ls='--', lw=0.6); ax.axhline(1, color='k', lw=0.6); ax.text(0.53, 1.5, 'target\nwindow', fontsize=7, va='center')
    ax.set_xlabel('Graft size (≈ GRWR, %)'); ax.set_ylabel('Load per lobule, z = h/g\n(at donor outflow resistance)'); ax.set_ylim(0, 7); ax.set_xlim(0.5, 2.4); ax.legend(fontsize=7); save(fig, 'load_curve')

def interventions():
    iv = pd.read_csv(f'{DATA}/interventions.tsv', sep='\t')
    items = [(f"{r.manoeuvre} ({r.study}, n={r.n})", r.delta_PVP_mmHg, f"PVF {r.delta_PVF_pct:+.0f}% ({r.flow_source})" if pd.notna(r.delta_PVF_pct) else '') for _, r in iv.iterrows()]
    fig, ax = plt.subplots(figsize=(5.2, 2.8)); y = np.arange(len(items))[::-1]
    for yi, (n, dp, note) in zip(y, items):
        if not np.isnan(dp): ax.barh(yi, dp, 0.55, color=C_SIN, alpha=0.85)
        ax.text(0.4, yi, n, fontsize=6.5, va='center'); 
        if note: ax.text((dp if not np.isnan(dp) else 0) - 0.3, yi, note, fontsize=6.5, ha='right', va='center', color=C_IN)
    ax.set_yticks([]); ax.axvline(0, color='k', lw=0.6); ax.spines['left'].set_visible(False); ax.set_xlim(-11, 13); ax.set_ylim(-0.7, len(items) - 0.3); ax.set_xticks([-10, -5, 0]); ax.set_xlabel('Change in final PVP (mmHg)'); save(fig, 'interventions')

# ------------------------------------------------------------------ model plots
def network_2d(R=4.5, b=0.75):
    from model import network as lc
    G = lc.build_liver_graph(R=R); gin, gout = lc.split_liver(G)
    tree = min((lc.best_sinusoid_tree(gin, b, seed=s) for s in range(3)), key=lambda s: s['D'])
    fig, ax = plt.subplots(figsize=(4.5, 4.5)); xy = G['xy']; cen, cor, pairs, Lh = lc.hex_lobules(R, 1.0)
    for (x, y) in cen: th = np.pi / 6 + np.arange(7) * np.pi / 3; ax.plot(x + Lh * np.cos(th), y + Lh * np.sin(th), color='0.82', lw=0.4, zorder=0)
    w, keep = tree['w'], tree['keep']; r = np.zeros(gin.m); r[keep] = (w[keep] * gin.L[keep]) ** 0.25; rm = r[keep].max()
    for k in keep: i, j = gin.E[k]; ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]], color=C_SIN if gin.kind[k] == 'sin' else C_IN, lw=0.3 + 3.2 * r[k] / rm, solid_capstyle='round', zorder=2)
    outl = lc.tree_allocation(gout, np.arange(gout.m), b); ro = (outl['w'] * gout.L) ** 0.25
    for k in range(gout.m): i, j = gout.E[k]; ax.plot([xy[i, 0], xy[j, 0]], [xy[i, 1], xy[j, 1]], color=C_OUT, lw=0.2 + 2.2 * ro[k] / ro.max(), alpha=0.45, solid_capstyle='round', zorder=1)
    ax.plot(*xy[gin.src], 'ko', ms=5, zorder=5); ax.set_aspect('equal'); ax.axis('off')
    ax.legend(handles=[Line2D([], [], color=C_IN, lw=2, label='portal tree'), Line2D([], [], color=C_SIN, lw=2, label='sinusoids (one per lobule)'), Line2D([], [], color=C_OUT, lw=2, alpha=0.5, label='hepatic venous tree'), Line2D([], [], marker='o', color='k', lw=0, label='hilum')], loc='lower center', fontsize=6.5, ncol=2, bbox_to_anchor=(0.5, -0.08)); save(fig, 'network_2d')

def network_3d(R=3.4, b=0.75):
    from model import network as lc, network3d as l3
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    G3 = l3.build_liver3d(R=R); g3in, g3out = lc.split_liver(G3); ft = l3.FastTree(g3in, b)
    best = min((ft.solve(seed=s) for s in range(3)), key=lambda s: s['D']); sol = lc.tree_allocation(g3in, best['keep'], b)
    r3 = np.zeros(g3in.m); r3[sol['keep']] = (sol['w'][sol['keep']] * g3in.L[sol['keep']]) ** 0.25; X = G3['xy']
    segs, cols, lws = [], [], []
    for k in sol['keep']: i, j = g3in.E[k]; segs.append([X[i], X[j]]); cols.append(C_SIN if g3in.kind[k] == 'sin' else C_IN); lws.append(0.35 if g3in.kind[k] == 'sin' else 0.3 + 2.6 * r3[k] / r3.max())
    o3 = lc.tree_allocation(g3out, np.arange(g3out.m), b); ro3 = (o3['w'] * g3out.L) ** 0.25
    for k in range(g3out.m): i, j = g3out.E[k]; segs.append([X[i], X[j]]); cols.append('#7fbf7f'); lws.append(0.15 + 1.4 * ro3[k] / ro3.max())
    fig = plt.figure(figsize=(5, 4)); ax = fig.add_subplot(111, projection='3d'); ax.add_collection3d(Line3DCollection(segs, colors=cols, linewidths=lws, alpha=0.85))
    ax.scatter([0], [0], [0], color='k', s=18, depthshade=False)
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, np.pi / 2, 30)); Rh = R + 0.15
    ax.plot_surface(Rh * np.cos(u) * np.sin(v), Rh * np.sin(u) * np.sin(v), Rh * np.cos(v), color='0.55', alpha=0.10, linewidth=0, shade=False)
    ax.set_xlim(-Rh, Rh); ax.set_ylim(-Rh, Rh); ax.set_zlim(0, Rh); ax.view_init(elev=24, azim=-58); ax.set_box_aspect((1, 1, 0.55), zoom=1.2); ax.set_axis_off()
    ax.set_title(f'Hemispheric liver, {G3["n_lob"]} lobules, b = {b:g}', fontsize=8); save(fig, 'network_3d')

def scaling():
    g3 = pd.read_csv(f'{DATA}/graph3d_all.tsv', sep='\t'); g2 = pd.read_csv(f'{DATA}/graph_scaling_all.tsv', sep='\t'); from model import scaling as sc
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2)); ax = axs[0]
    for b, gr in g2.groupby('b'): Dt = gr.D_in_perF2 + gr.D_out / gr.N ** 2; e = np.polyfit(np.log(gr.N), np.log(Dt), 1)[0]; ax.loglog(gr.N, Dt / Dt.iloc[0], 'o-', label=f'b = {b:g}, exponent {e:.2f}')
    ax.set_xlabel('Number of lobules N (2D disc)'); ax.set_ylabel('D*/F² (normalised)'); ax.legend(fontsize=6.5); ax.set_title('2D graph', fontsize=8)
    ax = axs[1]
    for b, gr in g3.groupby('b'):
        eg = np.polyfit(np.log(gr.N.values[-3:]), np.log(gr.D_tot_perF2.values[-3:]), 1)[0]
        ana = [sc.canopy_optimum(g, b, 3) for g in range(4, 9)]; Na = np.array([a['N'] for a in ana]); Da = np.array([a['D'] for a in ana]); ec = np.polyfit(np.log(Na[-3:]), np.log(Da[-3:]), 1)[0]
        ax.loglog(gr.N, gr.D_tot_perF2 / gr.D_tot_perF2.iloc[0], 's-', label=f'b = {b:.3g}: graph {eg:.2f}, closed form {ec:.2f}')
    ax.set_xlabel('Number of lobules N (3D hemisphere)'); ax.set_ylabel('D*/F² (normalised)'); ax.legend(fontsize=6.5); ax.set_title('3D graph vs closed form', fontsize=8); fig.tight_layout(); save(fig, 'scaling')

def allometry():
    from model import scaling as sc
    gs = np.arange(8, 19); rows = []
    for f, flab in ((0.77, 'portal flow ~ M^0.77'), (0.88, 'total hepatic flow ~ M^0.88')):
        for l, n in ((0.0, 0.89), (0.10, 0.59)):
            for b in (1.0, 0.9, 0.85, 0.8, 0.75, 2 / 3, 0.5):
                Rr = [sc.canopy_optimum(g, b, 3) for g in gs]; N = np.array([r['N'] for r in Rr])
                eD = np.polyfit(np.log(N), np.log([r['D'] for r in Rr]), 1)[0]; eV = np.polyfit(np.log(N), np.log([r['V'] for r in Rr]), 1)[0]
                c = b * (f + eD * n + 3 * l) / 2; rows.append(dict(flow=flab, l=l, b=b, tau0_exp=(2 * f - (2 + b) * c / b + eD * n + 3 * l) / 2, V_exp=c / b + eV * n))
    cl = pd.DataFrame(rows); os.makedirs(OUT, exist_ok=True); cl.to_csv(f'{OUT}/allometric_closure.tsv', sep='\t', index=False, float_format='%.4f')
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2))
    for (fl, l), gr in cl.groupby(['flow', 'l']):
        axs[0].plot(gr.b, gr.tau0_exp, '-' if l > 0 else ':', marker='o', ms=3, label=f'{fl}, lobule ~ M^{l:.2f}'); axs[1].plot(gr.b, gr.V_exp, '-' if l > 0 else ':', marker='s', ms=3, label=f'{fl}, lobule ~ M^{l:.2f}')
    axs[0].axhline(0, color='k', lw=0.6); axs[0].axhspan(-0.05, 0.05, color='0.9', lw=0); axs[0].set_xlabel('Maintenance exponent b'); axs[0].set_ylabel('Exponent of the shear set-point with body mass'); axs[0].legend(fontsize=6)
    axs[1].axhspan(0.84, 0.88, color=C_OUT, alpha=0.25, lw=0, label='hepatic blood volume ~ M^0.86'); axs[1].set_xlabel('Maintenance exponent b'); axs[1].set_ylabel('Exponent of vascular mass with body mass'); axs[1].legend(fontsize=6)
    fig.tight_layout(); save(fig, 'allometry')

def shear_profile():
    fig, ax = plt.subplots(figsize=(4, 3.2))
    for name, (rd, rl, n, g) in {'hepatic artery (0.74, n 2.76)': (0.74, 0.66, 2.76, 10), 'portal vein (0.70, n 2.80)': (0.70, 0.72, 2.80, 10), 'hepatic vein (0.59, n 3.22)': (0.59, 0.64, 3.22, 8)}.items():
        per = (1 / n) / rd ** 3; i = np.arange(g + 1); ax.semilogy(i, per ** i, 'o-', ms=3, label=f'{name}: ×{per**g:.2f}')
    ax.axhline(1, color='k', lw=0.6); ax.set_xlabel('Branching generation'); ax.set_ylabel('Wall shear relative to trunk'); ax.legend(fontsize=6.5); ax.set_title('Measured diameter ratios (Debbaut), equal flow split', fontsize=8); save(fig, 'shear_profile')

def certificate():
    ce = pd.read_csv(f'{DATA}/certificate_small.tsv', sep='\t'); fig, ax = plt.subplots(figsize=(3.6, 3.0))
    ax.bar(ce.b.astype(str), ce.frac_seeds_at_opt * 100, color=C_IN, width=0.5)
    for i, r in ce.iterrows(): ax.text(i, r.frac_seeds_at_opt * 100 + 2, f'gap {abs(r.rel_gap):.0e}', ha='center', fontsize=6.5)
    ax.set_ylim(0, 115); ax.set_xlabel('Maintenance exponent b'); ax.set_ylabel('Random starts reaching the exact optimum (%)'); ax.set_title('7 lobules, 279 936 trees enumerated', fontsize=8); save(fig, 'certificate')

PLOTS = dict(nomogram=nomogram, within_study=within_study, forest=forest, centred=centred, risk_change=risk_change, series_pressure=series_pressure, series_flow=series_flow, plane=plane, resection=resection, load_curve=load_curve,
             interventions=interventions, network_2d=network_2d, network_3d=network_3d, scaling=scaling, allometry=allometry, shear_profile=shear_profile, certificate=certificate)

if __name__ == '__main__':
    names = sys.argv[1:]
    if not names: print(__doc__); sys.exit(0)
    for n in (PLOTS if names == ['all'] else names): PLOTS[n]()
