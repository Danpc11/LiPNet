import numpy as np, pandas as pd, pickle, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, liver_canopy as lc, warnings; warnings.filterwarnings('ignore')
OUT='results'
G = lc.build_liver_graph(R=4.5); gin, gout = lc.split_liver(G)
st = pickle.load(open(f'{OUT}/partB_states.pkl','rb'))
dfB = pd.read_csv(f'{OUT}/partB_campaign.tsv', sep='\t'); dfA = pd.read_csv(f'{OUT}/partA_infer_b.tsv', sep='\t')
dfE = pd.read_csv(f'{OUT}/partA_efficiency.tsv', sep='\t')
xy = G['xy']
def draw(ax, g, w, keep, title, cmap='viridis'):
    r = np.zeros(g.m); r[keep] = (w[keep]*g.L[keep])**0.25; rm = r[keep].max()
    # lobule hexagons
    cen, cor, pairs, Lh = lc.hex_lobules(4.5, 1.0)
    for (x,y) in cen:
        th = np.pi/6 + np.arange(7)*np.pi/3
        ax.plot(x+Lh*np.cos(th), y+Lh*np.sin(th), color='0.85', lw=0.5, zorder=0)
    for k in keep:
        i,j = g.E[k]; c = 'tab:red' if g.kind[k]=='sin' else 'tab:blue'
        ax.plot([xy[i,0],xy[j,0]],[xy[i,1],xy[j,1]], color=c, lw=0.3+4.5*r[k]/rm, alpha=0.9, solid_capstyle='round')
    ax.plot(*xy[g.src], 'ko', ms=6); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title, fontsize=9)
fig, axs = plt.subplots(1,3, figsize=(13,4.4))
b=0.75
k0='0.7500_0.0'; k4='0.7500_4.0'
draw(axs[0], gin, st[k0]['lor_w'], np.arange(gin.m), 'a  Lorente canopy: ratios 3^-1/3, all 6 sinusoids/lobule\n(b = 3/4, allocation not optimal)')
draw(axs[1], gin, st[k0]['tree_w'], st[k0]['tree_keep'], f'b  best tree, steady demand (σ=0), b = 3/4\n{len(st[k0]["tree_keep"])} active edges, 1 sinusoid/lobule')
r=dfB[(dfB.b==0.75)&(dfB.sigma==4.0)].iloc[0]
draw(axs[2], gin, st[k4]['ad_w'], st[k4]['ad_keep'], f'c  local adaptation, fluctuating demand σ=4, b = 3/4\n{int(r.active_adapted)} active, {r.sin_per_lobule:.2f} sinusoids/lobule, β={int(r.beta_adapted)} loops')
plt.tight_layout(); plt.savefig(f'{OUT}/fig1_liver_networks.png', dpi=170); plt.close()

fig, axs = plt.subplots(2,2, figsize=(10,7.6))
ax=axs[0,0]
for b,gr in dfB.groupby('b'):
    ax.plot(gr.sigma, gr.sin_per_lobule, 'o-', label=f'b={b:.3g}')
ax.axhline(1, color='0.5', ls=':'); ax.set_xlabel('demand fluctuation σ'); ax.set_ylabel('active sinusoids per lobule'); ax.legend(fontsize=8); ax.set_title('a  redundancy kept by fluctuating demand', fontsize=10)
ax=axs[0,1]
for b,gr in dfB.groupby('b'):
    ax.plot(gr.sigma, gr.eta_adapted, 'o-', label=f'b={b:.3g}')
ax.axhline(1, color='0.5', ls=':'); ax.set_xlabel('demand fluctuation σ'); ax.set_ylabel(r'$\langle D\rangle_{best\ tree}/\langle D\rangle_{adapted}$'); ax.set_title('b  adapted (looped) network vs best tree at equal budget', fontsize=10); ax.legend(fontsize=8)
ax=axs[1,0]
for name,gr in dfE.groupby('geometry'):
    gr=gr.sort_values('b'); ax.plot(gr.b, gr.eta, 'o-', label=name)
ax.set_xlabel('maintenance exponent b'); ax.set_ylabel(r'$\eta_b = D^*/D$ (own support, own budget)'); ax.set_yscale('log'); ax.legend(fontsize=7); ax.set_title('c  generation model: efficiency of measured geometry', fontsize=10)
ax=axs[1,1]
y=np.arange(len(dfA))
ax.errorbar(dfA.b_hat, y, xerr=[dfA.b_hat-dfA.b_lo, dfA.b_hi-dfA.b_hat], fmt='o', color='k', capsize=3, label='measured n')
ax.errorbar(dfA.b_hat_n3, y+0.15, xerr=[dfA.b_hat_n3-dfA.b_lo_n3, dfA.b_hi_n3-dfA.b_hat_n3], fmt='s', color='tab:gray', capsize=3, label='n = 3')
for v,l,c in ((1,'Murray','tab:blue'),(0.75,'Kleiber','tab:green'),(2/3,'Rubner','tab:red')): ax.axvline(v, color=c, ls='--', lw=1); ax.text(v, len(dfA)-0.55, l, color=c, fontsize=8, ha='center')
ax.set_yticks(y); ax.set_yticklabels(dfA.tree); ax.set_xlabel('inferred maintenance exponent b'); ax.set_title('d  b from Lorente/Debbaut (d, L, n) ratios, 16–84% band', fontsize=10); ax.legend(fontsize=8, loc='lower right'); ax.set_xlim(-0.1, 2.0); ax.set_ylim(-0.5, len(dfA)-0.2)
plt.tight_layout(); plt.savefig(f'{OUT}/fig2_liver_summary.png', dpi=170); plt.close()
print(dfB[['b','sigma','eta_lorente','eta_adapted','excess_adapted_pct','sin_per_lobule','beta_adapted','active_tree','active_adapted','C_in_share','eta_organ_lorente','eta_organ_adapted']].round(3).to_string(index=False))
