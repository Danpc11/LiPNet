import numpy as np, pandas as pd, glob, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import scaling as sc
df = pd.concat([pd.read_csv(f, sep='\t') for f in glob.glob('results/graph_scaling_R*.tsv')]).sort_values(['b','N'])
df.to_csv('results/graph_scaling_all.tsv', sep='\t', index=False)
fig, axs = plt.subplots(1,3, figsize=(13,4))
ax=axs[0]
rows=[]
for b,gr in df.groupby('b'):
    Dtot = gr.D_in_perF2 + gr.D_out/gr.N**2
    e = np.polyfit(np.log(gr.N), np.log(Dtot),1)[0]
    ana = [sc.canopy_optimum(g,b,2) for g in range(3,8)]
    Na=np.array([a['N'] for a in ana]); Da=np.array([a['D'] for a in ana])
    ea = np.polyfit(np.log(Na[-3:]), np.log(Da[-3:]),1)[0]
    ax.loglog(gr.N, Dtot/Dtot.iloc[0], 'o-', label=f'b={b:g}: grafo {e:.2f}, simétrico 2D {ea:.2f}')
    rows.append(dict(b=b, e_graph_2D=e, e_symmetric_2D=ea))
ax.set_xlabel('N lobulillos'); ax.set_ylabel(r'$D^*/F^2$ (normalizado, C=1)'); ax.set_title('a  exponente de D* con N: grafo explícito vs canopy simétrico', fontsize=9); ax.legend(fontsize=7)
pd.DataFrame(rows).to_csv('results/graph_vs_symmetric_exponents.tsv', sep='\t', index=False, float_format='%.3f')
ax=axs[1]
al = pd.read_csv('results/scaling_allometry.tsv', sep='\t'); al3 = al[al.d==3]
ax.plot(al3.b, al3.gamma_F_vs_N, 'o-', label=r'$\gamma$: $F\propto N^{\gamma}$ (V∝N, $\tau_0$ fijo)')
ax.axhline(0.75/0.87, color='k', ls='--', label='empírico 0.75/0.87 = 0.86')
ax.set_xlabel('b'); ax.set_ylabel(r'$\gamma$'); ax.set_title('b  3D: flujo de entrada vs tamaño a set‑point fijo', fontsize=9); ax.legend(fontsize=7)
ax=axs[2]
ax.plot(al3.b, al3.liver_exp_M_for_invariant_tau0, 's-', color='tab:red', label=r'exponente hígado∝$M^{y}$ predicho ($F\propto M^{3/4}$, $\tau_0$ fijo)')
ax.axhspan(0.85, 0.90, color='0.85', label='empírico 0.85–0.90 (Stahl)')
ax2=ax.twinx(); ax2.plot(al3.b, al3.tau0_exp_M_given_Kleiber_and_0p87, 'd:', color='tab:blue', label=r'exponente de $\tau_0$ con M si hígado∝$M^{0.87}$'); ax2.axhline(0,color='tab:blue',lw=0.5); ax2.set_ylabel(r'$\tau_0 \propto M^{x}$', color='tab:blue')
ax.set_xlabel('b'); ax.set_ylabel('y'); ax.set_title('c  3D: cierre alométrico', fontsize=9)
h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels(); ax.legend(h1+h2,l1+l2,fontsize=7, loc='upper left')
plt.tight_layout(); plt.savefig('results/fig3_scaling.png', dpi=170)
print(pd.DataFrame(rows).round(3).to_string(index=False))
