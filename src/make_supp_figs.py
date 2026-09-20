"""Supplementary figures S1-S3 (run from repo root with PYTHONPATH=src after run_all.sh copied data/ to results/)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Liberation Sans','font.size':7,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
def panel(ax,l): ax.set_title(l,loc='left',fontweight='bold',fontsize=14,pad=4)
OUT='results'
g3=pd.read_csv(f'{OUT}/graph3d_all.tsv',sep='\t'); ex=pd.read_csv(f'{OUT}/graph3d_exponents.tsv',sep='\t'); g2=pd.read_csv(f'{OUT}/graph_scaling_all.tsv',sep='\t')
fig,axs=plt.subplots(1,2,figsize=(7.2,3.2)); ax=axs[0]; panel(ax,'A')
for b,gr in g2.groupby('b'):
    Dt=gr.D_in_perF2+gr.D_out/gr.N**2; e=np.polyfit(np.log(gr.N),np.log(Dt),1)[0]; ax.loglog(gr.N,Dt/Dt.iloc[0],'o-',label=f'b = {b:g}, exponent {e:.2f}')
ax.set_xlabel('Number of lobules N (2D disc)'); ax.set_ylabel('D*/F² (normalised)'); ax.legend(fontsize=6)
ax=axs[1]; panel(ax,'B')
for b,gr in g3.groupby('b'):
    e=ex.iloc[(ex.b-b).abs().argmin()]; ax.loglog(gr.N,gr.D_tot_perF2/gr.D_tot_perF2.iloc[0],'s-',label=f'b = {b:.3g}: graph {e.e_graph3D_last3:.2f}, closed form {e.e_sym3D_sameN:.2f}')
ax.set_xlabel('Number of lobules N (3D hemisphere)'); ax.set_ylabel('D*/F² (normalised)'); ax.legend(fontsize=6)
plt.tight_layout(); plt.savefig(f'{OUT}/JHEP_FigS1.png',dpi=300); plt.savefig(f'{OUT}/JHEP_FigS1.pdf'); plt.close()
cl=pd.read_csv(f'{OUT}/closure_kruepunga.tsv',sep='\t'); lab={'porta 0.77 (no-PAH)':'portal flow ~ M^0.77','hepatico total 0.88':'total hepatic flow ~ M^0.88'}
fig,axs=plt.subplots(1,3,figsize=(7.2,3.0)); ax=axs[0]; panel(ax,'A')
for n,m,d in [('mouse',0.025,210),('rat',0.3,330),('human',70,385),('pig',60,620)]: ax.plot(m,d,'o',color='#1f5fbf'); ax.annotate(n,(m,d),xytext=(4,4),textcoords='offset points',fontsize=6)
mm=np.logspace(-2,3,50); ax.plot(mm,330*(mm/0.3)**0.10,'k--',lw=0.8,label='exponent 0.10 [S6]'); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(100,1000); ax.set_xlabel('Body mass (kg)'); ax.set_ylabel('Portocentral distance (µm)'); ax.legend(fontsize=6)
ax=axs[1]; panel(ax,'B')
for (fl,l),gr in cl.groupby(['flujo','l_lobulo']): ax.plot(gr.b,gr.tau0_exp,'-' if l>0 else ':',marker='o',ms=3,label=f'{lab[fl]}, lobule ~ M^{l:.2f}')
ax.axhline(0,color='k',lw=0.6); ax.axhspan(-0.05,0.05,color='0.9',lw=0); ax.set_xlabel('Maintenance exponent b'); ax.set_ylabel('Exponent of the shear set-point with body mass'); ax.legend(fontsize=5)
ax=axs[2]; panel(ax,'C')
for (fl,l),gr in cl.groupby(['flujo','l_lobulo']): ax.plot(gr.b,gr.V_exp,'-' if l>0 else ':',marker='s',ms=3,label=f'{lab[fl]}, lobule ~ M^{l:.2f}')
ax.axhspan(0.84,0.88,color='#2e7d32',alpha=0.25,lw=0,label='hepatic blood volume ~ M^0.86'); ax.set_xlabel('Maintenance exponent b'); ax.set_ylabel('Exponent of vascular mass with body mass'); ax.legend(fontsize=5)
plt.tight_layout(); plt.savefig(f'{OUT}/JHEP_FigS2.png',dpi=300); plt.savefig(f'{OUT}/JHEP_FigS2.pdf'); plt.close()
ce=pd.read_csv(f'{OUT}/certificate_small.tsv',sep='\t')
fig,axs=plt.subplots(1,2,figsize=(7.2,3.0)); ax=axs[0]; panel(ax,'A'); ax.bar(ce.b.astype(str),ce.frac_seeds_at_opt*100,color='#1f5fbf',width=0.5)
for i,r in ce.iterrows(): ax.text(i,r.frac_seeds_at_opt*100+2,f'gap {abs(r.rel_gap):.0e}',ha='center',fontsize=6)
ax.set_ylim(0,115); ax.set_xlabel('Maintenance exponent b'); ax.set_ylabel('Random starts reaching the exact optimum (%)')
ax=axs[1]; panel(ax,'B')
for name,(rd,rl,n,g) in {'hepatic artery':(0.74,0.66,2.76,10),'portal vein':(0.70,0.72,2.80,10),'hepatic vein':(0.59,0.64,3.22,8)}.items():
    per=(1/n)/rd**3; i=np.arange(g+1); ax.semilogy(i,per**i,'o-',ms=3,label=f'{name}: ×{per**g:.2f} over {g} generations')
ax.axhline(1,color='k',lw=0.6); ax.set_xlabel('Branching generation'); ax.set_ylabel('Wall shear relative to trunk'); ax.legend(fontsize=6)
plt.tight_layout(); plt.savefig(f'{OUT}/JHEP_FigS3.png',dpi=300); plt.savefig(f'{OUT}/JHEP_FigS3.pdf'); plt.close(); print('supplementary figures done')
