"""Build sfss_calculator.html (repository root) from src/app_template.html and the results.

    python src/build_app.py

Outcome. The odds ratio for a change of gradient uses the pooled slope mu_beta of the primary analysis
(results/bayes_main.json) with its credible interval. A published series used as reference contributes its level at its
own mean load; curve and band come from the same posterior draws, with the pooled slope (alphas_pooled_slope,
pred_band_pooled_slope). The descriptive view without a reference shows the groups centred within their cohort
(results/meta_slope.json).

Resistance load. The page computes, for the flow load and CVP entered, the zR a rigid network gives a partial graft
(model/liver.two_inflows over the networks of model/compartments.geometries) and the range of the three passive
relations of model/distension.py, with the same constants as the Python code.
"""
import json, os, sys
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'src'))
from model import distension as ds, compartments as cm

B = json.load(open(f'{ROOT}/results/bayes_main.json'))['primary (SFSS or early dysfunction)']
L = json.load(open(f'{ROOT}/results/within_study_fit.json'))
C = json.load(open(f'{ROOT}/results/meta_slope.json'))['centred']
lab = lambda st: st.replace(' / SFSS_or_dysfunction', ' (SFSS)').replace(' / mortality_or_graft_loss', ' (mortality)')

refs = {lab(g): dict(alpha=round(B['alphas_pooled_slope'][g], 4), band=B['pred_band_pooled_slope'][g], n=B['group_counts'][g]['n'],
                     events=B['group_counts'][g]['events'], zmin=B['group_counts'][g]['zmin'], zmax=B['group_counts'][g]['zmax']) for g in B['strata']}
centred = [dict(s=lab(p['stratum']), zc=round(p['zP_c'], 4), lc=round(p['log_odds_c'], 4), n=p['n']) for p in C['points']]
net = dict(geoms=[[float(w[:k].sum()), float(rem)] for _, w, k, rem in cm.geometries()], a=cm.A_RATIO, g0=5.0, cube=ds.CUBE,
           doubling=list(ds.DOUBLING_MMHG), f=list(ds.FLOW_FIXED_SHARE), cvps=list(ds.CVPS), measured=max(ds.MEASURED_UP_TO.values()))

rep = {'__STUDIES__': json.dumps([lab(s) for s in B['strata']]), '__REFS__': json.dumps(refs), '__ZGRID__': json.dumps(B['zgrid']),
       '__CENTRED__': json.dumps(centred), '__NET__': json.dumps(net),
       '__BETA__': f"{B['mu_beta']:.4f}", '__BLO__': f"{B['mu_beta_ci'][0]:.4f}", '__BHI__': f"{B['mu_beta_ci'][1]:.4f}",
       '__GN__': f"{L['normal_gradient_mmHg']:g}", '__CVP__': f"{L['default_cvp_mmHg']:g}",
       '__OR_ZP__': f"{B['OR_per_zP']:.2f}", '__OR_LO__': f"{B['OR_per_zP_ci'][0]:.2f}", '__OR_HI__': f"{B['OR_per_zP_ci'][1]:.2f}", '__OR_MM__': f"{B['OR_per_mmHg']:.2f}",
       '__REFOPTS__': ''.join(f'<option value="{k}">{k} — {100 * r["events"] / r["n"]:.0f}% overall</option>' for k, r in refs.items()),
       '__NSTUDIES__': str(len(B['strata'])), '__NGROUPS__': str(B['groups'])}
t = open(f'{ROOT}/src/app_template.html').read()
for k, v in rep.items(): t = t.replace(k, v)
assert '__' not in t.replace('__proto__', ''), 'unfilled placeholder: ' + t[t.index('__') - 40:t.index('__') + 40]
open(f'{ROOT}/sfss_calculator.html', 'w').write(t); print('sfss_calculator.html written')
