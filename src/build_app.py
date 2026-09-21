"""Build sfss_calculator.html (repository root) from results/within_study_fit.json.

    python src/build_app.py

The within-study slope and its bootstrap interval, the study intercepts and groups, the normal gradient and
the default CVP are injected from the analysis outputs, so a re-run with new data updates the calculator.
The calculator reports the effect of CHANGING the gradient (an odds ratio) and converts it to an absolute risk
only with a baseline rate supplied by the user; the baseline level is not transferable between centres.
"""
import json, os
ROOT = os.path.join(os.path.dirname(__file__), '..')
L = json.load(open(f'{ROOT}/results/within_study_fit.json'))
M = json.load(open(f'{ROOT}/results/meta_slope.json')); C = M['centred']
lab = lambda st: st.replace(' / SFSS_or_dysfunction', ' (SFSS)').replace(' / mortality_or_graft_loss', ' (mortality)')
groups = [dict(s=lab(g['stratum']), z=round(g['zP'], 2), p=round(g['pct'], 1), n=g['n']) for g in L['study_groups']]
t = open(f'{ROOT}/src/app_template.html').read()
refs = {lab(g): dict(alpha=round(L['alphas'][g], 4), n=0, events=0, zmin=9e9, zmax=-9e9) for g in L['studies']}
for g in L['study_groups']:
    k = lab(g['stratum']); r = refs[k]; r['n'] += g['n']; r['events'] += g['events']
    r['zmin'] = min(r['zmin'], g['zP']); r['zmax'] = max(r['zmax'], g['zP'])
rep = {'__GROUPS__': json.dumps(groups), '__STUDIES__': json.dumps([lab(s) for s in L['studies']]), '__ALPHAS__': json.dumps({lab(k): round(v, 4) for k, v in L['alphas'].items()}),
       '__BETA__': f"{C['beta']:.4f}", '__BLO__': f"{C['beta_ci'][0]:.4f}", '__BHI__': f"{C['beta_ci'][1]:.4f}",
       '__R2__': f"{C['r2_weighted']:.2f}", '__OR_STRAT__': f"{L['OR_per_zP']:.2f}", '__OR_META__': f"{M['meta']['OR_per_zP']:.2f}", '__I2__': f"{M['meta']['I2']:.0f}",
       '__GN__': f"{L['normal_gradient_mmHg']:g}", '__CVP__': f"{L['default_cvp_mmHg']:g}",
       '__ZMIN__': f"{L['zP_min']:.3f}", '__ZMAX__': f"{L['zP_max']:.3f}",
       '__OR_ZP__': f"{C['OR_per_zP']:.2f}", '__OR_LO__': f"{C['OR_per_zP_ci'][0]:.2f}", '__OR_HI__': f"{C['OR_per_zP_ci'][1]:.2f}",
       '__OR_MM__': f"{C['OR_per_mmHg']:.2f}", '__NGROUPS__': str(L['groups']), '__REFS__': json.dumps(refs), '__REFOPTS__': ''.join(f'<option value="{k}">{k} — {100*r["events"]/r["n"]:.0f}% overall</option>' for k, r in refs.items()), '__NSTUDIES__': str(len({s.split(' / ')[0] for s in L['studies']})),
       '__EVENTS__': str(L['events']), '__PATIENTS__': str(L['patients'])}
for k, v in rep.items(): t = t.replace(k, v)
assert '__' not in t.replace('__proto__', ''), 'unfilled placeholder: ' + t[t.index('__') - 40:t.index('__') + 40]
open(f'{ROOT}/sfss_calculator.html', 'w').write(t); print('sfss_calculator.html written')
