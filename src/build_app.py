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
groups = [dict(s=g['study'], z=round(g['zP'], 2), p=round(g['pct'], 1), n=g['n']) for g in L['study_groups']]
t = open(f'{ROOT}/src/app_template.html').read()
rep = {'__GROUPS__': json.dumps(groups), '__STUDIES__': json.dumps(L['studies']), '__ALPHAS__': json.dumps({k: round(v, 4) for k, v in L['alphas'].items()}),
       '__BETA__': f"{L['beta']:.4f}", '__BLO__': f"{L['beta_ci'][0]:.4f}", '__BHI__': f"{L['beta_ci'][1]:.4f}",
       '__GN__': f"{L['normal_gradient_mmHg']:g}", '__CVP__': f"{L['default_cvp_mmHg']:g}",
       '__ZMIN__': f"{L['zP_min']:.3f}", '__ZMAX__': f"{L['zP_max']:.3f}",
       '__OR_ZP__': f"{L['OR_per_zP']:.2f}", '__OR_LO__': f"{L['OR_per_zP_ci'][0]:.2f}", '__OR_HI__': f"{L['OR_per_zP_ci'][1]:.2f}",
       '__OR_MM__': f"{L['OR_per_mmHg']:.2f}", '__NGROUPS__': str(L['groups']), '__NSTUDIES__': str(len(L['studies'])),
       '__EVENTS__': str(L['events']), '__PATIENTS__': str(L['patients'])}
for k, v in rep.items(): t = t.replace(k, v)
assert '__' not in t.replace('__proto__', ''), 'unfilled placeholder: ' + t[t.index('__') - 40:t.index('__') + 40]
open(f'{ROOT}/sfss_calculator.html', 'w').write(t); print('sfss_calculator.html written')
