"""Build sfss_calculator.html (repository root) from results/logit_zP.json.

    python src/build_app.py

The coefficients, the bootstrap band, the fitted z_P range, the normal gradient, the default CVP and the
provenance counts are injected from the analysis outputs, so a re-run with new data updates the app.
"""
import json, os
ROOT = os.path.join(os.path.dirname(__file__), '..')
L = json.load(open(f'{ROOT}/results/logit_zP.json'))
groups = [dict(s=f"{g['study']} {g['group'][:22]}", z=round(g['zP'], 2), p=round(100 * g['events'] / g['n'], 1), n=g['n']) for g in L['groups_used']]
studies = ', '.join(sorted({g['study'] for g in L['groups_used']}))
assumed = sum(1 for g in L['groups_used'] if g['CVP_basis'] != 'reported')
t = open(f'{ROOT}/src/app_template.html').read()
rep = {'__CURVE__': json.dumps({k: L[k] for k in ('z', 'p', 'lo', 'hi')}), '__GROUPS__': json.dumps(groups),
       '__B0__': f"{L['b0']:.4f}", '__B1__': f"{L['b1']:.4f}", '__B0S__': f"{L['b0']:.2f}".replace('-', '−'), '__B1S__': f"{L['b1']:.2f}",
       '__GN__': f"{L['normal_gradient_mmHg']:g}", '__CVP__': f"{L['default_cvp_mmHg']:g}", '__ZMIN__': f"{L['zP_min']:.3f}", '__ZMAX__': f"{L['zP_max']:.3f}",
       '__PATIENTS__': str(L['patients']), '__EVENTS__': str(L['events']), '__STUDIES__': studies, '__ASSUMED__': str(assumed), '__NGROUPS__': str(L['groups'])}
for k, v in rep.items(): t = t.replace(k, v)
assert '__' not in t.split('<script>')[0].replace('__proto__', ''), 'unfilled placeholder'
open(f'{ROOT}/sfss_calculator.html', 'w').write(t); print('sfss_calculator.html written')
