"""Build sfss_calculator.html (repository root) from results/logit_zP.json and data/series.tsv.

    python src/build_app.py

The app is a single self-contained HTML file; the fitted coefficients, the bootstrap band and the
published groups are injected from the analysis outputs so that a re-run with new data updates the app.
"""
import json, os, pandas as pd
ROOT = os.path.join(os.path.dirname(__file__), '..')
L = json.load(open(f'{ROOT}/results/logit_zP.json'))
groups = [dict(s=f"{g['study']} {g['group'][:22]}", z=round(g['zP'], 2), p=round(100 * g['events'] / g['n'], 1), n=g['n']) for g in L['groups_used']]
studies = ', '.join(sorted({g['study'] for g in L['groups_used']}))
t = open(f'{ROOT}/src/app_template.html').read()
t = (t.replace('__CURVE__', json.dumps({k: L[k] for k in ('z', 'p', 'lo', 'hi')})).replace('__GROUPS__', json.dumps(groups))
      .replace('__B0__', f"{L['b0']:.4f}").replace('__B1__', f"{L['b1']:.4f}").replace('__B0S__', f"{L['b0']:.2f}".replace('-', '−')).replace('__B1S__', f"{L['b1']:.2f}")
      .replace('__PATIENTS__', str(L['patients'])).replace('__EVENTS__', str(L['events'])).replace('__STUDIES__', studies))
open(f'{ROOT}/sfss_calculator.html', 'w').write(t); print('sfss_calculator.html written')
