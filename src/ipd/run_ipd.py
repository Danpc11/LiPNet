"""Command line for the patient-level layer. Run from the repository root:

    python src/ipd/run_ipd.py patients      # the published patient-level table
    python src/ipd/run_ipd.py ecological    # slope for a patient vs for a group mean
    python src/ipd/run_ipd.py figures
    python src/ipd/run_ipd.py all

Results go to results/ipd_*.json, figures to results/ipd_*.pdf and .png.
"""
import argparse, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np, pandas as pd
import indices
from ipd import ecological as ec
from ipd import patient_level as pl

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
OUT = os.path.join(ROOT, 'results')
SIGMAS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2)


def primary():
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    if 'in_main_analysis' in d:
        d = d[d.in_main_analysis == True]
    return d[d.outcome_type == 'SFSS_or_dysfunction']


def dump(name, obj, quiet=False):
    os.makedirs(OUT, exist_ok=True)
    p = os.path.join(OUT, name + '.json')
    json.dump(obj, open(p, 'w'), indent=1, sort_keys=True)
    if not quiet:
        print(json.dumps(obj, indent=1, sort_keys=True)[:4000])
    print('wrote', os.path.relpath(p, ROOT))


def cmd_ecological(args):
    d = primary()
    out = dict(reported_within_group_sd_mmHg=sorted(
        float(x) for x in d.gradient_sd.dropna().unique()),
        pooled_sigma_zP=ec.pooled_sigma(d))
    out['sweep'] = ec.sweep(d, SIGMAS, chains=args.chains, warmup=args.warmup,
                            samples=args.samples)
    out['reported_sigma_fit'] = ec.fit(d, sigma=None, chains=args.chains,
                                       warmup=args.warmup, samples=args.samples)
    out['reported_sigma_fit'].pop('sigma_by_group', None)
    dump('ipd_ecological', out, args.quiet)
    return out


def cmd_patients(args):
    out = pl.summary()
    out['verification_table'] = pl.verify().to_dict('records')
    out['chan_outcome_contrast'] = pl.chan_outcome_contrast()
    out['ogura_deaths'] = pl.ogura_deaths()
    dump('ipd_patient_level', out, args.quiet)
    return out


def cmd_figures(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': ['Liberation Sans', 'Arial', 'DejaVu Sans'],
                         'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False})
    p = os.path.join(OUT, 'ipd_ecological.json')
    out = json.load(open(p)) if os.path.exists(p) else cmd_ecological(args)
    sw = out['sweep']
    x = np.array([s['sigma_mmHg'] for s in sw])
    y = np.array([s['mu_beta'] for s in sw])
    lo = np.array([s['mu_beta_ci'][0] for s in sw])
    hi = np.array([s['mu_beta_ci'][1] for s in sw])
    fig, ax = plt.subplots(figsize=(3.6, 2.7))
    ax.axhline(0, color='0.7', lw=0.8)
    ax.fill_between(x, lo, hi, color='#1f5fbf', alpha=0.15)
    ax.plot(x, y, 'o-', color='#1f5fbf', lw=1.4, ms=3.5)
    ax.axhline(y[0], color='#c62828', lw=0.9, ls=':')
    for sd in out['reported_within_group_sd_mmHg']:
        ax.axvline(sd, color='0.6', lw=0.6, ls='--')
    ax.set_xlabel('assumed within-group SD of the gradient [mmHg]')
    ax.set_ylabel(r'slope for a patient, $\mu_\beta$')
    ax.set_title('Ecological bias is small and runs one way', fontsize=8)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, 'ipd_ecological.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(OUT, 'ipd_ecological.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved ipd_ecological')

    # the published patient-level table, one point per recipient
    d = pl.load()
    paired = d.dropna(subset=['gradient_before', 'gradient_after'])
    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    for _, r in paired.iterrows():
        bad = bool(r.sfss or r.graft_loss_or_death)
        ax.plot([0, 1], [r.gradient_before, r.gradient_after],
                '-o', ms=3, lw=1.6 if bad else 0.9,
                color='#c62828' if bad else '0.55', zorder=3 if bad else 2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['shunt clamped', 'shunt open'])
    ax.set_ylabel('portocaval gradient [mmHg]')
    ax.set_xlim(-0.25, 1.35)
    ax.axhline(10, color='#1f5fbf', lw=0.8, ls=':')
    ax.text(1.05, 10.4, r'$z_P=2$', color='#1f5fbf', fontsize=6)
    ax.text(1.05, 15, 'SFSS,\ngraft loss', color='#c62828', fontsize=6, va='center')
    ax.set_title('Each recipient of a hemiportocaval shunt', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'ipd_botha_patients.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(OUT, 'ipd_botha_patients.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('saved ipd_botha_patients')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('what', nargs='*', default=['all'])
    ap.add_argument('--chains', type=int, default=6)
    ap.add_argument('--warmup', type=int, default=8000)
    ap.add_argument('--samples', type=int, default=25000)
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args(argv)
    what = ['patients', 'ecological', 'figures'] if 'all' in a.what else a.what
    jobs = dict(patients=cmd_patients, ecological=cmd_ecological, figures=cmd_figures)
    for w in what:
        if w not in jobs:
            raise SystemExit('unknown target %r' % w)
        print('\n=== %s ===' % w)
        jobs[w](a)


if __name__ == '__main__':
    main()
