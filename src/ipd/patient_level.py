"""Patient-level data that was published and never used as such.

Some series print one row per recipient. Botha et al. (Liver Transpl 2010;16:649-657)
Table 5 lists all sixteen recipients of a left lobe graft with a hemiportocaval shunt, with
the portocaval gradient before and after the shunt was unclamped, the graft-to-recipient
weight ratio, encephalopathy, small-for-size syndrome, follow-up and graft loss. Nothing
has to be reconstructed, digitised or imputed: the authors already did the work, and the
meta-analysis then collapsed it into one row with a median.

What that row throws away is the whole point of the model. The shunt is a deliberate,
within-patient manipulation of the resistance load: the flow through the graft is barely
touched, the gradient falls from a median of 18 to 5 mmHg, and the one recipient whose
gradient stayed high is the one who developed small-for-size syndrome and lost the graft.

Two things this module does with it:

transcription control
    Every summary statistic the paper reports is recomputed from the transcribed table.
    Agreement is the evidence that the transcription is right; the one disagreement is
    reported rather than quietly corrected.

measured dispersion
    The within-group standard deviation of the gradient, which the ecological model in
    ecological.py otherwise has to assume for this group, is computed directly.
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
BOTHA = os.path.join(ROOT, 'data', 'patient_level', 'botha2010.tsv')
MMHG_PER_ZP = 5.0

# what the paper states in its own tables and text, for the transcription control
BOTHA_REPORTED = dict(n=16, before_n=11, before_median=18.0, before_range=(12.0, 25.0),
                      after_median=5.0, after_range=(1.0, 15.0),
                      gwrw_median=0.67, gwrw_range=(0.51, 1.06), pct_slv_median=28.5,
                      encephalopathy=10, sfss=1, graft_loss_or_death=3)


def load(path=BOTHA):
    d = pd.read_csv(path, sep='\t', comment='#')
    d['zP_before'] = d.gradient_before / MMHG_PER_ZP
    d['zP_after'] = d.gradient_after / MMHG_PER_ZP
    return d


def verify(d=None, reported=None):
    """Recompute every statistic the paper reports. Returns one row per check."""
    d = load() if d is None else d
    r = reported or BOTHA_REPORTED
    b, a = d.gradient_before.dropna(), d.gradient_after.dropna()
    checks = [
        ('patients', len(d), r['n']),
        ('gradients before the shunt', len(b), r['before_n']),
        ('median gradient before', float(b.median()), r['before_median']),
        ('minimum gradient before', float(b.min()), r['before_range'][0]),
        ('maximum gradient before', float(b.max()), r['before_range'][1]),
        ('median gradient after', float(a.median()), r['after_median']),
        ('minimum gradient after', float(a.min()), r['after_range'][0]),
        ('maximum gradient after', float(a.max()), r['after_range'][1]),
        ('median GW/RW', round(float(d.gwrw.median()), 2), r['gwrw_median']),
        ('median %SLV', round(float(d.pct_slv.median()), 1), r['pct_slv_median']),
        ('encephalopathy', int(d.encephalopathy.sum()), r['encephalopathy']),
        ('small-for-size syndrome', int(d.sfss.sum()), r['sfss']),
        ('graft loss or death', int(d.graft_loss_or_death.sum()), r['graft_loss_or_death']),
    ]
    out = pd.DataFrame(checks, columns=['quantity', 'transcribed', 'reported'])
    out['agrees'] = np.isclose(out.transcribed.astype(float),
                               out.reported.astype(float), atol=0.011)
    return out


def shunt_effect(d=None):
    """The within-patient manipulation: what the shunt does to each load."""
    d = load() if d is None else d
    paired = d.dropna(subset=['gradient_before', 'gradient_after'])
    before, after = paired.gradient_before.values, paired.gradient_after.values
    return dict(n_paired=int(len(paired)),
                before_median=float(np.median(before)),
                after_median=float(np.median(after)),
                zP_before_median=float(np.median(before)) / MMHG_PER_ZP,
                zP_after_median=float(np.median(after)) / MMHG_PER_ZP,
                drop_median=float(np.median(before - after)),
                all_fell=bool(np.all(after < before)),
                zP_after_min=float(after.min()) / MMHG_PER_ZP,
                zP_after_max=float(after.max()) / MMHG_PER_ZP)


def event_rank_test(d=None, outcome='sfss'):
    """Where the adverse outcomes sit in the ranking of the achieved pressure load.

    With one event among ten measured recipients no model can be fitted, so the question
    is asked exactly: under the null that the outcome is independent of the load, the
    probability that the events occupy the highest ranks is a hypergeometric tail. This is
    a small-sample statement and is reported as one.
    """
    d = load() if d is None else d
    sub = d.dropna(subset=['gradient_after'])
    g = sub.gradient_after.values
    ev = sub[outcome].values.astype(int)
    n, k = len(g), int(ev.sum())
    if k == 0:
        return dict(n=n, events=0, note='no events among the recipients with a gradient')
    order = np.argsort(-g)                      # highest gradient first
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    event_ranks = sorted(ranks[ev == 1].tolist())
    from math import comb
    # P(all k events fall in the top k of n) under random assignment
    p_top = 1.0 / comb(n, k)
    return dict(n=n, events=k,
                gradients=[float(x) for x in g],
                event_gradients=[float(x) for x in g[ev == 1]],
                event_ranks=event_ranks,
                highest_gradient=float(g.max()),
                events_are_the_top_k=bool(event_ranks == list(range(1, k + 1))),
                p_events_in_top_k=float(p_top))


def group_row(d=None):
    """The group-level row this series should contribute, computed rather than read off a
    median: the mean achieved gradient and its measured within-group dispersion."""
    d = load() if d is None else d
    a = d.gradient_after.dropna()
    return dict(study='Botha 2010', n=int(len(d)),
                gradient_mmHg=float(a.mean()), gradient_sd=float(a.std(ddof=1)),
                gradient_n_measured=int(len(a)),
                gradient_basis='computed from the published patient-level table',
                zP=float(a.mean()) / MMHG_PER_ZP,
                sigma_zP=float(a.std(ddof=1)) / MMHG_PER_ZP,
                events=int(d.sfss.sum()), outcome='SFSS')


def summary(d=None):
    d = load() if d is None else d
    v = verify(d)
    out = dict(source='Botha 2010, Liver Transpl 16:649-657, Table 5',
               checks_passed=int(v.agrees.sum()), checks_total=int(len(v)),
               disagreements=[{'quantity': r.quantity, 'transcribed': r.transcribed,
                               'reported': r.reported}
                              for r in v[~v.agrees].itertuples()])
    out.update(shunt_effect(d))
    out['sfss'] = event_rank_test(d, 'sfss')
    out['graft_loss_or_death'] = event_rank_test(d, 'graft_loss_or_death')
    out['group_row'] = group_row(d)
    return out


# ------------------------------------------------------------------ Chan 2011
CHAN = os.path.join(ROOT, 'data', 'patient_level', 'chan2011_groups.tsv')
CHAN_DONOR_PVF = 81.0          # mL/min/100 g, right liver of the donor, same paper
CHAN_GRADIENT_MEDIAN = 6.0     # mmHg after graft implantation, whole cohort
CHAN_GRADIENT_RANGE = (-1.0, 11.0)


def sd_from_range(lo, hi, n):
    """Standard deviation estimated from a published range (Wan et al., BMC Med Res
    Methodol 2014): SD ~ (hi - lo) / xi(n), with xi(n) = 2 Phi^-1((n - 0.375)/(n + 0.25)).
    It is an approximation that assumes approximate normality, and it is reported as
    such, but it turns a range into one more check on the within-group spread instead of
    leaving that spread unmeasured."""
    from scipy.stats import norm
    xi = 2.0 * norm.ppf((n - 0.375) / (n + 0.25))
    return float((hi - lo) / xi)


def chan_outcome_contrast(path=CHAN):
    """The loads of the complicated and uncomplicated recipients of one series.

    Chan et al. report graft portal flow and portal pressure separately for the recipients
    who had a Clavien grade 2+ complication and those who did not. That is a within-cohort
    outcome contrast with 12 events in 46 recipients, four times the two in-hospital deaths
    the main extraction scores, and it needs no digitising. Its result is a null, which is
    the reason to include it.
    """
    d = pd.read_csv(path, sep='\t', comment='#')
    d['zF'] = d.pvf_per_100g / CHAN_DONOR_PVF
    d['zP'] = CHAN_GRADIENT_MEDIAN / MMHG_PER_ZP        # reported for the whole cohort only
    d['zR'] = d.zP / d.zF
    d['pvp_sd_est'] = [sd_from_range(r.pvp_lo, r.pvp_hi, r.n) for r in d.itertuples()]
    out = dict(groups=d.to_dict('records'),
               n_total=int(d.n.sum()), events_total=int(d.events.sum()),
               zF_uncomplicated=float(d.zF.iloc[0]), zF_complicated=float(d.zF.iloc[1]),
               zP_both=float(d.zP.iloc[0]),
               zR_uncomplicated=float(d.zR.iloc[0]), zR_complicated=float(d.zR.iloc[1]),
               pressure_separation_mmHg=float(d.pvp_after.iloc[1] - d.pvp_after.iloc[0]),
               usable_in_the_slope_model=bool(d.zP.nunique() > 1),
               cohort_gradient_sd_est=sd_from_range(CHAN_GRADIENT_RANGE[0],
                                                    CHAN_GRADIENT_RANGE[1], 46),
               cohort_pressure_sd_est=sd_from_range(10.0, 26.0, 46),
               note=('the two outcome groups have the same median portal pressure, so this '
                     'series contributes no variation in the predictor and cannot enter a '
                     'slope model; it stands as an outcome contrast with no pressure '
                     'separation'))
    # which group held the extremes of pressure
    out['pressure_range_uncomplicated'] = [float(d.pvp_lo.iloc[0]), float(d.pvp_hi.iloc[0])]
    out['pressure_range_complicated'] = [float(d.pvp_lo.iloc[1]), float(d.pvp_hi.iloc[1])]
    out['extremes_were_uncomplicated'] = bool(d.pvp_hi.iloc[0] > d.pvp_hi.iloc[1]
                                              and d.pvp_lo.iloc[0] < d.pvp_lo.iloc[1])
    return out


# ------------------------------------------------------------------ Ogura 2010
OGURA = os.path.join(ROOT, 'data', 'patient_level', 'ogura2010_deaths.tsv')
OGURA_N = {'below_15': 86, 'at_or_above_15': 43}
OGURA_DEATHS_REPORTED = {'below_15': 4, 'at_or_above_15': 12}


def ogura_deaths(path=OGURA):
    """Every recipient death of one series, with its final portal pressure.

    Ogura et al. Table 3 lists all sixteen deaths individually. Combined with the two
    denominators the paper reports, that is a patient-level outcome data set for the
    events and a group-level one for the survivors: enough to rebuild the full
    pressure-by-outcome table without a single imputed event.

    It is also the exact check that licenses digitising that paper's Fig. 4, where living
    recipients are open circles and dead recipients are crosses: any extraction of the
    crosses must reproduce these sixteen pressures.
    """
    d = pd.read_csv(path, sep='\t', comment='#')
    below = int((d.final_pvp < 15).sum())
    above = int((d.final_pvp >= 15).sum())
    checks = dict(deaths_below_15=(below, OGURA_DEATHS_REPORTED['below_15']),
                  deaths_at_or_above_15=(above, OGURA_DEATHS_REPORTED['at_or_above_15']),
                  n_deaths=(int(len(d)), 16))
    return dict(pressures=[int(x) for x in sorted(d.final_pvp)],
                n_deaths=int(len(d)),
                deaths_below_15=below, deaths_at_or_above_15=above,
                denominator_below_15=OGURA_N['below_15'],
                denominator_at_or_above_15=OGURA_N['at_or_above_15'],
                rate_below_15=below / OGURA_N['below_15'],
                rate_at_or_above_15=above / OGURA_N['at_or_above_15'],
                checks={k: {'transcribed': v[0], 'reported': v[1], 'agrees': v[0] == v[1]}
                        for k, v in checks.items()},
                note=('the crosses in that paper Fig. 4 are these sixteen recipients; any '
                      'digitisation of that figure must reproduce these pressures'))
