"""Individual patient data recovered from a published Kaplan-Meier curve.

A survival curve is not a summary. A Kaplan-Meier step function, together with the numbers
at risk printed under it, is close to a complete record: every step is an event, the width
between steps and the drop in the risk table fix how many patients were censored in
between, and the total event count closes the system. Guyot et al. (BMC Med Res Methodol
2012;12:9) showed that the event and censoring times can be solved back out of those three
ingredients, with an error that is small compared with the digitisation itself.

That is a different thing from inventing patients. Nothing here creates information: the
algorithm inverts a function that the authors of the original paper published. It is the
accepted route from "the series reports 1-year graft survival of 86%" to a survival data
set that can be pooled, stratified or regressed on a covariate.

What this module needs, per curve, is what a reader can measure off the figure:

    points   (t, S) read off the step function, in order of increasing t
    risk     (t, n) the numbers at risk printed under the axis
    events   the total number of events, when the paper states it

and what it returns is one row per patient, (time, status), plus a check that the
reconstructed curve reproduces the digitised one.

None of the LiPNet series has been digitised yet. The reconstruction is validated here
against simulated cohorts whose truth is known (tests/test_ipd.py), so that when the curves
are read off the published figures the method is already checked rather than trusted.
"""
import numpy as np


def kaplan_meier(time, status):
    """KM estimate from individual data. Returns (times, survival, n_at_risk)."""
    time = np.asarray(time, dtype=float)
    status = np.asarray(status, dtype=int)
    order = np.argsort(time, kind='mergesort')
    time, status = time[order], status[order]
    ts = np.unique(time[status == 1])
    n = len(time)
    s, out_s, out_n = 1.0, [], []
    for t in ts:
        at_risk = int(np.sum(time >= t))
        d = int(np.sum((time == t) & (status == 1)))
        s *= (1 - d / at_risk) if at_risk else 1.0
        out_s.append(s)
        out_n.append(at_risk)
    return ts, np.array(out_s), np.array(out_n)


def reconstruct(points, risk, total_events=None, tol=1e-6, max_iter=200):
    """Guyot reconstruction.

    points : sequence of (t, S) read off the published curve, t increasing, S decreasing
    risk   : sequence of (t, n) from the numbers-at-risk row, t increasing
    total_events : the paper's reported number of events, if given

    Returns dict with time, status (one entry per patient), and the fit of the
    reconstructed curve to the digitised points.
    """
    pts = np.asarray(sorted(points, key=lambda p: p[0]), dtype=float)
    t_d, s_d = pts[:, 0], np.minimum.accumulate(pts[:, 1])
    rk = np.asarray(sorted(risk, key=lambda p: p[0]), dtype=float)
    t_r, n_r = rk[:, 0], rk[:, 1]
    N = len(t_d)
    M = len(t_r)

    # which digitised points fall in which risk interval
    lower, upper = [], []
    for i in range(M):
        hi = t_r[i + 1] if i + 1 < M else np.inf
        idx = np.flatnonzero((t_d >= t_r[i]) & (t_d < hi))
        lower.append(idx[0] if len(idx) else -1)
        upper.append(idx[-1] if len(idx) else -2)

    d = np.zeros(N)
    cen = np.zeros(N)
    n_hat = np.zeros(N + 1)
    n_hat[0] = n_r[0]
    s_last = 1.0
    k0_done = 0

    for i in range(M):
        lo, hi = lower[i], upper[i]
        if lo < 0 or hi < lo:
            continue
        n_start = n_hat[lo] if lo == k0_done else n_hat[lo]
        n_end_target = n_r[i + 1] if i + 1 < M else None
        n_in_interval = hi - lo + 1

        def run(cen_total):
            """Walk the interval with a given total censored, return (d, cen, n_end)."""
            dd = np.zeros(n_in_interval)
            cc = np.zeros(n_in_interval)
            # censoring spread evenly over the interval
            if cen_total > 0:
                base = int(cen_total) // n_in_interval
                cc[:] = base
                rem = int(round(cen_total)) - base * n_in_interval
                if rem > 0:
                    step = max(n_in_interval // max(rem, 1), 1)
                    for j in range(rem):
                        cc[min(j * step, n_in_interval - 1)] += 1
            nn = n_start
            s_prev = s_last
            for j in range(n_in_interval):
                k = lo + j
                if s_prev <= 0 or nn <= 0:
                    dd[j] = 0.0
                else:
                    dd[j] = max(round(nn * (1.0 - s_d[k] / s_prev)), 0.0)
                    if dd[j] > 0:
                        s_prev = s_prev * (1.0 - dd[j] / nn)
                nn = nn - dd[j] - cc[j]
            return dd, cc, nn, s_prev

        if n_end_target is None:
            # last interval: censoring rate carried over, or forced by the event total
            if total_events is not None:
                already = d[:lo].sum()
                target = max(total_events - already, 0)
                lo_c, hi_c = 0.0, max(n_start, 1.0)
                for _ in range(max_iter):
                    mid = 0.5 * (lo_c + hi_c)
                    dd, cc, nn, sp = run(mid)
                    if dd.sum() > target:
                        lo_c = mid
                    else:
                        hi_c = mid
                    if hi_c - lo_c < 1e-3:
                        break
                cen_total = 0.5 * (lo_c + hi_c)
            else:
                cen_total = max(n_start - 0.0, 0.0)
                dd, cc, nn, sp = run(0.0)
                cen_total = max(nn, 0.0)
                dd, cc, nn, sp = run(0.0)
                cc[-1] += max(nn, 0.0)
                d[lo:hi + 1] = dd
                cen[lo:hi + 1] = cc
                s_last = sp
                n_hat[lo + 1:hi + 2] = n_start - np.cumsum(dd + cc)
                break
            dd, cc, nn, sp = run(cen_total)
            cc[-1] += max(nn, 0.0)
        else:
            # choose the total censored so the interval ends on the printed number at risk
            lo_c, hi_c = 0.0, max(n_start - n_end_target, 0.0) + 1.0
            for _ in range(max_iter):
                mid = 0.5 * (lo_c + hi_c)
                dd, cc, nn, sp = run(mid)
                if nn > n_end_target:
                    lo_c = mid
                else:
                    hi_c = mid
                if hi_c - lo_c < 1e-3:
                    break
            dd, cc, nn, sp = run(0.5 * (lo_c + hi_c))
            cc[-1] += max(nn - n_end_target, 0.0)
            nn = n_end_target

        d[lo:hi + 1] = dd
        cen[lo:hi + 1] = cc
        s_last = sp
        n_hat[lo + 1:hi + 2] = n_start - np.cumsum(dd + cc)
        n_hat[hi + 1] = nn
        k0_done = hi + 1

    # expand to one row per patient
    time, status = [], []
    for k in range(N):
        for _ in range(int(round(d[k]))):
            time.append(t_d[k])
            status.append(1)
        for _ in range(int(round(cen[k]))):
            time.append(t_d[k])
            status.append(0)
    time = np.array(time, dtype=float)
    status = np.array(status, dtype=int)

    fit = _fit_quality(time, status, t_d, s_d)
    return dict(time=time, status=status,
                n_patients=int(len(time)), n_events=int(status.sum()),
                n_risk_reported=float(n_r[0]),
                events_reported=None if total_events is None else int(total_events),
                max_abs_error=fit['max_abs_error'], rmse=fit['rmse'])


def _fit_quality(time, status, t_d, s_d):
    if len(time) == 0:
        return dict(max_abs_error=float('nan'), rmse=float('nan'))
    ts, ss, _ = kaplan_meier(time, status)
    rec = np.array([ss[ts <= t][-1] if np.any(ts <= t) else 1.0 for t in t_d])
    err = np.abs(rec - s_d)
    return dict(max_abs_error=float(err.max()), rmse=float(np.sqrt((err ** 2).mean())))


def digitise(time, status, n_points=40, t_max=None):
    """Turn a known cohort into the (t, S) points a reader would measure off its figure.
    Used to test the reconstruction against a cohort whose truth is known."""
    ts, ss, _ = kaplan_meier(time, status)
    if len(ts) == 0:
        return []
    t_max = t_max or float(np.max(time))
    grid = np.linspace(float(ts[0]), t_max, n_points)
    out = []
    for t in grid:
        s = ss[ts <= t][-1] if np.any(ts <= t) else 1.0
        out.append((float(t), float(s)))
    return out


def risk_table(time, n_times=5, t_max=None):
    """The numbers-at-risk row a paper prints under the curve."""
    time = np.asarray(time, dtype=float)
    t_max = t_max or float(np.max(time))
    grid = np.linspace(0.0, t_max, n_times)
    return [(float(t), float(np.sum(time >= t))) for t in grid]


# ------------------------------------------------------------------ reading a digitised curve
def load_curve_csv(path):
    """Read one digitised curve in the format of data/curves/TEMPLATE_study_group.csv.

    Metadata lines begin with '#' and carry 'key value'; 'risk' holds the numbers-at-risk
    row as 'time:n, time:n, ...'. Everything after the header line 'time,survival' is the
    step function as it was measured off the figure.
    """
    meta, pts = {}, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                parts = line.lstrip('#').strip().split(None, 1)
                if len(parts) == 2:
                    meta[parts[0]] = parts[1].strip()
                continue
            if line.lower().startswith('time'):
                continue
            a, b = line.split(',')[:2]
            pts.append((float(a), float(b)))
    risk = []
    for item in meta.get('risk', '').split(','):
        if ':' in item:
            t, n = item.split(':')
            risk.append((float(t), float(n)))
    meta['points'] = pts
    meta['risk'] = risk
    return meta


def reconstruct_csv(path):
    """Reconstruct the patients behind one digitised curve file."""
    m = load_curve_csv(path)
    if not m['risk']:
        raise ValueError('%s has no numbers-at-risk row; the cohort size cannot be fixed' % path)
    tot = int(m['total_events']) if 'total_events' in m else None
    r = reconstruct(m['points'], m['risk'], total_events=tot)
    r['study'] = m.get('study')
    r['group'] = m.get('group')
    r['gradient_mmHg'] = float(m['gradient_mmHg']) if 'gradient_mmHg' in m else None
    r['time_unit'] = m.get('time_unit')
    return r
