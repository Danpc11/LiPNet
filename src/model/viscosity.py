"""Blood viscosity: the resistance of the same network perfused with the blood of the recipient.

The loads are normalised to a healthy liver perfused with normal blood (90 mL/min/100 g, 5 mmHg). The resistance of a
network is proportional to the viscosity of the blood in each of its vessels, and a recipient is anaemic and diluted
when portal flow and pressure are measured. This module gives the factor by which the resistance of the network changes
when the haematocrit goes from the reference to that of the recipient.

law     Pries, Secomb, Gessner et al. (Circ Res 1994;75:904): apparent viscosity of blood in vivo relative to plasma, as
        a function of vessel diameter and discharge haematocrit. Applied to each generation of the calibrated human
        network (`human.liver`) at its own diameter and weighted by the share of the pressure drop of that generation.
        Sinusoids take the law at their diameter (fem.params). No parameter is fitted here.
check   Kamikado, Shibamoto, Hongo and Koyama (Am J Physiol Heart 2004;286:H121): isolated rat livers perfused through
        the portal vein with blood at haematocrit 30, 20, 10 and 0% at five flows; total and segmental resistances
        digitised from their Fig. 4 (data/kamikado2004.tsv; the one set the text gives, 30% at 15 mL/min, is recovered
        within 0.002). The law was not derived from these livers.
ranges  recipient 0.25 to 0.35: in 332 liver recipients the haematocrit was 36.8% before the operation and 31.6% (SD
        4.6) after it, with red cells given below 25% during it (Choi et al., Int J Med Sci 2016;13:235). The series
        of data/series.tsv available in full text do not report it at the time of the measurement. Reference 0.40
        to 0.45, normal blood.

    python src/model/viscosity.py     -> results/model_viscosity.json
"""
import json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from model import human
from fem.params import Flow

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
HCT_REFERENCE = (0.40, 0.45)
HCT_RECIPIENT = (0.25, 0.35)
SINUSOID_UM = Flow().d_sin * 1e6
B_HUMAN = (0.6, 2 / 3, 0.75, 0.85, 0.93)        # the calibrated human livers of the rigid range (compartments.geometries)
RAT_SCALE = (10.0 / 1500.0) ** (1 / 3)          # linear scale of a 10 g liver against the 1500 g one


def relative_viscosity(D, H):
    """Apparent viscosity in vivo relative to plasma; D in micrometres, H the discharge haematocrit (Pries et al. 1994)."""
    D = np.asarray(D, float); H = np.asarray(H, float)
    e45 = 6 * np.exp(-0.085 * D) + 3.2 - 2.44 * np.exp(-0.06 * D ** 0.645)
    k = 1 / (1 + 1e-11 * D ** 12); C = (0.8 + np.exp(-0.075 * D)) * (-1 + k) + k; g = (D / (D - 1.1)) ** 2
    return (1 + (e45 - 1) * ((1 - H) ** C - 1) / ((1 - 0.45) ** C - 1) * g) * g


def _network(b, scale=1.0):
    h = human.liver(b); D = 2000.0 * np.asarray(h['r_mm']) * scale; k = int(h['g']); D[k] = SINUSOID_UM
    return D, h['dP_mmHg'] / h['dP_mmHg'].sum(), k


def factor(H, H_ref, b=0.75, scale=1.0, segments=False):
    """Resistance of the network at haematocrit H over that at H_ref. With `segments`, also the factor of the portal
    tree, the sinusoids and the hepatic venous tree."""
    D, w, k = _network(b, scale); r = relative_viscosity(D, H) / relative_viscosity(D, H_ref); total = float((w * r).sum())
    if not segments: return total
    return dict(total=total, portal=float((w[:k] * r[:k]).sum() / w[:k].sum()), sinusoidal=float(r[k]), hepatic_venous=float((w[k + 1:] * r[k + 1:]).sum() / w[k + 1:].sum()))


def factor_range(recipient=HCT_RECIPIENT, reference=HCT_REFERENCE, bs=B_HUMAN):
    """Range of the factor over the two haematocrit ranges and the calibrated networks."""
    z = [factor(H, Hr, b) for H in recipient for Hr in reference for b in bs]
    return float(min(z)), float(max(z))


def draw(rng, n, recipient=HCT_RECIPIENT, reference=HCT_REFERENCE, bs=B_HUMAN):
    """n draws of the factor: both haematocrits uniform in their ranges, one of the calibrated networks each."""
    nets = [_network(b) for b in bs]; H = rng.uniform(*recipient, n); Hr = rng.uniform(*reference, n); i = rng.integers(0, len(nets), n); out = np.empty(n)
    for j, (D, w, k) in enumerate(nets):
        m = i == j
        if m.any(): out[m] = (w[None, :] * relative_viscosity(D[None, :], H[m, None]) / relative_viscosity(D[None, :], Hr[m, None])).sum(1)
    return out


def published_check():
    """The law against the perfused rat livers of Kamikado et al. 2004: resistance at each haematocrit over that at
    30%, at each of the five flows, total and by segment; and what the law gives for a large vessel, for the human
    network and for the same network at the size of a rat liver."""
    k = pd.read_csv(os.path.join(ROOT, 'data', 'kamikado2004.tsv'), sep='\t'); out = dict(states={}, source='Kamikado et al. 2004, Fig. 4, digitised')
    seg = {'total': 'R_total', 'portal': 'R_portal', 'sinusoidal': 'R_sinusoidal', 'hepatic_venous': 'R_hepatic_venous'}; err = []
    for H in (0.20, 0.10, 0.0):
        obs = {name: [round(float(k[(k.flow_mL_min == f) & np.isclose(k.haematocrit, H)][col].iloc[0] / k[(k.flow_mL_min == f) & np.isclose(k.haematocrit, 0.30)][col].iloc[0]), 3)
                      for f in sorted(k.flow_mL_min.unique())] for name, col in seg.items()}
        pred = {'large vessel': round(float(relative_viscosity(1e4, H) / relative_viscosity(1e4, 0.30)), 3),
                'human network': {n: round(v, 3) for n, v in factor(H, 0.30, segments=True).items()},
                'network at rat size': {n: round(v, 3) for n, v in factor(H, 0.30, scale=RAT_SCALE, segments=True).items()}}
        out['states'][f'haematocrit {H:g} over 0.30'] = dict(observed_by_flow=obs, observed_mean={n: round(float(np.mean(v)), 3) for n, v in obs.items()}, predicted=pred)
        err.append(np.mean(obs['total']) / pred['network at rat size']['total'])
    out['flows_mL_min'] = [int(f) for f in sorted(k.flow_mL_min.unique())]
    hs = np.round(np.linspace(0.0, 0.30, 31), 3); out['curve'] = dict(haematocrit=[float(h) for h in hs], total_at_rat_size=[round(factor(float(h), 0.30, scale=RAT_SCALE), 4) for h in hs])
    out['observed_over_predicted_total_at_rat_size'] = [round(float(e), 3) for e in err]
    out['text_check'] = dict(given_in_the_text=dict(R_portal=0.251, R_sinusoidal=0.212, R_hepatic_venous=0.105),
                             digitised={c: float(k[(k.flow_mL_min == 15) & np.isclose(k.haematocrit, 0.30)][c].iloc[0]) for c in ('R_portal', 'R_sinusoidal', 'R_hepatic_venous')})
    out['note'] = 'the highest haematocrit perfused was 30%: between the recipient and the reference (up to 45%) the law is applied beyond these livers'
    return out


def clinical(d, draws=20000):
    """The six groups that report both loads against the model perfused with the blood of the recipient: the regions
    of `distension.bracket` with the resistance of the rigid network multiplied by the factor, and the joint draws of
    `distension.joint` with both haematocrits drawn in their ranges."""
    from model import distension as ds
    v = factor_range(); plain = ds.bracket(d); visc = ds.bracket(d, viscosity=v); j0 = ds.joint(d, draws=draws); j1 = ds.joint(d, draws=draws, haematocrit=(HCT_RECIPIENT, HCT_REFERENCE))
    rows = []
    for a, b, x, y in zip(plain['groups'], visc['groups'], j0['groups'], j1['groups']):
        rows.append(dict(study=a['study'], group=a['group'], zF=a['zF'], zR_observed=a['zR_observed'], zR_observed_ci=a['zR_observed_ci'],
                         normal_blood=dict(rigid=a['zR_rigid'], keeps_distending=a['zR_keeps_distending'], stops_distending=a['zR_stops_distending'], placement=a['placement']),
                         blood_of_the_recipient=dict(rigid=b['zR_rigid'], keeps_distending=b['zR_keeps_distending'], stops_distending=b['zR_stops_distending'], placement=b['placement']),
                         joint=dict(normal_blood={m: z['observed_over_model'] for m, z in x['models'].items()}, blood_of_the_recipient={m: z['observed_over_model'] for m, z in y['models'].items()},
                                    prob_below_rigid=dict(normal_blood=x['models']['rigid network']['prob_observed_below'], blood_of_the_recipient=y['models']['rigid network']['prob_observed_below']),
                                    prob_below_every_relation=dict(normal_blood=x['prob_observed_below_every_relation'], blood_of_the_recipient=y['prob_observed_below_every_relation']))))
    obs = np.array([r['zR_observed'] for r in rows]); rig = np.array([np.mean(r['blood_of_the_recipient']['rigid']) for r in rows]); rig0 = np.array([np.mean(r['normal_blood']['rigid']) for r in rows])
    return dict(factor=[round(v[0], 3), round(v[1], 3)], haematocrit_recipient=list(HCT_RECIPIENT), haematocrit_reference=list(HCT_REFERENCE), groups=rows,
                counts=dict(normal_blood=plain['counts'], blood_of_the_recipient=visc['counts']),
                rmse_rigid_midrange=dict(normal_blood=round(float(np.sqrt(np.mean((rig0 - obs) ** 2))), 3), blood_of_the_recipient=round(float(np.sqrt(np.mean((rig - obs) ** 2))), 3)),
                note='the haematocrit of the recipients is a range taken from other recipients, not a measurement of these groups')


def summary():
    import indices
    d = indices.compute(pd.read_csv(os.path.join(ROOT, 'data', 'series.tsv'), sep='\t'))
    table = {f'{Hr:g}': {f'{H:g}': [round(min(factor(H, Hr, b) for b in B_HUMAN), 3), round(max(factor(H, Hr, b) for b in B_HUMAN), 3)] for H in (0.20, 0.25, 0.30, 0.35)} for Hr in (0.40, 0.45)}
    return dict(law='Pries et al. 1994, in vivo', sinusoid_um=SINUSOID_UM, large_vessel_relative_viscosity={f'{H:g}': round(float(relative_viscosity(1e4, H)), 3) for H in (0.0, 0.25, 0.30, 0.35, 0.40, 0.45)},
                factor_by_haematocrit=table, segments_at_0p30_over_0p42={n: round(x, 3) for n, x in factor(0.30, 0.42, segments=True).items()}, published_check=published_check(), clinical=clinical(d))


if __name__ == '__main__':
    s = summary(); os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    json.dump(s, open(os.path.join(ROOT, 'results', 'model_viscosity.json'), 'w'), indent=1)
    print('factor on the resistance of the network, recipient over reference:', s['clinical']['factor'])
    for k, v in s['published_check']['states'].items(): print(f"  rat livers, {k}: observed {v['observed_mean']['total']} ({min(v['observed_by_flow']['total'])} to {max(v['observed_by_flow']['total'])}), law at rat size {v['predicted']['network at rat size']['total']}, large vessel {v['predicted']['large vessel']}")
    for r in s['clinical']['groups']: print(f"  {r['study']} / {r['group']}: observed {r['zR_observed']} {r['zR_observed_ci']}; rigid {r['normal_blood']['rigid']} -> {r['blood_of_the_recipient']['rigid']}; relations {r['blood_of_the_recipient']['keeps_distending']} and {r['blood_of_the_recipient']['stops_distending']}")
    print('  results/model_viscosity.json written')
