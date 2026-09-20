"""Scaling of the optimal canopy-to-canopy network with size N, inflow F and budget C.

Optimal tree (rank-one loading, concave cost): D* = C^{-1/alpha} S^{1+1/alpha},
    S = sum_e a_e^{1/(1+alpha)} f_e^{2alpha/(1+alpha)},  alpha = b/2,  a_e = l_e^{3b/2}.
Set-point scale: on every active edge D_e/c_e = lambda*alpha (Eq. 4 of the Letter), so
    D = lambda*alpha*C  and  tau0 = sqrt(lambda*alpha) = sqrt(D*/C).
Vessel mass: c_e = m_e^b (with a_e = l^{3b/2}, alpha=b/2: a w^alpha = (r^2 l)^b), and at the
optimum c_e ∝ s_e (the summand of S), so m_e = (C s_e/S)^{1/b}, V = sum_e m_e.

Symmetric space-filling canopy in d dimensions: n=3 branching, g generations, lobule
(terminal) length L_h fixed, so l_i = L_h * 3^{(g-1-i)/d}, N = 3^{g-1} lobules; the lobule
layer is N channels of length L_h and flow F/N; the outlet tree mirrors the inlet tree.
"""
import numpy as np, pandas as pd

def canopy_S_terms(g, b, d, F=1.0, Lh=1.0, n=3.0, Lsin=None):
    alpha = b / 2; q = 1 / (1 + alpha); i = np.arange(g)
    cnt = n ** i; f = F / cnt; l = Lh * n ** ((g - 1 - i) / d)
    N = n ** (g - 1)
    Lsin = Lh if Lsin is None else Lsin
    cnt = np.concatenate([cnt, [N], cnt[::-1]]); f = np.concatenate([f, [F / N], f[::-1]])
    l = np.concatenate([l, [Lsin], l[::-1]])
    a = l ** (1.5 * b)
    s = cnt * a ** q * f ** (2 * alpha * q)          # per-generation contribution to S
    return dict(N=N, cnt=cnt, f=f, l=l, s=s, S=s.sum())

def canopy_optimum(g, b, d, F=1.0, C=1.0, **kw):
    alpha = b / 2; t = canopy_S_terms(g, b, d, F, **kw)
    D = C ** (-1 / alpha) * t['S'] ** (1 + 1 / alpha)
    c = C * t['s'] / t['S']                           # cost per generation (all edges)
    m_e = (c / t['cnt']) ** (1 / b); V = np.sum(t['cnt'] * m_e)
    # radii from m_e = r^2 l ; shear tau = f / r^3
    r = np.sqrt(m_e / t['l']); tau = t['f'] / r ** 3
    k = len(r) // 2; tree = np.r_[tau[:k], tau[k + 1:]]         # tree generations only (lobule layer excluded)
    return dict(N=t['N'], D=D, tau0=np.sqrt(D / C), V=V, r_root=r[0], r_sin=r[k],
                tau_root=tau[0], tau_sin=tau[k], tau_ratio=tree.max() / tree.min(), tau_sin_over_root=tau[k] / tau[0])

def fit_exp(x, y):
    return np.polyfit(np.log(x), np.log(y), 1)[0]

if __name__ == '__main__':
    import sys, os
    OUT = sys.argv[1] if len(sys.argv) > 1 else 'results'; os.makedirs(OUT, exist_ok=True)
    rows = []
    gs = np.arange(6, 15)                             # N = 3^5 ... 3^13 ≈ 1.6e6 lobules
    for d in (2, 3):
        for b in (1.0, 0.75, 2 / 3, 0.5):
            R = [canopy_optimum(g, b, d) for g in gs]
            N = np.array([r['N'] for r in R])
            D = np.array([r['D'] for r in R]); V = np.array([r['V'] for r in R]); tr = np.array([r['tau_ratio'] for r in R])
            hi = slice(4, None)                        # large-N slope
            # exponents at fixed F, C:
            eD = fit_exp(N[hi], D[hi]); eV = fit_exp(N[hi], V[hi]); eTr = fit_exp(N[hi], tr[hi])
            # analytic: S ∝ N^{(2-b)/(2+b)} when terminal generation dominates
            alpha = b / 2
            eS_pred = (2 - b) / (2 + b) if (2 - b - 3 * b / d) > 0 else 3 * b / (d * (2 + b))
            eD_pred = eS_pred * (1 + 1 / alpha)
            rows.append(dict(d=d, b=b, e_D_N=eD, e_D_pred=eD_pred, e_V_N=eV, e_tauratio_N=eTr,
                             tau_ratio_1e6=canopy_optimum(14, b, d)['tau_ratio']))
    df = pd.DataFrame(rows); df.to_csv(f'{OUT}/scaling_exponents.tsv', sep='\t', index=False, float_format='%.4f')
    print(df.round(3).to_string(index=False))

    # ---------- allometric closure ----------
    # D* ∝ F^2 C^{-2/b} N^{eD};  tau0^2 = D*/C ∝ F^2 C^{-(2+b)/b} N^{eD};  V ∝ C^{1/b} N^{eV}
    # Closure 1 (blood volume ∝ liver volume): V ∝ N  =>  C ∝ N^{b(1-eV)}
    # Closure 2 (shear set-point invariant): tau0 = const  =>  F ∝ C^{(2+b)/(2b)} N^{-eD/2}
    print('\nAllometric closure (V ∝ N, tau0 invariant):  F ∝ N^gamma ;  empirical Kleiber+liver: F ∝ M^0.75, N ∝ M^0.87 => gamma_emp = 0.86')
    out = []
    for _, r in df.iterrows():
        b = r.b; eC = b * (1 - r.e_V_N)
        gamma = eC * (2 + b) / (2 * b) - r.e_D_N / 2
        # alternative: impose F ∝ M^0.75, N ∝ M^0.87, V ∝ N; predict tau0 ∝ M^x
        x = 0.75 - (2 + b) / (2 * b) * eC * 0.87 + r.e_D_N / 2 * 0.87
        # and the liver-mass exponent that WOULD make tau0 invariant given F ∝ M^0.75:
        # 0 = 0.75 - [(2+b)/(2b) eC - eD/2] * y  =>  y = 0.75 / gamma
        out.append(dict(d=r.d, b=b, C_exp_N=eC, gamma_F_vs_N=gamma, tau0_exp_M_given_Kleiber_and_0p87=x,
                        liver_exp_M_for_invariant_tau0=0.75 / gamma if gamma != 0 else np.nan))
    do = pd.DataFrame(out); do.to_csv(f'{OUT}/scaling_allometry.tsv', sep='\t', index=False, float_format='%.4f')
    print(do.round(3).to_string(index=False))
