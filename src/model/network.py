"""Liver as a source -> extensive network -> single outlet transport system.

Formulation (same conventions as vascular-shear-setpoint/src/flows.py):
    conductance  w_e = r_e^4 / l_e          (Poiseuille, viscosity absorbed)
    shear        tau_e = |f_e| / r_e^3
    cost         C_b = sum_e a_e w_e^alpha,  alpha = b/2,  a_e = l_e^{3b/2}
                 (proportional to sum_e (r_e^2 l_e)^b: vessel mass to the power b)
    loading      one source (hilum: HA+PV inlet), one sink (hepatic vein outlet):
                 I = e_s - e_t, D = I^T L^+ I = R_st  (rank-one loading)

Geometry from Lorente, Hautefeuille & Sanchez-Cedillo, Sci. Rep. 10, 16194 (2020):
    splitting number n = 3, d_{i+1}/d_i = L_{i+1}/L_i = 3^{-1/3} (theory);
    measured (Debbaut et al.) HA 0.74/0.66, PV 0.70/0.72, HV 0.59/0.64; Ma et al. d ratio 0.79;
    lobules: hexagonal, 6 portal triads at the corners feeding one central vein,
    sinusoid path length L_h (corner to centre), about 20 generations.

Part A: analytic series-parallel model by generations (closed-form optimum on a
        fixed tree, Appendix A of the metabolic-cost paper).
Part B: explicit planar canopy-to-canopy graph with a hexagonal lobule mesh
        (loops), Laplacian dissipation, numerical optimum on the fixed support
        (two-stage softmax / L-BFGS-B + prune + refine, as in the cost-convexity paper).
"""
import numpy as np, networkx as nx
from scipy.linalg import solve
from scipy.optimize import minimize

LN3 = np.log(3.0)

# --------------------------------------------------------------------------
# Part A: symmetric generation model
# --------------------------------------------------------------------------
def gen_tree(g, rho_d, rho_l, n=3.0, d1=1.0, L1=1.0):
    """generations 1..g: (count, flow fraction, diameter, length). Flow fraction
    per branch = n^{-(i-1)}, branches = n^{i-1}."""
    i = np.arange(g)
    cnt = n ** i
    return dict(cnt=cnt, f=1.0 / cnt, d=d1 * rho_d ** i, l=L1 * rho_l ** i)

def canopy_series_parallel(inlet, outlet, lob):
    """Edge list (as arrays) for inlet tree -> N_lob parallel lobule paths -> outlet tree.
    lob = dict(N, d, l): N parallel identical channels of diameter d and length l
    (the 6 sinusoid sectors of each lobule folded into one equivalent channel)."""
    cnt = np.concatenate([inlet['cnt'], [lob['N']], outlet['cnt'][::-1]])
    f = np.concatenate([inlet['f'], [1.0 / lob['N']], outlet['f'][::-1]])
    d = np.concatenate([inlet['d'], [lob['d']], outlet['d'][::-1]])
    l = np.concatenate([inlet['l'], [lob['l']], outlet['l'][::-1]])
    return dict(cnt=cnt, f=f, d=d, l=l)

def sp_dissipation(net):
    """D = sum over generations of cnt * f^2 / w  (unit total flow)."""
    r = net['d'] / 2.0; w = r ** 4 / net['l']
    return float(np.sum(net['cnt'] * net['f'] ** 2 / w))

def sp_cost(net, b):
    r = net['d'] / 2.0; w = r ** 4 / net['l']; a = net['l'] ** (1.5 * b)
    return float(np.sum(net['cnt'] * a * w ** (b / 2.0)))

def sp_optimum(net, b, C0):
    """closed form on a fixed series-parallel support: w_e ∝ (f_e^2/a_e)^{1/(alpha+1)}."""
    alpha = b / 2.0; a = net['l'] ** (1.5 * b)
    u = (net['f'] ** 2 / a) ** (1.0 / (alpha + 1.0))
    scale = (C0 / np.sum(net['cnt'] * a * u ** alpha)) ** (1.0 / alpha)
    w = u * scale
    D = float(np.sum(net['cnt'] * net['f'] ** 2 / w))
    r = (w * net['l']) ** 0.25
    return dict(D=D, w=w, r=r, d=2 * r, tau=net['f'] / r ** 3)

def efficiency(net, b):
    """eta_b = D*(C_b(net)) / D(net): allocation efficiency of the given geometry
    on its own support and at its own budget (unique optimum on a fixed tree)."""
    C = sp_cost(net, b); Dstar = sp_optimum(net, b, C)['D']
    return Dstar / sp_dissipation(net)

def optimal_diameter_ratio(b, rho_l, n=3.0):
    """for a symmetric n-ary tree with length ratio rho_l, the b-optimal diameter
    ratio: r_{i+1}/r_i = n^{-1/(2+b)} rho_l^{(1-b)/(2(2+b))}. b=1 gives Hess-Murray."""
    return n ** (-1.0 / (2.0 + b)) * rho_l ** ((1.0 - b) / (2.0 * (2.0 + b)))

def infer_b(rho_d, rho_l, n=3.0):
    """invert optimal_diameter_ratio for b given measured (rho_d, rho_l, n)."""
    num = -np.log(n) + 0.5 * np.log(rho_l) - 2.0 * np.log(rho_d)
    den = np.log(rho_d) + 0.5 * np.log(rho_l)
    return num / den

def infer_b_mc(rho_d, sd_d, rho_l, sd_l, n, sd_n, N=20000, seed=0):
    """propagate the generation-to-generation scatter of Table 1 (Lorente) as
    uncertainty of the mean ratio: sd of the mean = sd / sqrt(#generations)."""
    rng = np.random.default_rng(seed)
    rd = rng.normal(rho_d, sd_d, N); rl = rng.normal(rho_l, sd_l, N); nn = rng.normal(n, sd_n, N)
    ok = (rd > 0.05) & (rl > 0.05) & (nn > 1.2)
    bs = infer_b(rd[ok], rl[ok], nn[ok])
    return np.percentile(bs, [16, 50, 84])

# --------------------------------------------------------------------------
# Part B: explicit planar canopy-to-canopy graph
# --------------------------------------------------------------------------
def hex_lobules(R=5.5, s=1.0):
    """triangular lattice of lobule centres inside a disk of radius R (spacing s);
    hexagon corners (portal triads) shared by neighbouring lobules.
    Returns centres, corners, and (corner, centre) sinusoid pairs, L_h = s/sqrt3."""
    Lh = s / np.sqrt(3.0)
    centres = []
    for j in range(-int(R / (s * 0.866)) - 2, int(R / (s * 0.866)) + 3):
        for i in range(-int(R / s) - 2, int(R / s) + 3):
            x = s * (i + 0.5 * (j % 2)); y = s * 0.866025 * j
            if x * x + y * y <= R * R: centres.append((x, y))
    centres = np.array(centres)
    corners = {}; pairs = []
    for c, (x, y) in enumerate(centres):
        for k in range(6):
            th = np.pi / 6 + k * np.pi / 3
            p = (round(x + Lh * np.cos(th), 6), round(y + Lh * np.sin(th), 6))
            if p not in corners: corners[p] = len(corners)
            pairs.append((corners[p], c))
    corner_xy = np.array(list(corners.keys()))
    return centres, corner_xy, pairs, Lh

def ternary_tree(leaf_xy, root_xy, k=3, seed=0):
    """recursive k-means clustering of the leaves into k groups; internal nodes at
    cluster centroids; root at root_xy. Returns node coordinates (root=0, leaves
    at the end in leaf order), edges (parent, child), and depth of each edge."""
    rng = np.random.default_rng(seed)
    nodes = [tuple(root_xy)]; edges = []; depth = []
    nleaf = len(leaf_xy); leaf_id = {}
    def kmeans(idx, k):
        P = leaf_xy[idx]
        if len(idx) <= k: return [[i] for i in idx]
        cen = P[rng.choice(len(P), k, replace=False)]
        for _ in range(50):
            lab = np.argmin(((P[:, None, :] - cen[None]) ** 2).sum(-1), 1)
            new = np.array([P[lab == j].mean(0) if np.any(lab == j) else cen[j] for j in range(k)])
            if np.allclose(new, cen): break
            cen = new
        groups = [[idx[t] for t in np.flatnonzero(lab == j)] for j in range(k)]
        return [g for g in groups if g]
    def split(idx, k):
        groups = kmeans(idx, k)
        if len(groups) < min(k, len(idx)):     # empty k-means cluster: split along principal axis
            P = leaf_xy[idx]; c = P - P.mean(0)
            ax = np.linalg.svd(c, full_matrices=False)[2][0]; t = c @ ax
            q = np.quantile(t, np.linspace(0, 1, k + 1)[1:-1])
            lab = np.searchsorted(q, t)
            groups = [[idx[u] for u in np.flatnonzero(lab == j)] for j in range(k)]
            groups = [g for g in groups if g]
        return groups
    def build(parent, idx, d):
        for g in split(idx, k):
            if len(g) == 1:
                nid = len(nodes); nodes.append(tuple(leaf_xy[g[0]])); leaf_id[g[0]] = nid
                edges.append((parent, nid)); depth.append(d); continue
            nid = len(nodes); nodes.append(tuple(leaf_xy[g].mean(0)))
            edges.append((parent, nid)); depth.append(d)
            build(nid, g, d + 1)
    build(0, list(range(nleaf)), 1)
    return np.array(nodes), edges, np.array(depth), leaf_id

def build_liver_graph(R=5.5, s=1.0, seed=0):
    """nodes: [source, inlet internal..., triads..., centres..., outlet internal..., sink]
    edges carry (i, j, length, kind, depth) with kind in {'in','sin','out'}."""
    centres, corners, pairs, Lh = hex_lobules(R, s)
    origin = np.zeros(2)
    nin, ein, din, lid_in = ternary_tree(corners, origin, seed=seed)
    nout, eout, dout, lid_out = ternary_tree(centres, origin, seed=seed + 1)
    xy = [nin[0]]; E = []
    off_in = 1                       # inlet internal + leaves (leaves == triads)
    for (p, c), d in zip(ein, din):
        E.append((p, c, 'in', d))
    xy = list(map(tuple, nin))
    n_in = len(nin)
    tri = {k: lid_in[k] for k in range(len(corners))}            # triad -> node
    # outlet tree: its own node set; root renamed to a new sink node; leaves are centres
    base = n_in
    cen = {}
    for k in range(len(centres)): cen[k] = base + lid_out[k]
    xy += list(map(tuple, nout))
    sink = base                       # root of outlet tree = sink node
    for (p, c), d in zip(eout, dout):
        E.append((base + p, base + c, 'out', d))
    for (q, c) in pairs:
        E.append((tri[q], cen[c], 'sin', 0))
    xy = np.array(xy)
    L = np.array([np.hypot(*(xy[i] - xy[j])) if kind != 'sin' else Lh for i, j, kind, _ in E])
    L = np.maximum(L, 1e-3)
    return dict(xy=xy, E=[(i, j) for i, j, _, _ in E], kind=np.array([k for _, _, k, _ in E]),
                depth=np.array([d for _, _, _, d in E]), L=L, src=0, sink=sink,
                n=len(xy), m=len(E), n_lob=len(centres), n_tri=len(corners), Lh=Lh)


class Loaded:
    """A subgraph with prescribed demand.  mu: node demand (sum 0), sigma: independent
    fluctuation sd on the demand nodes (balanced at the source), as in the Letter:
    Q = mu mu^T + Cov,  <D> = tr(L^+ Q),  d<D>/dw_e = -b_e^T L^+ Q L^+ b_e."""
    def __init__(self, xy, E, L, kind, depth, src, demand_nodes, sigma=0.0):
        self.xy, self.E, self.L, self.kind, self.depth = xy, list(E), np.asarray(L), np.asarray(kind), np.asarray(depth)
        self.n, self.m, self.src = len(xy), len(E), src
        self.Ei = np.array([e[0] for e in E]); self.Ej = np.array([e[1] for e in E])
        self.B = np.zeros((self.n, self.m)); self.B[self.Ei, np.arange(self.m)] = 1; self.B[self.Ej, np.arange(self.m)] = -1
        self.dem = np.array(sorted(demand_nodes)); N = len(self.dem)
        self.mu = np.zeros(self.n); self.mu[self.dem] = -1.0; self.mu[src] = N
        self.N = N; self.set_sigma(sigma)
    def set_sigma(self, sigma):
        self.sigma = sigma
        Q = np.outer(self.mu, self.mu)
        if sigma > 0:
            for c in self.dem:
                v = np.zeros(self.n); v[c] = -1; v[self.src] = 1
                Q += sigma ** 2 * np.outer(v, v)
        self.Q = Q
    def a(self, b): return self.L ** (1.5 * b)
    def cost(self, w, b, keep=None):
        keep = np.arange(self.m) if keep is None else keep
        return float(np.sum(self.a(b)[keep] * w[keep] ** (b / 2)))
    def support_ok(self, keep):
        H = nx.Graph(); H.add_edges_from(self.E[k] for k in keep)
        if self.src not in H: return False
        comp = nx.node_connected_component(H, self.src)
        return all(c in comp for c in self.dem)
    def src_component(self, keep):
        """edges of keep lying in the connected component of the source (drops dangling pieces)."""
        H = nx.Graph(); H.add_edges_from(self.E[k] for k in keep)
        comp = nx.node_connected_component(H, self.src)
        return np.array([k for k in keep if self.E[k][0] in comp])
    def solve(self, w, keep=None):
        """returns <D>, grad wrt w on keep, and mean edge flows f (for shear)."""
        keep = np.arange(self.m) if keep is None else keep
        nodes = np.unique(np.concatenate([self.Ei[keep], self.Ej[keep]]))
        mask = nodes != self.src; red = nodes[mask]
        B = self.B[np.ix_(red, keep)]
        Lr = (B * w[keep]) @ B.T
        Qr = self.Q[np.ix_(red, red)]
        Y = solve(Lr, Qr, assume_a='pos')                      # L^+ Q
        D = float(np.trace(Y))
        M = solve(Lr, Y.T, assume_a='pos')                     # L^+ Q L^+  (symmetric)
        MB = M @ B
        grad = -np.einsum('ik,ik->k', MB, B)                    # -b_e^T M b_e
        phi = solve(Lr, self.mu[red], assume_a='pos')
        f = np.zeros(self.m); f[keep] = w[keep] * (B.T @ phi)
        return D, grad, f
    def dissipation(self, w, keep=None): return self.solve(w, keep)[0]
    def cycle_rank(self, keep):
        H = nx.Graph(); H.add_edges_from(self.E[k] for k in keep)
        return H.number_of_edges() - H.number_of_nodes() + nx.number_connected_components(H)
    def d_norm(self, w, b, keep=None):
        """budget-invariant D C^{1/alpha}."""
        return self.dissipation(w, keep) * self.cost(w, b, keep) ** (2 / b)

def split_liver(G):
    """the canopy graph as two loaded subproblems sharing the budget:
    inlet = source + inlet tree + triads + sinusoids, unit sinks at the lobule centres;
    outlet = unit sources at the centres, outlet tree, single sink (as a source with N in reverse)."""
    E, kind, depth, L = G['E'], G['kind'], G['depth'], G['L']
    sel_in = np.flatnonzero(kind != 'out'); sel_out = np.flatnonzero(kind == 'out')
    centres = sorted({j for k in np.flatnonzero(kind == 'sin') for j in [E[k][1]]})
    inlet = Loaded(G['xy'], [E[k] for k in sel_in], L[sel_in], kind[sel_in], depth[sel_in], G['src'], centres)
    outlet = Loaded(G['xy'], [E[k] for k in sel_out], L[sel_out], kind[sel_out], depth[sel_out], G['sink'], centres)
    inlet.idx, outlet.idx = sel_in, sel_out
    return inlet, outlet

def murray_radii(g, w_ref=None, rho=None, sin_rule='flow'):
    """baseline radii. Tree edges: r = rho^(depth-1) (Lorente generation ratio) if rho is
    given, else Murray from the actual mean flow r ∝ f^{1/3}. Sinusoids: r ∝ f^{1/3}
    with f the mean flow (Murray uniform shear) — evaluated from a uniform-sinusoid solve."""
    w0 = np.ones(g.m) if w_ref is None else w_ref
    _, _, f = g.solve(w0)
    f = np.abs(f); fN = f / g.N
    r = fN ** (1 / 3)                                   # Murray: r^3 ∝ f, root r = 1
    if rho is not None:
        tree = g.kind != 'sin'
        r[tree] = rho ** (g.depth[tree] - 1)
    return r

def to_budget(g, r, b, C0=1.0):
    w = r ** 4 / g.L
    return w * (C0 / g.cost(w, b)) ** (2 / b)

def optimize(g, b, C0, w0_list, eps_cut=1e-5, maxiter=4000, verbose=False):
    alpha = b / 2; a = g.a(b)
    def unpack(z, idx):
        z = z - z.max(); p = np.exp(z); p /= p.sum()
        w = np.zeros(g.m); w[idx] = (C0 * p / a[idx]) ** (1 / alpha)
        return w, p
    def obj(z, idx):
        w, p = unpack(z, idx)
        D, gw, _ = g.solve(w, idx)
        gp = gw * (1 / alpha) * w[idx] / p
        return D, p * (gp - np.dot(p, gp))
    best = None
    for w0 in w0_list:
        idx = np.arange(g.m)
        p0 = a * w0 ** alpha; p0 /= p0.sum()
        res = minimize(obj, np.log(np.maximum(p0, 1e-14)), args=(idx,), jac=True, method='L-BFGS-B',
                       options=dict(maxiter=maxiter, maxfun=2 * maxiter))
        w, p = unpack(res.x, idx)
        keep = np.flatnonzero(p > eps_cut)
        if not g.support_ok(keep): keep = idx
        res2 = minimize(obj, np.log(p[keep]), args=(keep,), jac=True, method='L-BFGS-B',
                        options=dict(maxiter=maxiter))
        w2, p2 = unpack(res2.x, keep); D2 = g.dissipation(w2, keep)
        if verbose: print(f'    start D={res.fun:.6g} -> {len(keep)} edges, D={D2:.6g}, beta={g.cycle_rank(keep)}')
        if best is None or D2 < best['D']:
            best = dict(D=D2, w=w2, keep=keep, n_active=len(keep), beta=g.cycle_rank(keep),
                        n_sin=int(np.sum(g.kind[keep] == 'sin')))
    return best

def combine_budget(A_in, A_out, b, C0):
    """each subproblem at optimum obeys D_i = A_i C_i^{-1/alpha}; minimizing the sum at
    C_in + C_out = C0 gives C_in/C_out = (A_in/A_out)^{alpha/(1+alpha)}."""
    alpha = b / 2; q = alpha / (1 + alpha)
    rho = (A_in / A_out) ** q; C_in = C0 * rho / (1 + rho); C_out = C0 - C_in
    return A_in * C_in ** (-1 / alpha) + A_out * C_out ** (-1 / alpha), C_in, C_out

def shear_fit(g, w, keep, b):
    """joint fit ln tau = A + s_r ln r + s_l ln l on the active tree edges (Appendix D)."""
    _, grad, f = g.solve(w, keep)
    frms = np.zeros(g.m); frms[keep] = np.sqrt(np.maximum(w[keep] ** 2 * (-grad), 0))   # sqrt<f^2>
    kk = keep[(g.kind[keep] != 'sin') & (frms[keep] > 1e-9)]
    r = (w[kk] * g.L[kk]) ** 0.25; tau = frms[kk] / r ** 3
    X = np.c_[np.log(r), np.log(g.L[kk]), np.ones(len(kk))]
    coef = np.linalg.lstsq(X, np.log(tau), rcond=None)[0]
    return dict(s_r=coef[0], s_l=coef[1], b_r=1 + coef[0], b_l=1 + 2 * coef[1],
                tau_ratio=tau.max() / tau.min(), r=r, tau=tau, l=g.L[kk])

# --------------------------------------------------------------------------
# exact tree optimum for prescribed demand (Banavar / Appendix A of the Letter)
# --------------------------------------------------------------------------
def tree_S(g, tree_keep, b, sigma=0.0):
    """closed-form optimum on a fixed tree (edges tree_keep spanning src and demand nodes):
    with downstream demand count k_e, h_e = k_e^2 + k_e sigma^2,
    S = sum_e a_e^{1/(1+alpha)} h_e^{alpha/(1+alpha)},  D* = C0^{-1/alpha} S^{1+1/alpha}.
    Returns S, k (per edge in tree_keep order)."""
    alpha = b / 2
    T = nx.Graph(); T.add_edges_from((g.E[k][0], g.E[k][1], {'k': k}) for k in tree_keep)
    parent = dict(nx.dfs_predecessors(T, g.src)); order = list(nx.dfs_preorder_nodes(T, g.src))
    dem = set(g.dem.tolist()); sub = {v: (1.0 if v in dem else 0.0) for v in order}
    for v in reversed(order):
        if v in parent: sub[parent[v]] += sub[v]
    k = np.array([sub[g.E[e][0]] if parent.get(g.E[e][0]) == g.E[e][1] else sub[g.E[e][1]] for e in tree_keep])
    h = k ** 2 + k * sigma ** 2
    S = float(np.sum(g.a(b)[tree_keep] ** (1 / (1 + alpha)) * h ** (alpha / (1 + alpha))))
    return S, k

def tree_allocation(g, tree_keep, b, C0=1.0, sigma=0.0):
    alpha = b / 2; S, k = tree_S(g, tree_keep, b, sigma)
    h = k ** 2 + k * sigma ** 2; a = g.a(b)[tree_keep]
    u = np.where(h > 0, (h / a) ** (1 / (1 + alpha)), 0.0)
    act = u > 0
    w = np.zeros(g.m); w[tree_keep[act]] = u[act] * (C0 / np.sum(a[act] * u[act] ** alpha)) ** (1 / alpha)
    keep = tree_keep[act]
    return dict(D=C0 ** (-1 / alpha) * S ** (1 + 1 / alpha), w=w, keep=keep, n_active=len(keep), beta=0,
                n_sin=int(np.sum(g.kind[keep] == 'sin')), S=S)

def best_sinusoid_tree(g, b, sigma=0.0, sweeps=20, seed=0, init=None):
    """the inlet candidate graph is a tree plus the sinusoid mesh; a spanning tree of the
    demand is the inlet tree plus ONE sinusoid per lobule. Coordinate descent over the
    choice of feeding triad, exact allocation at each step."""
    rng = np.random.default_rng(seed)
    sin = np.flatnonzero(g.kind == 'sin'); tree_e = np.flatnonzero(g.kind != 'sin')
    by_centre = {}
    for k in sin: by_centre.setdefault(g.E[k][1], []).append(k)
    centres = sorted(by_centre)
    choice = {c: (init[c] if init is not None else rng.choice(by_centre[c])) for c in centres}
    def keep_of(ch): return np.concatenate([tree_e, np.array([ch[c] for c in centres])])
    S = tree_S(g, keep_of(choice), b, sigma)[0]
    for _ in range(sweeps):
        moved = False
        for c in rng.permutation(centres):
            cur = choice[c]
            for k in by_centre[c]:
                if k == cur: continue
                choice[c] = k; Sk = tree_S(g, keep_of(choice), b, sigma)[0]
                if Sk < S - 1e-12: S, cur, moved = Sk, k, True
                else: choice[c] = cur
        if not moved: break
    sol = tree_allocation(g, keep_of(choice), b, 1.0, sigma); sol['choice'] = choice
    return sol

# --------------------------------------------------------------------------
# local shear set-point adaptation on the fixed liver support (Eq. 6 of the Letter)
# --------------------------------------------------------------------------
def rms_shear(g, w, keep, sigma):
    """sqrt(<f_e^2>)/r_e^3 with <f_e^2> = w_e^2 b_e^T L^+ Q L^+ b_e (Q includes fluctuations)."""
    D, grad, f = g.solve(w, keep)
    f2 = w[keep] ** 2 * (-grad)
    r = (w[keep] * g.L[keep]) ** 0.25
    return np.sqrt(np.maximum(f2, 0)) / r ** 3, D, f

def adapt(g, b, r0, sigma=0.0, tau0=1.0, kappa=0.2, delta=0.02, max_steps=4000, tol=1e-5,
          prune=1e-3, verbose=False):
    """dx_e = kappa (z_e - 1), x = ln r, z = |tau|_rms / (tau0 r^{b-1} l^{(b-1)/2});
    prune r < prune * max r. Returns rest point (radii, active set) and history."""
    g.set_sigma(sigma)
    r = r0.copy(); keep = np.arange(g.m); hist = []
    for it in range(max_steps):
        w = np.zeros(g.m); w[keep] = r[keep] ** 4 / g.L[keep]
        tau, D, f = rms_shear(g, w, keep, sigma)
        tset = tau0 * r[keep] ** (b - 1) * g.L[keep] ** ((b - 1) / 2)
        z = tau / tset
        dx = np.clip(kappa * (z - 1), -delta, delta)
        r[keep] *= np.exp(dx)
        cut = r[keep] < prune * r[keep].max()
        if cut.any():
            trial = keep[~cut]
            if g.support_ok(trial): keep = trial
            else:                                   # remove one at a time; clamp the ones that must stay
                for k in keep[cut]:
                    trial = keep[keep != k]
                    if g.support_ok(trial): keep = trial
                    else: r[k] = max(r[k], prune * r[keep].max())
                cut = r[keep] < prune * r[keep].max()
            keep = g.src_component(keep)
        hist.append((D, len(keep), np.abs(z - 1).max()))
        if verbose and it % 200 == 0: print(f'    it {it} D={D:.5g} active={len(keep)} max|z-1|={hist[-1][2]:.2e}')
        if hist[-1][2] < tol and not cut.any(): break
    w = np.zeros(g.m); w[keep] = r[keep] ** 4 / g.L[keep]
    return dict(r=r, keep=keep, w=w, n_active=len(keep), beta=g.cycle_rank(keep),
                n_sin=int(np.sum(g.kind[keep] == 'sin')), steps=it + 1, converged=hist[-1][2] < tol, hist=hist)
