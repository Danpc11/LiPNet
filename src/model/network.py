"""Liver as a perfusion network: hilum -> portal tree -> hexagonal lobules -> hepatic venous tree -> outlet.

Formulation:
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

"""
import numpy as np, networkx as nx
from scipy.linalg import solve

LN3 = np.log(3.0)

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
