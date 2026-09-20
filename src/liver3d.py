"""3D liver: hemisphere of hexagonal-prism lobules (Lorente's shape), hilum source at the
centre of the flat face, inlet/outlet trees by recursive ternary clustering in 3D,
lobule = prism (triangular lattice in xy, stacked in z), portal triads on the prism edges
(shared by 3 lobules), sinusoids triad -> central vein within each layer.
Fast incremental evaluation of the exact tree optimum when one lobule changes its feeding triad."""
import numpy as np, networkx as nx
import liver_canopy as lc

def hex_prisms(R=4.0, s=1.0):
    Lh = s / np.sqrt(3.0); dz = s
    centres, corners, pairs = [], {}, []
    zs = np.arange(dz / 2, R, dz)
    for z in zs:
        for j in range(-int(R / (s * 0.866)) - 2, int(R / (s * 0.866)) + 3):
            for i in range(-int(R / s) - 2, int(R / s) + 3):
                x = s * (i + 0.5 * (j % 2)); y = s * 0.866025 * j
                if x * x + y * y + z * z <= R * R: centres.append((x, y, z))
    centres = np.array(centres)
    for c, (x, y, z) in enumerate(centres):
        for k in range(6):
            th = np.pi / 6 + k * np.pi / 3
            p = (round(x + Lh * np.cos(th), 6), round(y + Lh * np.sin(th), 6), round(z, 6))
            if p not in corners: corners[p] = len(corners)
            pairs.append((corners[p], c))
    return centres, np.array(list(corners.keys())), pairs, Lh

def build_liver3d(R=4.0, s=1.0, seed=0):
    centres, corners, pairs, Lh = hex_prisms(R, s)
    origin = np.zeros(3)
    nin, ein, din, lid_in = lc.ternary_tree(corners, origin, seed=seed)
    nout, eout, dout, lid_out = lc.ternary_tree(centres, origin, seed=seed + 1)
    E = [(p, c, 'in', d) for (p, c), d in zip(ein, din)]
    xy = list(map(tuple, nin)); base = len(nin)
    E += [(base + p, base + c, 'out', d) for (p, c), d in zip(eout, dout)]
    xy += list(map(tuple, nout))
    tri = {k: lid_in[k] for k in range(len(corners))}; cen = {k: base + lid_out[k] for k in range(len(centres))}
    E += [(tri[q], cen[c], 'sin', 0) for (q, c) in pairs]
    xy = np.array(xy)
    L = np.array([np.linalg.norm(xy[i] - xy[j]) if kind != 'sin' else Lh for i, j, kind, _ in E]); L = np.maximum(L, 1e-3)
    return dict(xy=xy, E=[(i, j) for i, j, _, _ in E], kind=np.array([k for _, _, k, _ in E]),
                depth=np.array([d for _, _, _, d in E]), L=L, src=0, sink=base, n=len(xy), m=len(E),
                n_lob=len(centres), n_tri=len(corners), Lh=Lh)

class FastTree:
    """inlet problem: tree edges fixed, one sinusoid per lobule chosen. S is updated by path deltas."""
    def __init__(self, g, b, sigma=0.0):
        self.g, self.b, self.sigma = g, b, sigma
        self.alpha = b / 2; self.q = 1 / (1 + self.alpha); self.p = self.alpha * self.q
        self.tree = np.flatnonzero(g.kind != 'sin'); self.sin = np.flatnonzero(g.kind == 'sin')
        T = nx.Graph(); T.add_edges_from((g.E[k][0], g.E[k][1], {'k': k}) for k in self.tree)
        parent = dict(nx.dfs_predecessors(T, g.src))
        self.path = {}                                        # triad node -> array of tree edge ids to root
        for k in self.sin:
            t = g.E[k][0]
            if t in self.path: continue
            v, pth = t, []
            while v in parent:
                pth.append(T[v][parent[v]]['k']); v = parent[v]
            self.path[t] = np.array(pth)
        self.aq = g.a(b) ** self.q
        self.by_centre = {}
        for k in self.sin: self.by_centre.setdefault(g.E[k][1], []).append(k)
        self.centres = sorted(self.by_centre)
    def term(self, k_edges, kcount):
        h = kcount ** 2 + kcount * self.sigma ** 2
        return self.aq[k_edges] * h ** self.p
    def solve(self, seed=0, sweeps=10, init=None):
        g = self.g; rng = np.random.default_rng(seed)
        choice = {c: (init[c] if init else rng.choice(self.by_centre[c])) for c in self.centres}
        k = np.zeros(g.m)
        for c, e in choice.items(): k[self.path[g.E[e][0]]] += 1; k[e] = 1
        term = np.zeros(g.m); act = np.concatenate([self.tree, np.array(list(choice.values()))])
        term[act] = self.term(act, k[act]); S = term.sum()
        for _ in range(sweeps):
            moved = False
            for c in rng.permutation(self.centres):
                cur = choice[c]; pc = self.path[g.E[cur][0]]
                for e in self.by_centre[c]:
                    if e == cur: continue
                    pe = self.path[g.E[e][0]]
                    idx = np.unique(np.concatenate([pc, pe]))
                    knew = k[idx].copy(); knew[np.isin(idx, pc)] -= 1; knew[np.isin(idx, pe)] += 1
                    dS = self.term(idx, knew).sum() - term[idx].sum() + self.aq[e] * (1 + self.sigma ** 2) ** self.p - term[cur]
                    if dS < -1e-12 * S:
                        k[idx] = knew; term[idx] = self.term(idx, knew); term[cur] = 0; k[cur] = 0
                        k[e] = 1; term[e] = self.aq[e] * (1 + self.sigma ** 2) ** self.p
                        S += dS; choice[c] = cur = e; pc = pe; moved = True
            if not moved: break
        keep = np.flatnonzero(k > 0)
        D = S ** (1 + 1 / self.alpha)                        # C0 = 1
        return dict(D=D, S=S, keep=keep, n_active=len(keep), choice=choice, k=k)

if __name__ == '__main__':
    import sys, time, pandas as pd, os
    R = float(sys.argv[1]); out = sys.argv[2] if len(sys.argv) > 2 else 'results'
    t0 = time.time(); G = build_liver3d(R); gin, gout = lc.split_liver(G)
    print(f'R={R}: lobules {G["n_lob"]} triads {G["n_tri"]} nodes {G["n"]} edges {G["m"]} build {time.time()-t0:.0f}s', flush=True)
    rows = []
    for b in (1.0, 0.75, 2 / 3, 0.5):
        t = time.time(); ft = FastTree(gin, b)
        sols = [ft.solve(seed=s) for s in range(3)]; best = min(sols, key=lambda s: s['D'])
        spread = max(s['D'] for s in sols) / best['D'] - 1
        outl = lc.tree_allocation(gout, np.arange(gout.m), b)
        # check closed form against the Laplacian on one case
        chk = np.nan
        if gin.m < 6000:
            sol = lc.tree_allocation(gin, best['keep'], b); chk = gin.dissipation(sol['w'], best['keep']) / best['D'] - 1
        rows.append(dict(R=R, N=gin.N, b=b, D_in=best['D'], D_out=outl['D'], D_tot_perF2=(best['D'] + outl['D']) / gin.N ** 2,
                         active_in=best['n_active'], spread=spread, check=chk, t=time.time() - t))
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    os.makedirs(out, exist_ok=True); pd.DataFrame(rows).to_csv(f'{out}/graph3d_R{R}.tsv', sep='\t', index=False)
