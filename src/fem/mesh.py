"""Structured triangular mesh of one hexagonal lobule.

The domain is the hexagon of circumradius R minus the central vein disc of radius r_cv.
Nodes live on a polar (s, theta) lattice mapped onto the hexagon,

    r(i, j) = r_cv + s_i * (R_hex(theta_j) - r_cv),      s_i = (i/nr)^grade,

so the mesh is valid by construction (no Delaunay, no degenerate elements) and the angular
index is exactly periodic. nt is forced to a multiple of 6, which puts a node on every
corner: the portal triads are mesh nodes, not interpolated.

R_hex(theta) = a / cos((theta mod pi/3) - pi/6) with apothem a = R cos(pi/6): corners sit at
theta = 0, pi/3, ..., i.e. the six portal triads.
"""
import numpy as np


def hex_radius(theta, R):
    """Distance from the centre to the hexagon boundary along theta (corners at theta=0, pi/3...)."""
    a = R * np.cos(np.pi / 6.0)
    return a / np.cos(np.mod(theta, np.pi / 3.0) - np.pi / 6.0)


class Mesh(object):
    """P1 triangular mesh with the bookkeeping the solvers need.

    Attributes
    ----------
    xy    : (N, 2) node coordinates [m]
    tri   : (M, 3) node indices, counter-clockwise
    area  : (M,) element areas [m^2]
    grad  : (M, 3, 2) gradients of the three P1 basis functions on each element
    r, th : (N,) polar coordinates of the nodes
    zeta  : (N,) 0 at the portal boundary, 1 at the central vein (the zonation coordinate)
    inner : node ids on the central vein
    outer : node ids on the hexagon boundary
    portal: node ids inside a portal triad window
    """

    def __init__(self, geom):
        nr, nt = int(geom.nr), int(geom.nt)
        nt = nt + (-nt) % 6                                  # force a multiple of 6
        self.nr, self.nt, self.geom = nr, nt, geom
        R, r_cv = geom.R, geom.r_cv

        th = 2.0 * np.pi * np.arange(nt) / nt
        s = (np.arange(nr + 1) / float(nr)) ** geom.grade
        Rh = hex_radius(th, R)                               # (nt,)
        r = r_cv + s[:, None] * (Rh[None, :] - r_cv)         # (nr+1, nt)
        TH = np.broadcast_to(th[None, :], r.shape)

        self.r = r.ravel()
        self.th = TH.ravel()
        self.xy = np.column_stack([self.r * np.cos(self.th), self.r * np.sin(self.th)])
        self.n = self.xy.shape[0]
        self.zeta = np.clip(1.0 - (self.r - r_cv) / (np.broadcast_to(Rh[None, :], r.shape).ravel() - r_cv),
                            0.0, 1.0)

        idx = lambda i, j: i * nt + (j % nt)
        tri = []
        for i in range(nr):
            for j in range(nt):
                a_, b_, c_, d_ = idx(i, j), idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)
                tri.append((a_, b_, c_))
                tri.append((a_, c_, d_))
        self.tri = np.asarray(tri, dtype=np.int64)
        self._geometry()

        self.inner = np.arange(nt)                            # i = 0
        self.outer = np.arange(nr * nt, (nr + 1) * nt)        # i = nr
        self._rings()
        self.dth_corner = np.abs(((self.th[self.outer] + np.pi / 6.0) % (np.pi / 3.0)) - np.pi / 6.0)
        self.set_inlet(getattr(geom, 'inlet', 'perimeter'))

    def set_inlet(self, mode):
        """'perimeter': the whole hexagon boundary is fed by the interlobular venules
        (the usual porous-lobule boundary condition, and the one without stagnant corners).
        'triads': only six windows of half-width geom.triad_arc at the corners, i.e. the
        point-source idealisation that the LiPNet network graph uses."""
        if mode == 'perimeter':
            self.portal = self.outer.copy()
            self.outer_wall = np.array([], dtype=np.int64)
        elif mode == 'triads':
            self.portal = self.outer[self.dth_corner <= self.geom.triad_arc]
            self.outer_wall = self.outer[self.dth_corner > self.geom.triad_arc]
            if self.portal.size == 0:
                raise ValueError('triad_arc too small for this nt: no portal node selected')
        else:
            raise ValueError('inlet must be "perimeter" or "triads"')
        self.inlet_mode = mode
        return self

    def _geometry(self):
        p = self.xy[self.tri]                                 # (M, 3, 2)
        v1 = p[:, 1] - p[:, 0]
        v2 = p[:, 2] - p[:, 0]
        cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
        flip = cross < 0
        if np.any(flip):                                      # orient counter-clockwise
            self.tri[flip] = self.tri[flip][:, [0, 2, 1]]
            p = self.xy[self.tri]
            v1 = p[:, 1] - p[:, 0]
            v2 = p[:, 2] - p[:, 0]
            cross = v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]
        self.area = 0.5 * cross
        if np.any(self.area <= 0):
            raise RuntimeError('degenerate element in the lobule mesh')
        x, y = p[:, :, 0], p[:, :, 1]
        b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1)
        c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1)
        self.grad = np.stack([b, c], axis=2) / (2.0 * self.area)[:, None, None]
        self.centroid = p.mean(axis=1)
        self.h = np.sqrt(2.0 * self.area)                     # element size for SUPG
        self.lumped = np.bincount(self.tri.ravel(),
                                  weights=np.repeat(self.area / 3.0, 3), minlength=self.n)

    def _rings(self):
        """Boundary edges with their adjacent element, length and outward unit normal.

        Both rings are stored as (nodeA, nodeB, element, length, nx, ny). The element is
        known from the construction of the mesh, so no search is needed; the normal points
        out of the lobule (away from the centre on the hexagon, into the lumen on the
        central vein).
        """
        nr, nt = self.nr, self.nt
        idx = lambda i, j: i * nt + (j % nt)
        out, inn = [], []
        for j in range(nt):
            a, b = idx(nr, j), idx(nr, j + 1)
            out.append((a, b, 2 * ((nr - 1) * nt + j) + 1, +1.0))
            a2, b2 = idx(0, j), idx(0, j + 1)
            inn.append((a2, b2, 2 * (0 * nt + j) + 0, -1.0))

        def pack(rows):
            A = np.array([r[0] for r in rows]); B = np.array([r[1] for r in rows])
            E = np.array([r[2] for r in rows]); sgn = np.array([r[3] for r in rows])
            d = self.xy[B] - self.xy[A]
            L = np.hypot(d[:, 0], d[:, 1])
            nvec = np.column_stack([d[:, 1], -d[:, 0]]) / L[:, None] * sgn[:, None]
            return dict(a=A, b=B, elem=E, length=L, normal=nvec)

        self.ring_outer = pack(out)
        self.ring_inner = pack(inn)

    def portal_edges(self):
        """Mask of outer-ring edges whose two nodes are both inlet nodes."""
        is_in = np.zeros(self.n, dtype=bool)
        is_in[self.portal] = True
        r = self.ring_outer
        return is_in[r['a']] & is_in[r['b']]

    # -- helpers ----------------------------------------------------------
    def elem_mean(self, c):
        """Nodal field -> element average."""
        return c[self.tri].mean(axis=1)

    def elem_to_node(self, ce):
        """Element field -> nodal field, area-weighted."""
        w = np.bincount(self.tri.ravel(), weights=np.repeat(self.area, 3), minlength=self.n)
        v = np.bincount(self.tri.ravel(), weights=np.repeat(ce * self.area, 3), minlength=self.n)
        return v / np.maximum(w, 1e-300)

    def integrate(self, c):
        """Integral of a nodal field over the lobule cross-section [field * m^2]."""
        return float(np.sum(self.area * self.elem_mean(c)))

    @property
    def cross_area(self):
        return float(np.sum(self.area))

    def volume(self):
        return self.cross_area * self.geom.H

    def radial_profile(self, c, nbin=40):
        """Field averaged over theta against the zonation coordinate zeta."""
        edges = np.linspace(0.0, 1.0, nbin + 1)
        which = np.clip(np.digitize(self.zeta, edges) - 1, 0, nbin - 1)
        num = np.bincount(which, weights=c * self.lumped, minlength=nbin)
        den = np.bincount(which, weights=self.lumped, minlength=nbin)
        ok = den > 0
        return 0.5 * (edges[:-1] + edges[1:])[ok], (num[ok] / den[ok])


def build(geom):
    return Mesh(geom)
