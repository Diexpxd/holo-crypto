"""
holocrypto_engine_v1 - the Bar-Natan / van der Veen invariant Θ = (Δ, θ).

Reference: D. Bar-Natan, R. van der Veen, "A Fast, Strong, Topologically Meaningful
and Fun Knot Invariant", arXiv:2509.18456 (v4, May 2026), Section 2 ("The Main Theorem").
Formulas (1)-(6) are copied from the paper's TeX source. The knot data and the
reference value of θ(trefoil) come from van der Veen's Sage implementation
(rolandvdv.nl/Theta); they are used ONLY for verification, not for the computation.

Input: an "upright" diagram of a long knot with n crossings:
  * edges numbered 1..2n+1 along the knot;
  * each crossing is (s, i, j): sign s=±1, incoming over-edge i, incoming under-edge j;
  * rotation numbers φ_k of each edge (needed for the normalised Δ and for θ).

Mathematics implemented:
  A  = I - Σ_c ( T^s E_{i,i+1} + (1-T^s) E_{i,j+1} + E_{j,j+1} )          ((2n+1)² matrix)
  Δ  = T^{(-Σφ - w)/2} det A                                               (Alexander)
  G  = A⁻¹ ;  g_{ν,αβ} = G with T→T_ν ;  T_3 = T_1 T_2                     (Green)
  θ0 = Σ_c F1(c) + Σ_{c0,c1} F2(c0,c1) + Σ_k F3(k)
  θ  = Δ(T1) Δ(T2) Δ(T3) · θ0                                             (Laurent polynomial)

Hash-style output (7d): SHA-256(domain || q || T1 || T2 || θ(T1,T2)); random keys
reject knots with θ = 0.

Usage:  python holocrypto_engine_v1.py [--base] [--input] [--hardening] [--bench]
"""
from __future__ import annotations

import hashlib
import itertools
import random
import sys
import time
from dataclasses import dataclass
from fractions import Fraction


# ==========================================================================
# 1. Rings
# ==========================================================================
class Fp:
    """Finite field F_q (q prime). Same interface as Fraction for the θ computation."""
    __slots__ = ("v", "q")

    def __init__(self, v: int, q: int) -> None:
        self.v, self.q = v % q, q

    def _w(self, o) -> "Fp":
        return o if isinstance(o, Fp) else Fp(int(o), self.q)

    def __add__(self, o): return Fp(self.v + self._w(o).v, self.q)
    def __sub__(self, o): return Fp(self.v - self._w(o).v, self.q)
    def __mul__(self, o): return Fp(self.v * self._w(o).v, self.q)
    def __radd__(self, o): return self + o
    def __rmul__(self, o): return self * o
    def __rsub__(self, o): return self._w(o) - self
    def __neg__(self): return Fp(-self.v, self.q)
    def __truediv__(self, o): return self * self._w(o).inv()
    def __rtruediv__(self, o): return self._w(o) * self.inv()
    def inv(self): return Fp(pow(self.v, -1, self.q), self.q)
    def __pow__(self, k: int): return Fp(pow(self.v, k, self.q), self.q)
    def __eq__(self, o): return self.v == self._w(o).v
    def __hash__(self): return hash((self.v, self.q))
    def __repr__(self): return f"{self.v} (mod {self.q})"


class LaurentPoly:
    """Laurent polynomial in T with integer coefficients: {exponent: coefficient}."""
    __slots__ = ("c",)

    def __init__(self, c: dict[int, int] | None = None) -> None:
        self.c = {e: v for e, v in (c or {}).items() if v}

    @classmethod
    def mono(cls, e: int, v: int = 1) -> "LaurentPoly":
        return cls({e: v})

    def __add__(self, o: "LaurentPoly") -> "LaurentPoly":
        r = dict(self.c)
        for e, v in o.c.items():
            r[e] = r.get(e, 0) + v
        return LaurentPoly(r)

    def __neg__(self) -> "LaurentPoly":
        return LaurentPoly({e: -v for e, v in self.c.items()})

    def __sub__(self, o: "LaurentPoly") -> "LaurentPoly":
        return self + (-o)

    def __mul__(self, o: "LaurentPoly") -> "LaurentPoly":
        r: dict[int, int] = {}
        for e1, v1 in self.c.items():
            for e2, v2 in o.c.items():
                r[e1 + e2] = r.get(e1 + e2, 0) + v1 * v2
        return LaurentPoly(r)

    def shift(self, k: int) -> "LaurentPoly":
        return LaurentPoly({e + k: v for e, v in self.c.items()})

    def exact_div(self, d: "LaurentPoly") -> "LaurentPoly":
        """Exact division by d (in Z[T,T⁻¹]); AssertionError if it is not exact."""
        assert d.c, "division by zero"
        ed = max(d.c)
        lead = d.c[ed]
        rem, q = dict(self.c), {}
        while rem:
            er = max(rem)
            cq, r = divmod(rem[er], lead)
            assert r == 0, "inexact division"
            q[er - ed] = cq
            for e, v in d.c.items():
                rem[e + er - ed] = rem.get(e + er - ed, 0) - cq * v
                if rem[e + er - ed] == 0:
                    del rem[e + er - ed]
        return LaurentPoly(q)

    def __call__(self, t):
        return sum((v * t ** e for e, v in self.c.items()), t * 0)

    def __eq__(self, o): return isinstance(o, LaurentPoly) and self.c == o.c
    def __bool__(self): return bool(self.c)

    def __repr__(self) -> str:
        if not self.c:
            return "0"
        return " ".join(f"{v:+d}·T^{e}" for e, v in sorted(self.c.items()))


# ==========================================================================
# 2. Exact linear algebra
# ==========================================================================
def det_bareiss(M: list[list[LaurentPoly]]) -> LaurentPoly:
    """Fraction-free determinant (Bareiss) over Z[T,T⁻¹]. O(n³) ring operations."""
    n = len(M)
    M = [row[:] for row in M]
    sign, prev = 1, LaurentPoly.mono(0)
    for k in range(n - 1):
        if not M[k][k]:
            piv = next((r for r in range(k + 1, n) if M[r][k]), None)
            if piv is None:
                return LaurentPoly()
            M[k], M[piv] = M[piv], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]).exact_div(prev)
        prev = M[k][k]
    d = M[n - 1][n - 1]
    return d if sign == 1 else -d


def det_inv_field(M: list[list]) -> tuple[object, list[list]]:
    """(det, inverse) by Gauss-Jordan over a field (Fraction or Fp)."""
    n = len(M)
    zero = M[0][0] - M[0][0]
    one = zero + 1
    a = [row[:] + [one if i == j else zero for j in range(n)] for i, row in enumerate(M)]
    det = one
    for k in range(n):
        p = next((r for r in range(k, n) if not (a[r][k] == zero)), None)
        if p is None:
            raise ZeroDivisionError("singular matrix at this point (T is a root of Δ)")
        if p != k:
            a[k], a[p] = a[p], a[k]
            det = -det
        pv = a[k][k]
        det = det * pv
        a[k] = [x / pv for x in a[k]]
        for r in range(n):
            if r != k and not (a[r][k] == zero):
                f = a[r][k]
                a[r] = [x - f * y for x, y in zip(a[r], a[k])]
    return det, [row[n:] for row in a]


# ==========================================================================
# 3. Knot diagram
# ==========================================================================
@dataclass(frozen=True)
class KnotDiagram:
    """Upright diagram of a long knot (edges 1..2n+1)."""
    crossings: tuple[tuple[int, int, int], ...]     # (s, i, j), edges indexed from 1
    rotations: tuple[int, ...]                       # φ_1..φ_{2n+1}
    name: str = ""

    def __post_init__(self) -> None:
        n = len(self.crossings)
        assert len(self.rotations) == 2 * n + 1, "2n+1 rotation numbers expected"
        steps = sorted(p for _, i, j in self.crossings for p in (i, j))
        assert steps == list(range(1, 2 * n + 1)), \
            "each pass 1..2n must be the over- or under-pass of exactly one crossing"
        assert all(s in (1, -1) for s, _, _ in self.crossings)

    @property
    def n(self) -> int:
        return len(self.crossings)

    @property
    def writhe(self) -> int:
        return sum(s for s, _, _ in self.crossings)

    @classmethod
    def from_sage(cls, rows: list[list[int]], rot: list[int], name: str = "") -> "KnotDiagram":
        """van der Veen format: rows [s, over, under] with edges indexed from 0."""
        return cls(tuple((s, i + 1, j + 1) for s, i, j in rows), tuple(rot), name)

    def matrix_A(self, T, one):
        """A(T) as a list of lists, with entries in the ring of T (Fraction/Fp/LaurentPoly)."""
        m = 2 * self.n + 1
        zero = one - one
        A = [[one if r == c else zero for c in range(m)] for r in range(m)]
        for s, i, j in self.crossings:
            Ts = T ** s
            A[i - 1][i] = A[i - 1][i] - Ts
            A[i - 1][j] = A[i - 1][j] + Ts - one
            A[j - 1][j] = A[j - 1][j] - one
        return A

    def matrix_A_laurent(self) -> list[list[LaurentPoly]]:
        m = 2 * self.n + 1
        one = LaurentPoly.mono(0)
        A = [[one if r == c else LaurentPoly() for c in range(m)] for r in range(m)]
        for s, i, j in self.crossings:
            A[i - 1][i] = A[i - 1][i] - LaurentPoly.mono(s)
            A[i - 1][j] = A[i - 1][j] + LaurentPoly.mono(s) - one
            A[j - 1][j] = A[j - 1][j] - one
        return A

    @property
    def delta_exponent(self) -> int:
        e, r = divmod(-sum(self.rotations) - self.writhe, 2)
        assert r == 0, "Σφ + w must be even"
        return e


# ==========================================================================
# 4. Δ: Alexander polynomial
# ==========================================================================
def alexander(d: KnotDiagram) -> LaurentPoly:
    """Δ(K) = T^{(-φ(D)-w(D))/2} · det A   (equation (2) of the paper)."""
    if d.n == 0:
        return LaurentPoly.mono(0)
    return det_bareiss(d.matrix_A_laurent()).shift(d.delta_exponent)


# ==========================================================================
# 5. θ: evaluation at a point (T1, T2) of a field
# ==========================================================================
def theta_eval(d: KnotDiagram, T1, T2):
    """θ(D) evaluated at (T1, T2) = Δ(T1)Δ(T2)Δ(T3)·θ0, with T3 = T1·T2. Exact arithmetic."""
    one = T1 ** 0
    half = one / (one + one)
    T3 = T1 * T2
    Ts_ = {1: T1, 2: T2, 3: T3}
    G, Del = {}, 1
    for nu, T in Ts_.items():
        det, inv = det_inv_field(d.matrix_A(T, one))
        G[nu] = inv
        Del = Del * (T ** d.delta_exponent) * det
    if T2 == one:
        raise ZeroDivisionError("T2 = 1 is an apparent pole; evaluate at another point")

    def g(nu: int, a: int, b: int):
        return G[nu][a - 1][b - 1]

    def F1(c):
        s, i, j = c
        t1s, t2s, t3s = T1 ** s, T2 ** s, T3 ** s
        poly = (half - g(3, i, i) + t2s * g(1, i, i) * g(2, j, i) - t2s * g(3, j, j) * g(2, j, i)
                - (t2s - 1) * g(3, i, i) * g(2, j, i) + (t3s - 1) * g(2, j, i) * g(3, j, i)
                - g(1, i, i) * g(2, j, j) + 2 * g(3, i, i) * g(2, j, j)
                + g(1, i, i) * g(3, j, j) - g(2, i, i) * g(3, j, j))
        frac = ((t1s - 1) * t2s * (g(3, j, j) * g(1, j, i) - g(2, j, j) * g(1, j, i)
                                   + t2s * g(1, j, i) * g(2, j, i))
                + (t3s - 1) * g(3, j, i) * (1 - t2s * g(1, i, i) + g(2, i, j) + (t2s - 2) * g(2, j, j)
                                            - (t1s - 1) * (t2s + 1) * g(1, j, i)))
        return s * poly + s * frac / (t2s - 1)

    def F2(c0, c1):
        s0, i0, j0 = c0
        s1, i1, j1 = c1
        t2s0, t2s1 = T2 ** s0, T2 ** s1
        pref = s1 * (T1 ** s0 - 1) * (T3 ** s1 - 1) * g(1, j1, i0) * g(3, j0, i1) / (t2s1 - 1)
        return pref * (t2s0 * g(2, i1, i0) + g(2, j1, j0) - t2s0 * g(2, j1, i0) - g(2, i1, j0))

    def F3(k: int):
        return (g(3, k, k) - half) * d.rotations[k - 1]

    th0 = sum((F1(c) for c in d.crossings), one * 0)
    th0 = th0 + sum((F2(c0, c1) for c0 in d.crossings for c1 in d.crossings), one * 0)
    th0 = th0 + sum((F3(k) for k in range(1, 2 * d.n + 2)), one * 0)
    return Del * th0


# ==========================================================================
# 6. θ: full polynomial by (exact) interpolation
# ==========================================================================
def _solve(M: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    _, inv = det_inv_field(M)
    return [sum((inv[i][j] * b[j] for j in range(len(b))), Fraction(0)) for i in range(len(b))]


def theta_polynomial(d: KnotDiagram, D: int | None = None, check: bool = True
                     ) -> dict[tuple[int, int], int]:
    """θ(D) ∈ Z[T1^±,T2^±] by exact tensor interpolation; {(a,b): coef} of T1^a T2^b.

    Evaluates on a (2D+1)² grid of rationals and solves two Vandermonde systems.
    D bounds |exponent|; checked at extra points (raises an error if it is not enough).
    """
    D = D if D is not None else d.n + 2
    m = 2 * D + 1
    xs = [Fraction(k + 2) for k in range(m)]            # 2, 3, ..., 2D+2  (T2 ≠ 1)
    ys = [Fraction(k) + Fraction(1, 2) + 2 for k in range(m)]
    V = lambda pts: [[p ** e for e in range(m)] for p in pts]
    vals = [[theta_eval(d, x, y) * (x * y) ** D for y in ys] for x in xs]    # polynomial P(x,y)
    Vx_inv = det_inv_field(V(xs))[1]
    Vy_inv = det_inv_field(V(ys))[1]
    # coef[a][b] = Σ_i,k Vx_inv[a][i] vals[i][k] Vy_inv[b][k]
    tmp = [[sum((Vx_inv[a][i] * vals[i][k] for i in range(m)), Fraction(0)) for k in range(m)]
           for a in range(m)]
    coef = [[sum((tmp[a][k] * Vy_inv[b][k] for k in range(m)), Fraction(0)) for b in range(m)]
            for a in range(m)]
    out: dict[tuple[int, int], int] = {}
    for a in range(m):
        for b in range(m):
            c = coef[a][b]
            if c:
                if c.denominator != 1:
                    raise ValueError("non-integer coefficient: the bound D is insufficient or there is a bug")
                out[(a - D, b - D)] = int(c)
    if check:
        rng = random.Random(1)
        for _ in range(3):
            x, y = Fraction(rng.randint(5, 40), 7), Fraction(rng.randint(5, 40), 11)
            if y == 1:
                continue
            expected = theta_eval(d, x, y)
            obtained = sum((Fraction(v) * x ** a * y ** b for (a, b), v in out.items()), Fraction(0))
            if expected != obtained:
                raise ValueError("the interpolation does not reproduce θ at a new point: increase D")
    return out


def format_poly(poly: dict[tuple[int, int], int]) -> str:
    return " ".join(f"{v:+d}·T1^{a}T2^{b}" for (a, b), v in sorted(poly.items()))


# ==========================================================================
# 7. Verification data (van der Veen, knots.sage; edges indexed from 0)
# ==========================================================================
KNOTS = {
    "0_1": ([], [0]),
    "3_1": ([[-1, 3, 0], [-1, 5, 2], [-1, 1, 4]], [0, 0, 0, -1, 0, 0, 0]),
    "4_1": ([[1, 0, 3], [1, 4, 7], [-1, 2, 5], [-1, 6, 1]], [0, 0, 0, -1, 0, 0, -1, 0, 0]),
    "5_1": ([[-1, 5, 0], [-1, 7, 2], [-1, 9, 4], [-1, 1, 6], [-1, 3, 8]],
            [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0]),
    "5_2": ([[-1, 3, 0], [-1, 7, 2], [-1, 9, 4], [-1, 5, 8], [-1, 1, 6]],
            [0, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0]),
    "6_1": ([[-1, 3, 0], [-1, 9, 6], [1, 7, 2], [1, 1, 8], [-1, 11, 4], [-1, 5, 10]],
            [0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 1, -1, 0]),
    "6_2": ([[-1, 3, 0], [-1, 9, 4], [1, 7, 2], [1, 1, 8], [-1, 11, 6], [-1, 5, 10]],
            [0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0]),
    "6_3": ([[1, 0, 3], [1, 2, 7], [-1, 8, 11], [-1, 4, 9], [-1, 10, 5], [1, 6, 1]],
            [0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0]),
    "7_1": ([[-1, 7, 0], [-1, 9, 2], [-1, 11, 4], [-1, 13, 6], [-1, 1, 8], [-1, 3, 10], [-1, 5, 12]],
            [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0]),
    "8_19": ([[1, 0, 3], [1, 2, 7], [1, 13, 8], [1, 11, 4], [1, 5, 12], [1, 15, 10], [1, 9, 14], [1, 6, 1]],
             [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]),
    "8_20": ([[1, 0, 3], [1, 2, 7], [-1, 11, 4], [-1, 15, 12], [-1, 13, 8], [-1, 9, 14], [-1, 5, 10], [1, 6, 1]],
             [0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, -1, 0, -1, 0, 0, 0]),
}


def knot(name: str) -> KnotDiagram:
    rows, rot = KNOTS[name]
    return KnotDiagram.from_sage(rows, rot, name)


def paper_trefoil() -> KnotDiagram:
    """Figure 'upright knot diagram' of the paper: X = {(1,1,4),(1,5,2),(1,3,6)}, φ_4 = -1."""
    return KnotDiagram(((1, 1, 4), (1, 5, 2), (1, 3, 6)), (0, 0, 0, -1, 0, 0, 0), "3_1 (paper)")


def _L(*terms: tuple[int, int]) -> LaurentPoly:
    return LaurentPoly({e: c for e, c in terms})


# Known Δ (normalised, symmetric Alexander polynomials, Δ(1)=1)
ALEXANDER_KNOWN = {
    "0_1": _L((0, 1)),
    "3_1": _L((-1, 1), (0, -1), (1, 1)),
    "4_1": _L((-1, -1), (0, 3), (1, -1)),
    "5_1": _L((-2, 1), (-1, -1), (0, 1), (1, -1), (2, 1)),
    "5_2": _L((-1, 2), (0, -3), (1, 2)),
    "6_1": _L((-1, -2), (0, 5), (1, -2)),
    "6_2": _L((-2, -1), (-1, 3), (0, -3), (1, 3), (2, -1)),
    "6_3": _L((-2, 1), (-1, -3), (0, 5), (1, -3), (2, 1)),
    "7_1": _L((-3, 1), (-2, -1), (-1, 1), (0, -1), (1, 1), (2, -1), (3, 1)),
    "8_19": _L((-3, 1), (-2, -1), (0, 1), (2, -1), (3, 1)),
    "8_20": _L((-2, 1), (-1, -2), (0, 3), (1, -2), (2, 1)),
}

# θ(3_1) published by van der Veen for Xings=[[1,0,3],[1,4,1],[1,2,5]], rot=[0,0,0,-1,0,0,0]
# (t1^4 t2^4 - t1^4 t2^3 - t1^3 t2^4 + t1^4 t2^2 + t1^2 t2^4 - t1^3 t2 - t1 t2^3
#   + t1^2 + t2^2 - t1 - t2 + 1) / (t1^2 t2^2)
THETA_TREFOIL_REFERENCE = {(2, 2): 1, (2, 1): -1, (1, 2): -1, (2, 0): 1, (0, 2): 1,
                           (1, -1): -1, (-1, 1): -1, (0, -2): 1, (-2, 0): 1,
                           (-1, -2): -1, (-2, -1): -1, (-2, -2): 1}


# ==========================================================================
# 7b. Input: PD / DT / braids  ->  upright diagram with rotation numbers
# ==========================================================================
# Theory (derived here; the paper only says "rotate each crossing until it is upright,
# the different choices differ by instances of the Sw relation").
#
#  * Upright crossing: both strands point upwards. Rays from the centre (in degrees):
#      s=+1: under-in 315, over-out 45, under-out 135, over-in 225   (PD order i,j,k,l)
#      s=-1: under-in 225, over-in 315, under-out 45, over-out 135
#    The 4 slots advance +90° counter-clockwise, and each corner of a face measures 90°.
#  * Turning of an edge e (from tangent angle θ_s to θ_e): τ_e = θ_e - θ_s + 360·φ_e, with φ_e
#    an integer = the paper's rotation number (cups +1, caps -1).
#  * Each face (boundary traversed with the face on the left) sums to 360·ε, ε=+1 bounded, -1 outer:
#        Σ ± τ_e + 90·(number of corners) = 360·ε       (± depending on traversal direction)
#    => an INTEGER flow system on the dual graph: Σ ± φ_e = c_f. It is solved with a spanning
#    tree of the dual (the remaining edges set to 0); the other solutions differ by 2π turns at a
#    crossing, i.e. by instances of Sw.
#  * Long knot: the edge e0 is cut; the face on its left is the outer one (the closure at
#    infinity goes around on the right, φ = -1): φ_1 = 0, φ_{2n+1} = φ_{e0} + 1.
def _rays(s: int) -> tuple[int, int, int, int]:
    return (315, 45, 135, 225) if s > 0 else (225, 315, 45, 135)


class PlanarMap:
    """4-valent planar map of a PD (KnotTheory convention X[i,j,k,l]; edges 1..2n)."""

    def __init__(self, pd) -> None:
        self.pd = [tuple(x) for x in pd]
        self.n = n = len(self.pd)
        self.E = E = 2 * n
        assert E >= 4, "at least 2 crossings are needed (the sign of a crossing with 2 edges is ambiguous)"
        self.sign: list[int] = []
        self.tail: dict[int, tuple[int, int]] = {}      # edge -> (crossing, slot) it LEAVES through
        self.head: dict[int, tuple[int, int]] = {}    # edge -> (crossing, slot) it ENTERS through
        for c, (i, j, k, l) in enumerate(self.pd):
            if (j - l) % E == 1:
                s = 1
            elif (l - j) % E == 1:
                s = -1
            else:
                raise ValueError(f"crossing {c} {self.pd[c]}: j and l are not consecutive")
            assert (k - i) % E == 1, f"crossing {c}: the lower strand must go from i to i+1"
            self.sign.append(s)
            roles = ((True, i), (False, j), (False, k), (True, l)) if s > 0 else \
                    ((True, i), (True, j), (False, k), (False, l))
            for q, (is_head, e) in enumerate(roles):
                d = self.head if is_head else self.tail
                assert e not in d, f"edge {e} appears twice as {'input' if is_head else 'output'}"
                d[e] = (c, q)
        assert set(self.head) == set(self.tail) == set(range(1, E + 1)), "invalid edge labels"
        ray = lambda dart: _rays(self.sign[dart[0]])[dart[1]]
        self.th_s = {e: ray(self.tail[e]) for e in range(1, E + 1)}
        self.th_e = {e: (ray(self.head[e]) + 180) % 360 for e in range(1, E + 1)}
        state_of = {}                       # slot -> state (edge, ±1) when LEAVING through it
        for e in range(1, E + 1):
            state_of[self.tail[e]] = (e, 1)
            state_of[self.head[e]] = (e, -1)
        self.face_of: dict[tuple[int, int], int] = {}
        self.faces: list[list[tuple[int, int]]] = []
        for e0 in range(1, E + 1):
            for d0 in (1, -1):
                if (e0, d0) in self.face_of:
                    continue
                f, orbit, st = len(self.faces), [], (e0, d0)
                while st not in self.face_of:
                    self.face_of[st] = f
                    orbit.append(st)
                    e, d = st
                    c, q = self.head[e] if d > 0 else self.tail[e]     # arrival slot
                    st = state_of[(c, (q - 1) % 4)]                     # turn left
                self.faces.append(orbit)

    @property
    def is_planar(self) -> bool:
        return len(self.faces) == self.n + 2           # V - E + F = 2


def solve_rotations(mp: PlanarMap, e_cut: int, rng: random.Random | None = None,
                        twists: int = 0) -> dict[int, int]:
    """φ_e (integers) for each edge of the closed diagram, with the face to the left of e_cut being the outer one."""
    assert mp.is_planar, "the PD is not planar (sphere): F != n+2"
    rng = rng or random.Random(0)
    outer = mp.face_of[(e_cut, 1)]
    c_f = []
    for f, orbit in enumerate(mp.faces):
        total = sum(d * (mp.th_e[e] - mp.th_s[e]) for e, d in orbit) + 90 * len(orbit)
        num = 360 * (-1 if f == outer else 1) - total
        assert num % 360 == 0, "non-integer constraint: inconsistent diagram or convention"
        c_f.append(num // 360)
    left = {e: mp.face_of[(e, 1)] for e in range(1, mp.E + 1)}
    right = {e: mp.face_of[(e, -1)] for e in range(1, mp.E + 1)}
    adj: list[list[tuple[int, int]]] = [[] for _ in mp.faces]
    for e in range(1, mp.E + 1):
        adj[left[e]].append((e, right[e]))
        adj[right[e]].append((e, left[e]))
    root = rng.randrange(len(mp.faces))
    parent: dict[int, tuple[int, int]] = {}
    order, seen = [root], {root}
    for f in order:                                    # BFS with shuffled neighbours
        nbrs = adj[f][:]
        rng.shuffle(nbrs)
        for e, g in nbrs:
            if g not in seen:
                seen.add(g)
                parent[g] = (e, f)
                order.append(g)
    assert len(seen) == len(mp.faces), "the dual is not connected"
    phi = {e: 0 for e in range(1, mp.E + 1)}
    r = list(c_f)
    for f in reversed(order[1:]):
        e, g = parent[f]
        phi[e] = r[f] if f == left[e] else -r[f]
        r[g] -= phi[e] if g == left[e] else -phi[e]
    assert r[root] == 0, "the face constraints do not sum to 0"
    for _ in range(twists):                             # 2π turns at a crossing (Sw relation)
        c = rng.randrange(mp.n)
        k = rng.choice((-2, -1, 1, 2))
        for q in range(4):
            e = next(e for e in range(1, mp.E + 1) if mp.tail[e] == (c, q) or mp.head[e] == (c, q))
            phi[e] += k if mp.tail[e] == (c, q) else -k
    check_rotations(mp, phi, e_cut)
    return phi


def check_rotations(mp: PlanarMap, phi: dict[int, int], e_cut: int) -> None:
    outer = mp.face_of[(e_cut, 1)]
    for f, orbit in enumerate(mp.faces):
        turn = sum(d * (mp.th_e[e] - mp.th_s[e] + 360 * phi[e]) for e, d in orbit) + 90 * len(orbit)
        assert turn == (-360 if f == outer else 360), f"face {f}: total turning {turn}"


def pd_to_upright(pd, e_cut: int = 1, rng: random.Random | None = None, twists: int = 0,
                 name: str = "") -> KnotDiagram:
    """PD (KnotTheory) -> upright KnotDiagram of a long knot cut at edge e_cut."""
    mp = PlanarMap(pd)
    phi = solve_rotations(mp, e_cut, rng, twists)
    E = mp.E
    lab = lambda e: (e - e_cut) % E + 1
    crossings = []
    for c, (i, j, k, l) in enumerate(mp.pd):
        s = mp.sign[c]
        crossings.append((s, lab(l if s > 0 else j), lab(i)))
    rot = [0] * (E + 1)
    for e in range(1, E + 1):
        rot[lab(e) - 1] = phi[e]
    rot[0], rot[E] = 0, phi[e_cut] + 1
    return KnotDiagram(tuple(crossings), tuple(rot), name)


# --------------------------- DT -> PD (sign search) ---------------------------
def dt_to_pd(dt: list[int], max_n: int = 16) -> list[tuple[int, int, int, int]]:
    """DT code (signed pairs: sign + => the odd pass goes over) -> planar PD.

    The DT does not fix the chirality of each crossing: a sign vector whose map is planar
    (F = n+2) is found by brute force (2^(n-1) candidates, the global mirror image is discarded).
    Returns ONE of the two mirror knots. Only for small n (max_n)."""
    n, E = len(dt), 2 * len(dt)
    if n > max_n:
        raise ValueError(f"DT->PD exhaustive search is limited to n<={max_n}; use PD or braids")
    pairs = []
    for idx, d in enumerate(dt):
        a, b = 2 * idx + 1, abs(d)
        over, under = (a, b) if d > 0 else (b, a)
        pairs.append((under, over))
    assert sorted(p for pr in pairs for p in pr) == list(range(1, E + 1)), "invalid DT"
    for mask in range(1 << (n - 1)):
        signs = [1] + [1 if (mask >> t) & 1 else -1 for t in range(n - 1)]
        pd = []
        for (u, o), s in zip(pairs, signs):
            n_under, n_over = u % E + 1, o % E + 1
            pd.append((u, n_over, n_under, o) if s > 0 else (u, o, n_under, n_over))
        if PlanarMap(pd).is_planar:
            return pd
    raise ValueError("no sign vector gives a planar map: DT not realizable")


# --------------------------- braids -> upright / PD ---------------------------
def _walk_braid(word: list[int], m: int, p: int):
    steps, arcs, x, strands = [], 0, p, 0         # steps: (letter, goes over?, preceding closure arcs)
    while True:
        pos = x
        for c, g in enumerate(word):
            k = abs(g)
            if pos == k:
                steps.append((c, g > 0, arcs)); arcs = 0; pos = k + 1
            elif pos == k + 1:
                steps.append((c, g < 0, arcs)); arcs = 0; pos = k
        strands += 1
        if pos == p:
            break
        arcs += 1
        x = pos
    if strands != m:
        raise ValueError("the braid closure is not a knot (several components)")
    return steps


def braid_to_upright(word: list[int], m: int, p: int = 1, name: str = "") -> KnotDiagram:
    """Closure of the braid (σ_k^{±1} = ±k, m strands) cut at the closing arc of strand p.

    Direct, no planarity needed: the strands go up, each closing arc (on the right) contributes φ = -1."""
    steps = _walk_braid(word, m, p)
    over, under = {}, {}
    for idx, (c, is_over, _) in enumerate(steps, start=1):
        (over if is_over else under)[c] = idx
    crossings = tuple((1 if word[c] > 0 else -1, over[c], under[c]) for c in range(len(word)))
    rot = tuple(-a for _, _, a in steps) + (0,)
    return KnotDiagram(crossings, rot, name)


def braid_to_pd(word: list[int], m: int, p: int = 1) -> list[tuple[int, int, int, int]]:
    """PD of the closure (edges 1..2n; edge k enters at pass k)."""
    steps = _walk_braid(word, m, p)
    E = len(steps)
    over, under = {}, {}
    for idx, (c, is_over, _) in enumerate(steps, start=1):
        (over if is_over else under)[c] = idx
    pd = []
    for c, g in enumerate(word):
        u, o = under[c], over[c]
        n_under, n_over = u % E + 1, o % E + 1
        pd.append((u, n_over, n_under, o) if g > 0 else (u, o, n_under, n_over))
    return pd


def random_braid(n: int, m: int, rng: random.Random, require_theta_nonzero: bool = True,
                     attempts: int = 1000) -> list[int]:
    """Word of n letters on m strands whose closure is ONE knot (m-cyclic permutation).

    Hardening (7d): by default it rejects knots with θ = 0 (trivial knot, amphichiral ones)
    and draws again, up to `attempts` times."""
    if m < 2 or n < m - 1 or (n - (m - 1)) % 2:
        raise ValueError("an m-cycle needs n >= m-1 letters and n ≡ m-1 (mod 2)")
    for _ in range(attempts):
        w = _cyclic_word(n, m, rng)
        if not require_theta_nonzero or theta_nonzero(braid_to_upright(w, m), rng):
            return w
    raise ValueError(f"{attempts} consecutive words with θ = 0 (n={n}, m={m}): too few crossings")


def _cyclic_word(n: int, m: int, rng: random.Random) -> list[int]:
    while True:
        w = [rng.choice((-1, 1)) * rng.randint(1, m - 1) for _ in range(n)]
        final_pos = {}                      # permutation bottom position -> top position
        for start in range(1, m + 1):
            pos = start
            for g in w:
                k = abs(g)
                pos = k + 1 if pos == k else (k if pos == k + 1 else pos)
            final_pos[start] = pos
        x, cycle_len = 1, 0
        while True:
            x = final_pos[x]
            cycle_len += 1
            if x == 1:
                break
        if cycle_len == m:
            return w


# ==========================================================================
# 7c. Pointwise evaluation of θ in F_q with numpy (q < 2^31 so that a·b fits in int64)
# ==========================================================================
import numpy as np

Q31 = (1 << 31) - 1


class _Mod:
    """Vector/scalar modulo q with the 4 operations; lets F1, F2 be written as in the paper."""
    __slots__ = ("v", "q")
    __array_ufunc__ = None            # ndarray (op) _Mod  ->  calls _Mod's __rop__

    def __init__(self, v, q: int) -> None:
        self.v, self.q = v, q

    def _u(self, o):
        return o.v if isinstance(o, _Mod) else o

    def __add__(self, o): return _Mod((self.v + self._u(o)) % self.q, self.q)
    def __sub__(self, o): return _Mod((self.v - self._u(o)) % self.q, self.q)
    def __rsub__(self, o): return _Mod((self._u(o) - self.v) % self.q, self.q)
    def __mul__(self, o): return _Mod((self.v * self._u(o)) % self.q, self.q)
    def __neg__(self): return _Mod((-self.v) % self.q, self.q)
    __radd__, __rmul__ = __add__, __mul__


def _inverse_np(A: np.ndarray, q: int) -> tuple[int, np.ndarray]:
    """(det, A^-1) mod q by Gauss-Jordan touching only the non-zero rows/columns (A is very sparse)."""
    m = A.shape[0]
    a = np.concatenate([A, np.eye(m, dtype=np.int64)], axis=1)
    det = 1
    for k in range(m):
        nz = np.flatnonzero(a[k:, k])
        if nz.size == 0:
            raise ZeroDivisionError("A(T) singular: T is a root of Δ")
        p = k + int(nz[0])
        if p != k:
            a[[k, p]] = a[[p, k]]
            det = -det
        pv = int(a[k, k])
        det = det * pv % q
        cols = np.flatnonzero(a[k])
        a[k, cols] = a[k, cols] * pow(pv, -1, q) % q
        rows = np.flatnonzero(a[:, k])
        rows = rows[rows != k]
        if rows.size:
            f = a[rows, k][:, None]
            a[np.ix_(rows, cols)] = (a[np.ix_(rows, cols)] - f * a[k, cols][None, :]) % q
    return det % q, a[:, m:]


def _mulmod(A: np.ndarray, B: np.ndarray, q: int) -> np.ndarray:
    """A·B mod q (q < 2^31) with EXACT float64 BLAS (16-bit limbs => every sum < 2^53). Accepts stacks (B,m,m)."""
    a1, a0 = (A >> 16).astype(np.float64), (A & 0xFFFF).astype(np.float64)
    b1, b0 = (B >> 16).astype(np.float64), (B & 0xFFFF).astype(np.float64)
    hh = (a1 @ b1).astype(np.int64) % q
    hl = (np.concatenate([a1, a0], axis=-1) @ np.concatenate([b0, b1], axis=-2)).astype(np.int64) % q
    ll = (a0 @ b0).astype(np.int64)
    return (hh * pow(2, 32, q) + (hl << 16) + ll) % q


def _inverse_dense(M: np.ndarray, q: int) -> tuple[np.ndarray, np.ndarray]:
    """Batched dense Gauss-Jordan (stack (B,m,m)): few numpy calls per step, shared by the batch."""
    nb, m, _ = M.shape
    a = np.concatenate([M, np.broadcast_to(np.eye(m, dtype=np.int64), (nb, m, m))], axis=2)
    ar = np.arange(nb)
    det = [1] * nb
    for k in range(m):
        nz = a[:, k:, k] != 0
        p = k + nz.argmax(axis=1)
        if not nz[ar, p - k].all():
            raise ZeroDivisionError("singular block")
        if (p != k).any():
            fk, fp = a[ar, k].copy(), a[ar, p].copy()
            a[ar, k], a[ar, p] = fp, fk
            for b in np.nonzero(p != k)[0]:
                det[b] = -det[b]
        pv = a[:, k, k].tolist()
        ip = np.array([pow(x, -1, q) for x in pv], dtype=np.int64)
        det = [d * x % q for d, x in zip(det, pv)]
        a[:, k] = a[:, k] * ip[:, None] % q
        f = a[:, :, k].copy()
        f[:, k] = 0
        a = (a - f[:, :, None] * a[:, k][:, None, :]) % q
    return np.array(det, dtype=np.int64) % q, a[:, :, m:]


_BASE = 32


def _inverse_blocks(M: np.ndarray, q: int) -> tuple[np.ndarray, np.ndarray]:
    """(det, M^-1) of a stack (B,m,m) by recursive Schur complement: O(m^3), almost all in BLAS matmul."""
    m = M.shape[-1]
    if m <= _BASE:
        return _inverse_dense(M, q)
    h = m // 2
    A, B, C, D = M[:, :h, :h], M[:, :h, h:], M[:, h:, :h], M[:, h:, h:]
    dA, Ai = _inverse_blocks(A, q)
    CAi = _mulmod(C, Ai, q)
    dS, Si = _inverse_blocks((D - _mulmod(CAi, B, q)) % q, q)
    AiB = _mulmod(Ai, B, q)
    AiBSi = _mulmod(AiB, Si, q)
    inv = np.empty_like(M)
    inv[:, :h, :h] = (Ai + _mulmod(AiBSi, CAi, q)) % q
    inv[:, :h, h:] = (-AiBSi) % q
    inv[:, h:, :h] = (-_mulmod(Si, CAi, q)) % q
    inv[:, h:, h:] = Si
    return dA * dS % q, inv


def _inverse_fast(M: np.ndarray, q: int) -> tuple[np.ndarray, np.ndarray]:
    """Inverts a stack of matrices; if some leading block is singular, falls back to full pivoting matrix by matrix."""
    try:
        return _inverse_blocks(M, q)
    except ZeroDivisionError:
        res = [_inverse_np(x, q) for x in M]
        return np.array([r[0] for r in res], dtype=np.int64), np.stack([r[1] for r in res])


def theta_eval_np(d: KnotDiagram, T1: int, T2: int, q: int = Q31, profile: dict | None = None) -> int:
    """θ(D)(T1,T2) ∈ F_q, same as theta_eval but with vectorised linear algebra."""
    n, m = d.n, 2 * d.n + 1
    t0 = time.perf_counter()
    inv_ = lambda x: pow(x % q, -1, q)
    T3 = T1 * T2 % q
    if T2 % q == 1:
        raise ZeroDivisionError("T2 = 1 is an apparent pole")
    s = np.array([c[0] for c in d.crossings], dtype=np.int64)
    i0 = np.array([c[1] for c in d.crossings], dtype=np.int64) - 1
    j0 = np.array([c[2] for c in d.crossings], dtype=np.int64) - 1
    mats, Ts = [], (T1, T2, T3)
    for T in Ts:
        Tc = np.where(s > 0, T, inv_(T)).astype(np.int64)
        A = np.eye(m, dtype=np.int64)
        A[i0, i0 + 1] = (A[i0, i0 + 1] - Tc) % q
        A[i0, j0 + 1] = (A[i0, j0 + 1] + Tc - 1) % q
        A[j0, j0 + 1] = (A[j0, j0 + 1] - 1) % q
        mats.append(A)
    dets, invs = _inverse_fast(np.stack(mats), q)
    res = [(int(dets[k]), invs[k]) for k in range(3)]
    G, Del = [], 1
    for T, (det, inv) in zip(Ts, res):
        G.append(inv)
        Del = Del * pow(T, d.delta_exponent, q) % q * det % q
    t1 = time.perf_counter()
    G1, G2, G3 = G
    g = lambda M, a, b: _Mod(M[a, b], q)
    pw = lambda T: _Mod(np.where(s > 0, T, inv_(T)).astype(np.int64), q)      # T^s per crossing
    t1s, t2s, t3s = pw(T1), pw(T2), pw(T3)
    S = _Mod(s % q, q)
    inv_t2m1 = _Mod(np.where(s > 0, inv_(T2 - 1), inv_(inv_(T2) - 1)).astype(np.int64), q)
    half = (q + 1) // 2
    # ---- F1 (sum over crossings) ----
    g1ii, g2ii, g3ii = g(G1, i0, i0), g(G2, i0, i0), g(G3, i0, i0)
    g1ji, g2ji, g3ji = g(G1, j0, i0), g(G2, j0, i0), g(G3, j0, i0)
    g2ij = g(G2, i0, j0)
    g2jj, g3jj = g(G2, j0, j0), g(G3, j0, j0)
    poly = (half - g3ii + t2s * g1ii * g2ji - t2s * g3jj * g2ji - (t2s - 1) * g3ii * g2ji
            + (t3s - 1) * g2ji * g3ji - g1ii * g2jj + 2 * g3ii * g2jj + g1ii * g3jj - g2ii * g3jj)
    frac = ((t1s - 1) * t2s * (g3jj * g1ji - g2jj * g1ji + t2s * g1ji * g2ji)
            + (t3s - 1) * g3ji * (1 - t2s * g1ii + g2ij + (t2s - 2) * g2jj - (t1s - 1) * (t2s + 1) * g1ji))
    th0 = int(np.sum(S.v * ((poly + frac * inv_t2m1).v) % q)) % q
    # ---- F2 (sum over pairs of crossings: c0 on axis 0, c1 on axis 1) ----
    I0, J0, I1, J1 = i0[:, None], j0[:, None], i0[None, :], j0[None, :]
    col = lambda x: _Mod(x.v[:, None], q)
    row = lambda x: _Mod(x.v[None, :], q)
    interior = (col(t2s) * _Mod(G2[I1, I0], q) + _Mod(G2[J1, J0], q)
                - col(t2s) * _Mod(G2[J1, I0], q) - _Mod(G2[I1, J0], q))
    pref = col(t1s - 1) * row(S * (t3s - 1) * inv_t2m1) * _Mod(G1[J1, I0], q) * _Mod(G3[J0, I1], q)
    th0 += int(np.sum((pref * interior).v)) % q
    # ---- F3 (sum over edges) ----
    phi = np.array(d.rotations, dtype=np.int64) % q
    th0 += int(np.sum(phi * ((np.diag(G3) - half) % q) % q)) % q
    t2 = time.perf_counter()
    if profile is not None:
        profile["inversions"], profile["F1+F2+F3"] = t1 - t0, t2 - t1
    return th0 % q * Del % q


# ==========================================================================
# 7d. Hardening: θ ≠ 0 filter and SHA-256 output
# ==========================================================================
HASH_DOMAIN = b"holocrypto-v1/theta"


def theta_nonzero(d: KnotDiagram, rng: random.Random, points: int = 2, q: int = Q31) -> bool:
    """False if θ vanishes at any of `points` random points of F_q.

    If θ ≢ 0 has total degree D, it vanishes at a random point with probability ≤ D/q
    (Schwartz-Zippel; D = O(n), q ≈ 2^31), so a knot is wrongly rejected only with negligible
    probability. If θ ≡ 0 (trivial knot, amphichiral ones) it is always rejected."""
    done = 0
    while done < points:
        T1, T2 = rng.randrange(2, q), rng.randrange(2, q)
        try:
            v = theta_eval_np(d, T1, T2, q)
        except ZeroDivisionError:             # T is a root of det A: another point
            continue
        if v == 0:
            return False
        done += 1
    return True


def holocrypto_hash(d: KnotDiagram, T1: int, T2: int, q: int = Q31) -> str:
    """SHA-256( domain || q || T1 || T2 || θ(T1,T2) ), each integer in 8 bytes big-endian.

    θ is hashed as is, NOT |θ| nor θ²: with θ² a knot and its mirror (θ and -θ) would give the SAME
    hash, which is precisely a collision. With SHA-256 the relation θ(mirror) = -θ is no longer
    visible in the output: the two hashes have no algebraic relation."""
    th = theta_eval_np(d, T1, T2, q)
    msg = HASH_DOMAIN + b"".join(x.to_bytes(8, "big") for x in (q, T1 % q, T2 % q, th))
    return hashlib.sha256(msg).hexdigest()


# ==========================================================================
# 8. Verification
# ==========================================================================
def verify_base() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    print("== A. Paper input: matrix A and its inverse (trefoil of the 'upright' figure) ==")
    d = paper_trefoil()
    T = Fraction(3)
    det, inv = det_inv_field(d.matrix_A(T, Fraction(1)))
    D2 = T * T - T + 1
    expected_G = [[1, T, 1, T, 1, T, 1],
                  [0, 1, 1 / D2, T / D2, T / D2, T * T / D2, 1],
                  [0, 0, 1 / D2, T / D2, T / D2, T * T / D2, 1],
                  [0, 0, (1 - T) / D2, 1 / D2, 1 / D2, T / D2, 1],
                  [0, 0, (1 - T) / D2, (T - T * T) / D2, 1 / D2, T / D2, 1],
                  [0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 1]]
    ok = all(inv[r][c] == Fraction(expected_G[r][c]) for r in range(7) for c in range(7))
    print(f"  G = A⁻¹ at T=3 matches the matrix printed in the paper: {ok};  det A = T²-T+1 = {D2}: {det == D2}")
    assert ok and det == D2

    print("\n== B. Symbolic Δ (exact Bareiss) vs known Alexander polynomials ==")
    print(f"  {'knot':<5} {'n':>2} {'Δ(1)':>5} {'symmetric':>9} {'= known':>11}   Δ")
    all_ok = True
    for name in KNOTS:
        dg = knot(name)
        t0 = time.perf_counter()
        delta = alexander(dg)
        dt = time.perf_counter() - t0
        sym = all(delta.c.get(-e, 0) == v for e, v in delta.c.items())
        ok = delta == ALEXANDER_KNOWN[name]
        all_ok &= ok
        print(f"  {name:<5} {dg.n:>2} {delta(Fraction(1))!s:>5} {sym!s:>9} {ok!s:>11}   {delta}  ({dt * 1e3:.1f} ms)")
    assert all_ok

    print("\n== C. θ of the trefoil vs the value published by van der Veen ==")
    d_ref = KnotDiagram.from_sage([[1, 0, 3], [1, 4, 1], [1, 2, 5]], [0, 0, 0, -1, 0, 0, 0], "3_1 ref")
    poly_ref = theta_polynomial(d_ref)
    print(f"  θ computed (reference trefoil):\n    {format_poly(poly_ref)}")
    print(f"  matches the published polynomial EXACTLY: {poly_ref == THETA_TREFOIL_REFERENCE}")
    assert poly_ref == THETA_TREFOIL_REFERENCE
    print(f"  paper trefoil (same figure) identical: {theta_polynomial(paper_trefoil()) == poly_ref}")

    print("\n== D. Paper's conjecture: θ(mirror) = -θ(K)  (=> θ = 0 for amphichiral knots) ==")
    th31 = theta_polynomial(knot("3_1"))            # 3_1 of knots.sage = left-handed trefoil
    negated = {k: -v for k, v in poly_ref.items()}
    print(f"  θ(left trefoil) == -θ(right trefoil): {th31 == negated}")
    assert th31 == negated
    for name in ("4_1", "6_3", "8_19", "8_20", "5_1", "5_2"):
        t0 = time.perf_counter()
        p = theta_polynomial(knot(name))
        amph = name in ("4_1", "6_3")
        print(f"  θ({name}): {len(p):>3} monomials  ({time.perf_counter() - t0:.1f}s)"
              + ("   <- amphichiral: θ = 0 expected" if amph else ""))
        assert (len(p) == 0) == amph

    print("\n== E. The apparent poles at T2=1 cancel (θ is a Laurent polynomial) ==")
    for name in ("5_1", "5_2"):
        dg = knot(name)
        vals = [theta_eval(dg, Fraction(7, 3), 1 + Fraction(1, 10 ** k)) for k in (3, 6, 9)]
        print(f"  {name}: θ(7/3, 1+10⁻ᵏ) k=3,6,9 -> {[float(v) for v in vals]}")

    print("\n== F. Modular arithmetic: same θ in F_q (q = 2^61-1) ==")
    q = (1 << 61) - 1
    dg = knot("5_2")
    pm = theta_eval(dg, Fp(123456789, q), Fp(987654321, q))
    p52 = theta_polynomial(dg)
    pe = sum((Fp(v, q) * Fp(123456789, q) ** a * Fp(987654321, q) ** b for (a, b), v in p52.items()), Fp(0, q))
    print(f"  θ(5_2) in F_q: direct evaluation == exact polynomial reduced mod q: {pm == pe}")
    assert pm == pe


# ==========================================================================
# 9. Verification of the input module and of the pointwise evaluation
# ==========================================================================
DT_KNOWN = {                       # KnotInfo (negative sign: the even pass goes over)
    "3_1": [4, 6, 2], "4_1": [4, 6, 8, 2], "5_1": [6, 8, 10, 2, 4], "5_2": [4, 8, 10, 2, 6],
    "6_1": [4, 8, 12, 10, 2, 6], "6_2": [4, 8, 10, 12, 2, 6], "6_3": [4, 8, 10, 2, 12, 6],
    "7_1": [8, 10, 12, 14, 2, 4, 6], "8_19": [4, 8, -12, 2, -14, -16, -6, -10],
}
P_RATIONAL = (Fraction(7, 3), Fraction(5, 4))


def _theta_rational(d: KnotDiagram) -> Fraction:
    return theta_eval(d, *P_RATIONAL)


def verify_input() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    rng = random.Random(2026)

    print("== G. Braid -> direct upright: reproduces the paper's trefoil ==")
    d = braid_to_upright([1, 1, 1], 2)
    ok = d.crossings == paper_trefoil().crossings and d.rotations == paper_trefoil().rotations
    print(f"  σ1³ on 2 strands == X={{(1,1,4),(1,5,2),(1,3,6)}}, φ=(0,0,0,-1,0,0,0): {ok}")
    assert ok

    print("\n== H. PD -> upright (flow on the dual) == direct braid, for different cuts, trees and Sw twists ==")
    cases = [([1, 1, 1], 2), ([1, -2, 1, -2], 3), ([1, 1, 1, 1, 1], 2)]
    while len(cases) < 12:
        m = rng.choice((3, 4, 5))
        n_ = rng.randint(8, 13)
        cases.append((random_braid(n_ + (n_ - (m - 1)) % 2, m, rng), m))
    total = failures = 0
    for w, m in cases:
        base = _theta_rational(braid_to_upright(w, m))
        d_base = alexander(braid_to_upright(w, m))
        pd = braid_to_pd(w, m)
        assert PlanarMap(pd).is_planar
        for _ in range(6):
            e0 = rng.randint(1, 2 * len(w))
            du = pd_to_upright(pd, e0, random.Random(rng.random()), twists=rng.randint(0, 4))
            total += 1
            if _theta_rational(du) != base or alexander(du) != d_base:
                failures += 1
    print(f"  {len(cases)} knots x 6 variants (cut edge, dual tree and random 2π twists): "
          f"{total - failures}/{total} with Δ and θ identical to the direct braid")
    assert failures == 0

    print("\n== I. Reidemeister invariance (validates the rotations in kinks, R2 and R3) ==")
    ref = _theta_rational(braid_to_upright([1, 1, 1], 2))
    tests = {
        "Markov +σ2 (kink)": ([1, 1, 1, 2], 3), "Markov -σ2 (kink)": ([1, 1, 1, -2], 3),
        "R2 inserted": ([1, 1, 1, 2, 2, -2], 3), "R2 inserted (σ1σ1⁻¹)": ([1, 1, -1, 1, 1, 2], 3),
        "conjugated by σ2": ([2, 1, 1, 1, 2, -2], 3),
    }
    for name, (w, m) in tests.items():
        try:
            v = _theta_rational(braid_to_upright(w, m))
        except ValueError:
            print(f"  {name:<26} (word is not a knot, skipped)")
            continue
        print(f"  {name:<26} θ == θ(trefoil): {v == ref}")
        assert v == ref
    assert _theta_rational(braid_to_upright([-2, 1, -2, 1], 3)) == _theta_rational(braid_to_upright([1, -2, 1, -2], 3))
    made = 0
    while made < 6:                                   # R3: σ1σ2σ1 <-> σ2σ1σ2 inside a random word
        w = [rng.choice((-1, 1)) * rng.randint(1, 2) for _ in range(rng.choice((5, 7)))]
        t = rng.randrange(len(w) + 1)
        wa, wb = w[:t] + [1, 2, 1] + w[t:], w[:t] + [2, 1, 2] + w[t:]
        try:
            va, vb = _theta_rational(braid_to_upright(wa, 3)), _theta_rational(braid_to_upright(wb, 3))
        except ValueError:
            continue
        assert va == vb, (wa, wb)
        made += 1
    print("  R3 (σ1σ2σ1 <-> σ2σ1σ2) on 6 random knots: identical θ; conjugation of 4_1: identical")

    print("\n== J. DT -> PD -> upright (chirality is determined by planarity) ==")
    print(f"  {'knot':<5} {'DT':<34} {'Δ = Δ_known':>14}   θ relative to van der Veen's table")
    for name, dt in DT_KNOWN.items():
        t0 = time.perf_counter()
        pd = dt_to_pd(dt)
        dg = pd_to_upright(pd, 1, random.Random(1), name=name)
        dt_s = time.perf_counter() - t0
        ok_d = alexander(dg) == ALEXANDER_KNOWN[name]
        v, tab = _theta_rational(dg), _theta_rational(knot(name))
        rel = "= θ" if v == tab and v != 0 else ("= -θ (mirror image)" if v == -tab and v != 0
                                                  else ("= 0 (amphichiral)" if v == 0 == tab else "DIFFERENT"))
        print(f"  {name:<5} {str(dt):<34} {ok_d!s:>14}   {rel}   [{dt_s * 1e3:.0f} ms]")
        assert ok_d and not rel.startswith("DIFF")

    print("\n== K. numpy backend (F_q, q=2^31-1) == exact Fp backend ==")
    diagrams = [knot("3_1"), knot("5_2"), knot("8_20")]
    for _ in range(3):
        m = rng.choice((3, 5))
        diagrams.append(braid_to_upright(random_braid(2 * rng.randint(6, 9), m, rng), m))
    for dg in diagrams:
        T1, T2 = rng.randrange(2, Q31), rng.randrange(2, Q31)
        a = theta_eval_np(dg, T1, T2)
        b = theta_eval(dg, Fp(T1, Q31), Fp(T2, Q31)).v
        assert a == b, (dg.n, a, b)
    print(f"  {len(diagrams)} diagrams (n = {[dg.n for dg in diagrams]}): they agree")

    print("\n== L. Large knots: direct braid == PD with Sw twists (numpy evaluation in F_q) ==")
    for n, m in ((30, 5), (60, 7), (100, 9)):
        w = random_braid(n, m, rng)
        dg = braid_to_upright(w, m)
        T1, T2 = rng.randrange(2, Q31), rng.randrange(2, Q31)
        a = theta_eval_np(dg, T1, T2)
        du = pd_to_upright(braid_to_pd(w, m), rng.randint(1, 2 * n), random.Random(5), twists=8)
        b = theta_eval_np(du, T1, T2)
        print(f"  n={n:<4} m={m:<2} θ(braid) == θ(PD, random cut and twists): {a == b}")
        assert a == b


# ==========================================================================
# 10. Benchmark of the pointwise evaluation
# ==========================================================================
def benchmark(repetitions: int = 7) -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    rng = random.Random(7)
    print(f"== Evaluation of θ at a random point of F_q (q = 2^31-1), median of {repetitions} points ==")
    print("   (only theta_eval_np is timed; building the diagram does NOT count)\n")
    print(f"  {'n':>4} {'strands':>8} {'matrix':>7} {'median(ms)':>12} {'min':>8} {'max':>8}"
          f" {'inversions':>12} {'F1+F2+F3':>9}")
    median_100 = None
    for n in (20, 50, 100, 200):
        for m in (5, 11):
            w = random_braid(n, m, rng)
            d = braid_to_upright(w, m)
            theta_eval_np(d, 12345, 67890)                       # warm-up
            ts, profiles = [], []
            for _ in range(repetitions):
                T1, T2 = rng.randrange(2, Q31), rng.randrange(2, Q31)
                prof: dict = {}
                t0 = time.perf_counter()
                theta_eval_np(d, T1, T2, profile=prof)
                ts.append((time.perf_counter() - t0) * 1e3)
                profiles.append(prof)
            med = sorted(ts)[len(ts) // 2]
            pm = profiles[ts.index(med)]
            print(f"  {n:>4} {m:>8} {2 * n + 1:>6}² {med:>12.1f} {min(ts):>8.1f} {max(ts):>8.1f}"
                  f" {pm['inversions'] * 1e3:>10.1f}ms {pm['F1+F2+F3'] * 1e3:>7.1f}ms")
            if n == 100:
                median_100 = med if median_100 is None else max(median_100, med)
    print("\n  Reference, exact Python backend with q = 2^61-1 (Fp, no numpy):")
    for n, m in ((20, 5), (50, 5)):
        d = braid_to_upright(random_braid(n, m, rng), m)
        T1, T2 = Fp(rng.randrange(2, (1 << 61) - 1), (1 << 61) - 1), Fp(rng.randrange(2, (1 << 61) - 1), (1 << 61) - 1)
        t0 = time.perf_counter()
        theta_eval(d, T1, T2)
        print(f"    n={n:<4} {((time.perf_counter() - t0) * 1e3):>10.0f} ms")
    verdict = "MET" if median_100 is not None and median_100 < 50 else "NOT met"
    print(f"\n  Target < 50 ms at n=100 (worse of the two widths, median): {median_100:.1f} ms -> {verdict}")


# ==========================================================================
# 11. Verification of the hardening (7d)
# ==========================================================================
def verify_hardening() -> None:
    sys.stdout.reconfigure(encoding="utf-8")
    rng = random.Random(11)

    print("== M. θ ≠ 0 filter in key generation ==")
    cases = {"trivial knot (σ1, kink)": (braid_to_upright([1], 2), False),
             "4_1 (amphichiral)": (knot("4_1"), False), "6_3 (amphichiral)": (knot("6_3"), False),
             "3_1": (knot("3_1"), True), "5_2": (knot("5_2"), True), "8_19": (knot("8_19"), True)}
    for name, (dg, expected) in cases.items():
        v = theta_nonzero(dg, rng)
        print(f"  {name:<26} accepted: {v!s:<5} (expected {expected})")
        assert v == expected
    for n, m in ((4, 3), (5, 2), (6, 3)):
        raw = [random_braid(n, m, rng, require_theta_nonzero=False) for _ in range(300)]
        nulls = sum(not theta_nonzero(braid_to_upright(w, m), rng) for w in raw)
        filtered = [random_braid(n, m, rng) for _ in range(100)]
        zero = sum(not theta_nonzero(braid_to_upright(w, m), random.Random(rng.random()), points=4)
                   for w in filtered)
        print(f"  n={n} m={m}: without filter {nulls:>3}/300 with θ = 0;  with filter {zero}/100 with θ = 0")
        assert zero == 0

    print("\n== N. SHA-256 output: the mirror no longer shows a visible -θ, and no collision is created ==")
    pairs = equal = 0
    for _ in range(20):
        n, m = rng.choice(((20, 5), (30, 7), (16, 3)))
        w = random_braid(n, m, rng)
        dk, d_mir = braid_to_upright(w, m), braid_to_upright([-g for g in w], m)
        T1, T2 = rng.randrange(2, Q31), rng.randrange(2, Q31)
        a, b = theta_eval_np(dk, T1, T2), theta_eval_np(d_mir, T1, T2)
        assert b == (-a) % Q31 and a != 0                    # raw relation θ(mirror) = -θ
        assert a * a % Q31 == b * b % Q31                    # with θ² they would be equal: collision
        hk, he = holocrypto_hash(dk, T1, T2), holocrypto_hash(d_mir, T1, T2)
        assert len(hk) == 64 and hk != he
        du = pd_to_upright(braid_to_pd(w, m), rng.randint(1, 2 * n), random.Random(rng.random()), twists=4)
        equal += holocrypto_hash(du, T1, T2) == hk         # same knot, other diagram
        pairs += 1
    print(f"  {pairs} random knots (n = 16..30) and their mirrors:")
    print("    raw θ:   θ(mirror) = -θ in all             -> algebraic relation visible")
    print("    θ²:      θ(mirror)² = θ² in all            -> would be a collision (option discarded)")
    print("    SHA-256: hash(K) ≠ hash(mirror) in all     -> no visible relation, no collision")
    print(f"    same knot via another diagram (PD, cut and twists): same hash in {equal}/{pairs}")
    assert equal == pairs
    T1, T2 = 123456789, 987654321
    print(f"  example: hash(right trefoil)        = {holocrypto_hash(braid_to_upright([1, 1, 1], 2), T1, T2)}")
    print(f"           hash(left trefoil)         = {holocrypto_hash(braid_to_upright([-1, -1, -1], 2), T1, T2)}")
    print(f"           hash(3_1 from table, left) = {holocrypto_hash(knot('3_1'), T1, T2)}")


def main() -> None:
    args = set(sys.argv[1:])
    run_all = not args
    if run_all or "--base" in args:
        verify_base()
    if run_all or "--input" in args:
        verify_input()
    if run_all or "--hardening" in args:
        verify_hardening()
    if run_all or "--bench" in args:
        benchmark()


if __name__ == "__main__":
    main()
