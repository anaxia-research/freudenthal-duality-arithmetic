"""Lattices in Q^2 = Q(i) (or Q(sqrt D)) with rational bases; Hermite normal form, sums, intersections."""
from fractions import Fraction as Fr
from math import gcd

def lcm(a, b): return a // gcd(a, b) * b

def hnf_int(rows):
    """Row HNF of an integer matrix with 2 columns and full rank 2: returns ((a,b),(0,c)), a,c>0, 0<=b<c."""
    rows = [list(r) for r in rows if r[0] != 0 or r[1] != 0]
    # column 0 gcd elimination
    while sum(1 for r in rows if r[0] != 0) > 1:
        rows.sort(key=lambda r: (r[0] == 0, abs(r[0])))
        p = rows[0]
        for r in rows[1:]:
            if r[0] != 0:
                q = r[0] // p[0]
                r[0] -= q*p[0]; r[1] -= q*p[1]
    rows.sort(key=lambda r: (r[0] == 0, abs(r[0])))
    p = rows[0]
    if p[0] < 0: p = [-p[0], -p[1]]
    c = 0
    for r in rows[1:]:
        c = gcd(c, abs(r[1]))
    assert p[0] != 0 and c != 0, "not full rank"
    return ((p[0], p[1] % c), (0, c))

class Lat:
    __slots__ = ("a", "b", "c")
    def __init__(self, gens):
        d = 1
        for g in gens:
            for x in g:
                d = lcm(d, Fr(x).denominator)
        rows = [(int(Fr(g[0])*d), int(Fr(g[1])*d)) for g in gens]
        (a, b), (_, c) = hnf_int(rows)
        self.a, self.b, self.c = Fr(a, d), Fr(b, d), Fr(c, d)
    def basis(self): return ((self.a, self.b), (Fr(0), self.c))
    def key(self): return (self.a, self.b, self.c)
    def __eq__(self, o): return self.key() == o.key()
    def __hash__(self): return hash(self.key())
    def covol(self): return self.a * self.c
    def contains(self, v):
        x, y = Fr(v[0]), Fr(v[1])
        t = x / self.a
        if t.denominator != 1: return False
        r = y - t*self.b
        return (r / self.c).denominator == 1
    def contains_lat(self, o): return all(self.contains(v) for v in o.basis())
    def __add__(self, o): return Lat(list(self.basis()) + list(o.basis()))
    def dual(self):
        (a, b), (_, c) = self.basis()
        # B = [[a,b],[0,c]], B^{-T} rows
        det = a*c
        return Lat([(c/det, Fr(0)), (-b/det, a/det)])
    def __and__(self, o): return (self.dual() + o.dual()).dual()
    def __repr__(self): return f"Lat({self.a},{self.b};{self.c})"

def sublattices(L, d):
    """All sublattices of L of index d."""
    (a, b), (_, c) = L.basis()
    v1, v2 = (a, b), (Fr(0), c)
    out = []
    for p in range(1, d+1):
        if d % p: continue
        q = d // p
        for r in range(q):
            # rows p*v1 + r*v2, q*v2
            out.append(Lat([(p*v1[0] + r*v2[0], p*v1[1] + r*v2[1]), (q*v2[0], q*v2[1])]))
    return out
