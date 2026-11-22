"""Explicit cube criteria for U-equivalence of x and its F-dual.
non-BPS (Delta = -m^2): x~ ~ x  iff  every face form Q_j is properly equivalent to -Q_j  and
                        gcd(m x - T(x)) = gcd(m x + T(x)).
BPS (Delta = 4k^2):     x~ ~ x  iff  some face form has content k, and then x~ = Phi_j x with Phi_j in SL2(Z)."""
from math import gcd
from functools import reduce
from split import triples as striples
from splitcube import cube_from as scube
from bps import triples as btriples
from tricube import cube_from as bcube
from canonsq import canon_sq
from canon import canon
from cube import T, form, act, Det
import sys

def g(v): return reduce(gcd, [abs(t) for t in v])

def form_class(Q, m):
    from nonbps_check import form_class as fc
    return fc(Q, m)

M = int(sys.argv[1]) if len(sys.argv) > 1 else 10
K = int(sys.argv[2]) if len(sys.argv) > 2 else 12
bad = 0; tot = 0
import importlib
# local copy of the form-class routine (nonbps_check runs a loop on import)
from canonsq import isotropic, egcd
def fclass(Q, m):
    A, B, C = Q; best = None
    for v in isotropic(Q, m):
        v1, v2 = v
        g0, p, q = egcd(v2, v1)
        if g0 < 0: p, q = -p, -q
        w = (p, -q)
        Ev = lambda x, y: A*x*x + B*x*y + C*y*y
        A2 = Ev(*w); C2 = Ev(*v); B2 = Ev(w[0]+v[0], w[1]+v[1]) - A2 - C2
        key = (A2 % abs(B2), B2, 0)
        best = key if best is None or key < best else best
    return best
for m in range(1, M+1):
    for x in striples(m):
        if not x['fint']: continue
        a = scube(m, x['J1'], x['J2'], x['I3'])
        t = T(a)
        ad = tuple(-s // m for s in t)
        same = canon_sq(ad) == canon_sq(a)
        ci = all(fclass(form(a, j), m) == fclass(tuple(-c for c in form(a, j)), m) for j in range(3))
        cii = g([m*u - s for u, s in zip(a, t)]) == g([m*u + s for u, s in zip(a, t)])
        tot += 1
        bad += (same != (ci and cii))
print(f"non-BPS m <= {M}: {tot} F-integral classes, explicit criterion failures: {bad}")
bad = 0; tot = 0
for k in range(1, K+1):
    for x in btriples(k):
        if not x['fint']: continue
        a = bcube(k, x['J1'], x['J2'], x['I3'])
        ad = tuple(s // (2*k) for s in T(a))
        same = canon(ad) == canon(a)
        cont = [g(form(a, j)) for j in range(3)]
        crit = k in cont
        if crit:
            j = cont.index(k); p, q, r = form(a, j)
            Phi = ((-q//(2*k), -2*p//(2*k)), (2*r//(2*k), q//(2*k)))
            assert Phi[0][0]*Phi[1][1] - Phi[0][1]*Phi[1][0] == 1
            assert act(a, j, Phi) == ad
        tot += 1; bad += (same != crit)
print(f"BPS k <= {K}: {tot} F-integral positive classes, content criterion failures: {bad}; Phi_j integral gives x~ exactly")
