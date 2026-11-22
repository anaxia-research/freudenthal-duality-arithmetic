"""Completeness of the enumeration: every ST^2 charge with entries in [-R, R] and Delta = 4k^2 or -m^2 (k, m <= K)
lies in one of the enumerated orbits, with the same integrality and U-duality flags."""
import sys, itertools, numpy as np
from st2 import *
R = int(sys.argv[1]); K = int(sys.argv[2])
rng = np.arange(-R, R + 1, dtype=np.int64)
a, b, c, e, f, g = np.meshgrid(rng, rng, rng, rng, rng, rng, indexing='ij')
a, b, c, e, f, g = (t.ravel() for t in (a, b, c, e, f, g))
D = (-a*a*g*g + 4*a*b*f*g + 2*a*c*e*g - 4*a*c*f*f - 4*b*b*e*g + 4*b*c*e*f - c*c*e*e)
ok = np.zeros(D.shape, dtype=bool)
pos = D > 0; r = np.round(np.sqrt(np.abs(D)/4.0)).astype(np.int64)
ok |= pos & (4*r*r == D) & (r <= K)
r2 = np.round(np.sqrt(np.abs(D))).astype(np.int64)
ok |= (D < 0) & (r2*r2 == -D) & (r2 <= K)
idx = np.nonzero(ok)[0]
print("charges in the box [-%d, %d]^6 with Delta = 4k^2 or -m^2, k, m <= %d:" % (R, R, K), len(idx))
tables = {}
def flags(Dq):
    if Dq not in tables:
        O = orbits(Dq)
        tables[Dq] = {s: (dual(s) is not None, dual(s) is not None and canon(dual(s)) == s) for s in O}
    return tables[Dq]
nF = nU = 0
for i in idx:
    s = (int(a[i]), int(b[i]), int(c[i]), int(e[i]), int(f[i]), int(g[i]))
    Dq = Delta(s); assert Dq == int(D[i])
    t = flags(Dq); cs = canon(s)
    assert cs in t, s
    d = dual(s); fl = (d is not None, d is not None and canon(d) == cs)
    assert fl == t[cs], (s, fl, t[cs])
    nF += fl[0]; nU += fl[1]
print("all lie in enumerated orbits with matching flags; %d have an integral dual, %d are U-dual to it" % (nF, nU))
print("distinct values of Delta met:", len(tables))
