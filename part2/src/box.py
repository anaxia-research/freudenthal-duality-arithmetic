"""Completeness of the Hessian-reduction enumeration: every form with entries in [-B, B] and
Delta = 4k^2 (k <= K) or Delta = -m^2 (m <= K) reduces to an enumerated class, with the same flags."""
import sys, numpy as np
from cubic import *
from classes import classes_bps, classes_non, canon
B, K = int(sys.argv[1]), int(sys.argv[2])
ref = {}
for k in range(1, K + 1): ref[4*k*k] = classes_bps(k)
for m in range(1, K + 1): ref[-m*m] = classes_non(m)
r = np.arange(-B, B + 1, dtype=np.int64)
a, b, c, d = np.meshgrid(r, r, r, r, indexing='ij'); a, b, c, d = (x.ravel() for x in (a, b, c, d))
Dl = -(a*a*d*d - 6*a*b*c*d + 4*a*c**3 + 4*b**3*d - 3*b*b*c*c)
keys = np.array(sorted(ref)); hit = np.isin(Dl, keys)
met = {D: set() for D in ref}; n = 0
for i in np.nonzero(hit)[0]:
    f = (int(a[i]), int(b[i]), int(c[i]), int(d[i])); D = Delta(f)
    cf = canon(f); n += 1
    assert cf in ref[D], (f, D, cf)
    ft = dual(f)
    assert (ft is None) == (dual(cf) is None)
    if ft is not None: assert (canon(ft) == cf) == (canon(dual(cf)) == cf)
    met[D].add(cf)
for k in range(1, K + 1):
    print(f"Delta = {4*k*k:5d}: classes met {len(met[4*k*k]):3d} of {len(ref[4*k*k]):3d};   Delta = {-k*k:5d}: classes met {len(met[-k*k]):3d} of {len(ref[-k*k]):3d}")
print(n, "forms with entries in [-%d, %d] reduced; no form outside the enumerated classes; box check passed" % (B, B))
