"""non-BPS triality criterion: x~ in the (Gamma x| S3)-orbit of x iff gcd(x+x~) = gcd(x-x~) and there is an
involution pi in S3 (identity or a transposition) with -Q_j ~ Q_pi(j) for all j (proper equivalence)."""
import itertools, sys
from math import gcd
from functools import reduce
from canonsq import canon_sq
from cube import T, form
from split import triples as striples
from splitcube import cube_from as scube
from triality import perm, P
from criteria import fclass
def g(v): return reduce(gcd, [abs(t) for t in v])
M = int(sys.argv[1]) if len(sys.argv) > 1 else 12
bad = tot = 0
for m in range(1, M+1):
    for x in striples(m):
        if not x['fint']: continue
        a = scube(m, x['J1'], x['J2'], x['I3']); ad = tuple(-t // m for t in T(a))
        c = canon_sq(ad)
        same3 = any(canon_sq(perm(a, s)) == c for s in P)
        Q = [fclass(form(a, j), m) for j in range(3)]
        mQ = [fclass(tuple(-u for u in form(a, j)), m) for j in range(3)]
        invols = [(0, 1, 2), (1, 0, 2), (2, 1, 0), (0, 2, 1)]
        crit = g([u + v for u, v in zip(a, ad)]) == g([u - v for u, v in zip(a, ad)]) and \
               any(all(mQ[j] == Q[p[j]] for j in range(3)) for p in invols)
        tot += 1; bad += (crit != same3)
print(f"non-BPS m <= {M}: {tot} F-integral classes; triality criterion failures: {bad}")
