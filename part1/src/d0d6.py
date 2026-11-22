"""D0-D6 type charges x = e_000 + m e_111 (Delta = -m^2) and their F-duals; BPS analogue and gcd data."""
from cube import T, Det
from canonsq import canon_sq
from canon import canon
from math import gcd
from functools import reduce
from bps import triples
from tricube import cube_from
res = []
for m in range(1, 41):
    x = (1, 0, 0, 0, 0, 0, 0, m)
    xd = tuple(-t // m for t in T(x))
    assert xd == (1, 0, 0, 0, 0, 0, 0, -m)
    res.append((m, canon_sq(x) == canon_sq(xd)))
print("D0-D6 x=(1,0,...,0,m): F-dual (1,0,...,0,-m) U-equivalent for m in", [m for m, s in res if s],
      "; inequivalent for all other m <= 40")
g = lambda v: reduce(gcd, [abs(t) for t in v])
for k in range(8, 25):
    pairs = set()
    for x in triples(k):
        if x['fint'] and not x['equiv']:
            a = cube_from(k, x['J1'], x['J2'], x['I3'])
            ad = tuple(t // (2*k) for t in T(a))
            pairs.add((g(a), g(ad)))
    if pairs: print(f"k={k}: gcd(x), gcd(x~) over U-inequivalent F-dual pairs:", sorted(pairs))
