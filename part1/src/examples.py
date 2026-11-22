"""Explicit examples of integral Freudenthal pairs."""
from math import gcd
from functools import reduce
from bps import triples as btriples
from tricube import cube_from as bcube
from split import triples as striples
from splitcube import cube_from as scube
from cube import Det, T, form
from canon import canon
from canonsq import canon_sq
import itertools

def content(Q): return reduce(gcd, [abs(x) for x in Q])
def g(v): return reduce(gcd, [abs(x) for x in v])
def small(a, canonf):
    """search a smaller representative by brute force over a few SL2 words"""
    from cube import act
    gens = [((0, -1), (1, 0)), ((1, 1), (0, 1)), ((1, -1), (0, 1))]
    best = a; seen = {a}; frontier = [a]
    for _ in range(4):
        nf = []
        for b in frontier:
            for j in range(3):
                for h in gens:
                    c = act(b, j, h)
                    if c not in seen:
                        seen.add(c); nf.append(c)
                        if sum(x*x for x in c) < sum(x*x for x in best): best = c
        frontier = nf
    assert canonf(best) == canonf(a)
    return best

print("BPS, Delta = 256 (k = 8): F-dual integral and not U-dual (contents (4,4,4))")
for x in btriples(8):
    if x['fint'] and not x['equiv']:
        a = small(bcube(8, x['J1'], x['J2'], x['I3']), canon)
        ad = tuple(t // 16 for t in T(a))
        print("  x  =", a, " Delta =", -Det(a), " gcd =", g(a), " forms", [form(a, j) for j in range(3)])
        print("  x~ =", ad, " Delta =", -Det(ad), " gcd =", g(ad), " forms", [form(ad, j) for j in range(3)],
              " same orbit:", canon(a) == canon(ad))
for m in (2, 3):
    print(f"non-BPS, Delta = -{m*m}: F-dual integral and not U-dual")
    for x in striples(m):
        if not x['fint']: continue
        a = small(scube(m, x['J1'], x['J2'], x['I3']), canon_sq)
        ad = tuple(-t // m for t in T(a))
        if canon_sq(a) == canon_sq(ad): continue
        print("  x  =", a, " Delta =", -Det(a), " gcd =", g(a), " forms", [form(a, j) for j in range(3)])
        print("  x~ =", ad, " gcd =", g(ad), " forms", [form(ad, j) for j in range(3)])
        break
