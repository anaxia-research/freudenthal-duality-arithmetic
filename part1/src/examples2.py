"""More explicit pairs: a Delta = 576 pair with equal gcds and identical face forms that is not U-dual;
U-dual non-BPS pairs with an explicit U-duality; the unit charges Delta = 4, -1, -4."""
from math import gcd
from functools import reduce
import itertools
from bps import triples as btriples
from tricube import cube_from as bcube
from split import triples as striples
from splitcube import cube_from as scube
from canon import canon
from canonsq import canon_sq
from cube import T, form, act, Det
from examples import small, g

print("BPS k = 12 (Delta = 576), not U-dual:")
seen = set()
for x in btriples(12):
    if x['fint'] and not x['equiv']:
        a = small(bcube(12, x['J1'], x['J2'], x['I3']), canon)
        ad = tuple(t // 24 for t in T(a))
        key = (canon(a), canon(ad))
        print("  x =", a, "gcd", g(a), "gcd T", g(T(a)), "forms", [form(a, j) for j in range(3)])
        print("  x~=", ad, "gcd", g(ad), "gcd T", g(T(ad)), "forms", [form(ad, j) for j in range(3)])
        break
SM = [((p, q), (r, s)) for p, q, r, s in itertools.product(range(-4, 5), repeat=4) if p*s - q*r == 1]
def find_g(a, b):
    """search g = (g1, g2, g3) in SL2(Z)^3 with small entries and g a = b, using the face forms"""
    cands = []
    for j in range(3):
        Q = form(a, j); Qb = form(b, j)
        cj = [h for h in SM if form(act(a, j, h), j) == Qb]
        cands.append(cj)
    for h1 in cands[0]:
        a1 = act(a, 0, h1)
        for h2 in cands[1]:
            a2 = act(a1, 1, h2)
            for h3 in cands[2]:
                if act(a2, 2, h3) == b: return (h1, h2, h3)
    return None
for m in (2, 5):
    print(f"non-BPS m = {m}: U-dual pairs with explicit g")
    n = 0
    for x in striples(m):
        if not x['fint']: continue
        a = small(scube(m, x['J1'], x['J2'], x['I3']), canon_sq)
        ad = tuple(-t // m for t in T(a))
        if canon_sq(a) != canon_sq(ad): continue
        gg = find_g(a, ad)
        print("  x =", a, " x~ =", ad, " forms", [form(a, j) for j in range(3)], " g =", gg)
        n += 1
        if n >= 3: break
for a in [(0, 1, 1, 0, 1, 0, 0, 0)]:
    print("Delta =", -Det(a), "x =", a, "T =", T(a), "forms", [form(a, j) for j in range(3)])
