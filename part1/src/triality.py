"""Effect of the triality S3 (permutations of the three indices) on U-equivalence of F-dual pairs."""
import itertools
from canon import canon
from canonsq import canon_sq
from cube import T, idx
from bps import triples as btriples
from tricube import cube_from as bcube
from split import triples as striples
from splitcube import cube_from as scube

def perm(a, s):
    out = [0]*8
    for A, B, C in itertools.product(range(2), repeat=3):
        ind = (A, B, C)
        new = [0, 0, 0]
        for i in range(3): new[s[i]] = ind[i]
        out[idx(*new)] = a[idx(A, B, C)]
    return tuple(out)
P = list(itertools.permutations(range(3)))
for name, a, m, cf in (("Delta=256", (0, 1, 4, 0, 4, 0, 0, -4), 16, canon), ("Delta=576", (-4, -9, -2, 3, -2, 3, 0, 3), 24, canon)):
    ad = tuple(t // m for t in T(a))
    print(name, "x~ in Gamma x S3 orbit of x:", any(cf(perm(a, s)) == cf(ad) for s in P))
K = 16
for k in range(1, K+1):
    nF = nE = nE3 = 0
    for x in btriples(k):
        if not x['fint']: continue
        a = bcube(k, x['J1'], x['J2'], x['I3']); ad = tuple(t // (2*k) for t in T(a))
        c = canon(ad); nF += 1
        nE += canon(a) == c
        nE3 += any(canon(perm(a, s)) == c for s in P)
    if nE3 != nE: print(f"BPS k={k}: positive F-integral classes {nF}, U-dual {nE}, U-dual with triality {nE3}")
for m in range(1, 13):
    nF = nE = nE3 = 0
    for x in striples(m):
        if not x['fint']: continue
        a = scube(m, x['J1'], x['J2'], x['I3']); ad = tuple(-t // m for t in T(a))
        c = canon_sq(ad); nF += 1
        nE += canon_sq(a) == c
        nE3 += any(canon_sq(perm(a, s)) == c for s in P)
    if nE3 != nE: print(f"non-BPS m={m}: F-integral classes {nF}, U-dual {nE}, U-dual with triality {nE3}")
print("done")
