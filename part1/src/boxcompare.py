"""Completeness check that does not use Bhargava's theorem (BPS part): every SL2(Z)^3-orbit met by a cube
with entries in [-B, B] and Det = -4k^2 is one of the orbits produced from balanced triples, with the same
integrality and U-equivalence flags. The non-BPS part is boxcompare_non.py."""
import sys, numpy as np
from cube import Det, T, act
from canon import canon
from bps import triples as btriples
from tricube import cube_from as bcube
B = int(sys.argv[1]); KMAX = int(sys.argv[2])
rng = np.arange(-B, B+1)
grid = np.array(np.meshgrid(*([rng]*8), indexing='ij')).reshape(8, -1).astype(np.int64)
a0,a1,a2,a3,a4,a5,a6,a7 = grid
Dv = (a0*a0*a7*a7 + a1*a1*a6*a6 + a2*a2*a5*a5 + a3*a3*a4*a4
      - 2*(a0*a1*a6*a7 + a0*a2*a5*a7 + a0*a3*a4*a7 + a1*a2*a5*a6 + a1*a3*a4*a6 + a2*a3*a4*a5)
      + 4*(a0*a3*a5*a6 + a1*a2*a4*a7))
J = ((0, 1), (1, 0))
for k in range(1, KMAX+1):
    ref = {}
    for x in btriples(k):
        a = bcube(k, x['J1'], x['J2'], x['I3'])
        for b in (a, act(act(a, 0, J), 1, J), act(act(a, 0, J), 2, J), act(act(a, 1, J), 2, J)):
            ref[canon(b)] = (x['fint'], x['equiv'] if x['fint'] else None)
    seen = {}
    for t in np.nonzero(Dv == -4*k*k)[0]:
        a = tuple(int(v) for v in grid[:, t]); c = canon(a)
        if c in seen: continue
        fint = all(v % (2*k) == 0 for v in T(a))
        same = canon(tuple(v // (2*k) for v in T(a))) == c if fint else None
        seen[c] = (fint, same)
    bad = sum(1 for c, v in seen.items() if ref.get(c) != v)
    print(f"BPS k={k}: orbits met in box {len(seen)} of {len(ref)}; mismatches {bad}", flush=True)
    assert bad == 0
print("box check passed")
