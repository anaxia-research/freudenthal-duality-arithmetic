"""Parallel version of the non-BPS half of the box completeness check (same test as boxcompare.py)."""
import sys, numpy as np
from multiprocessing import Pool
from cube import T
from canonsq import canon_sq
from split import triples as striples
from splitcube import cube_from as scube
def work(a): return canon_sq(a), a
if __name__ == "__main__":
    B = int(sys.argv[1]); M0 = int(sys.argv[2]); M1 = int(sys.argv[3])
    rng = np.arange(-B, B+1)
    grid = np.array(np.meshgrid(*([rng]*8), indexing='ij')).reshape(8, -1).astype(np.int64)
    a0,a1,a2,a3,a4,a5,a6,a7 = grid
    Dv = (a0*a0*a7*a7 + a1*a1*a6*a6 + a2*a2*a5*a5 + a3*a3*a4*a4
          - 2*(a0*a1*a6*a7 + a0*a2*a5*a7 + a0*a3*a4*a7 + a1*a2*a5*a6 + a1*a3*a4*a6 + a2*a3*a4*a5)
          + 4*(a0*a3*a5*a6 + a1*a2*a4*a7))
    with Pool(4) as pool:
        for m in range(M0, M1+1):
            ref = {}
            for x in striples(m):
                a = scube(m, x['J1'], x['J2'], x['I3']); ref[canon_sq(a)] = x['fint']
            cubes = [tuple(int(v) for v in grid[:, t]) for t in np.nonzero(Dv == m*m)[0]]
            seen = {}
            for c, a in pool.imap_unordered(work, cubes, chunksize=200):
                if c not in seen: seen[c] = all(v % m == 0 for v in T(a))
            bad = sum(1 for c, v in seen.items() if ref.get(c) != v)
            print(f"non-BPS m={m}: cubes {len(cubes)}, orbits met in box {len(seen)} of {len(ref)}; mismatches {bad}", flush=True)
            assert bad == 0
    print("box check passed")
