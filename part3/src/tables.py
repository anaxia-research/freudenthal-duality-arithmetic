"""Enumeration of all SL2(Z)^2-orbits with Delta = 4k^2 and Delta = -m^2 (k, m <= K), integral duals and U-duality,
compared with the closed formulas.  Writes logs/table.txt (K = 40) or logs/table<K>.txt."""
import sys, time
from st2 import *
from counts import Fp, Ep, hp, Fm, Em, hm

K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
t0 = time.time()
rows = []
totF = totE = 0
for n in range(1, K + 1):
    O = orbits(4*n*n)
    Fl = [s for s in O if dual(s) is not None]
    El = [s for s in Fl if canon(dual(s)) == s]
    pos = lambda L: sum(1 for s in L if QT(s)[0] > 0)
    assert 2*pos(O) == len(O) and 2*pos(Fl) == len(Fl) and 2*pos(El) == len(El), n
    hP, FP, EP = len(O), len(Fl), len(El)
    O = orbits(-n*n)
    Fl = [s for s in O if dual(s) is not None]
    El = [s for s in Fl if canon(dual(s)) == s]
    hM, FM, EM = len(O), len(Fl), len(El)
    assert (hP, FP, EP) == (hp(n), Fp(n), Ep(n)), (n, hP, FP, EP, hp(n), Fp(n), Ep(n))
    assert (hM, FM, EM) == (hm(n), Fm(n), Em(n)), (n, hM, FM, EM, hm(n), Fm(n), Em(n))
    rows.append((n, hP, FP, EP, hM, FM, EM)); totF += FP + FM; totE += EP + EM
    if n % 10 == 0 or n == K: print("n <= %d checked (%.0f s)" % (n, time.time() - t0), flush=True)
print("h_+, F_+, E_+, h_-, F_-, E_- agree with the closed formulas for all n <=", K)
print("each sign of Q_T carries exactly half of the orbits with Delta > 0")
print("orbits with integral dual:", totF, " U-dual:", totE)
with open('../logs/table%s.txt' % ('' if K == 40 else K), 'w') as fh:
    fh.write("n h+ F+ E+ h- F- E-\n")
    for r in rows: fh.write(" ".join(map(str, r)) + "\n")
