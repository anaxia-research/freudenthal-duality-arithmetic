"""Tables of class numbers and Freudenthal counts (enumeration by Hessian reduction), k, m <= K."""
import sys
from cubic import *
from classes import classes_bps, classes_non, canon
from counts import Fp, Ep, Fm, Em
K = int(sys.argv[1])
rows = []
for k in range(1, K + 1):
    C = classes_bps(k); Fi = [f for f in C if dual(f) is not None]; Ei = [f for f in Fi if canon(dual(f)) == f]
    Cn = classes_non(k); Fn = [f for f in Cn if dual(f) is not None]; En = [f for f in Fn if canon(dual(f)) == f]
    assert (len(Fi), len(Ei), len(Fn), len(En)) == (Fp(k), Ep(k), Fm(k), Em(k)), k
    rows.append((k, len(C), len(Fi), len(Ei), len(Cn), len(Fn), len(En)))
    print(*rows[-1], flush=True)
with open('../logs/table%s.txt' % ('' if K == 40 else K), 'w') as fh:
    fh.write("k  h+(k) F+(k) E+(k) | m  h-(m) F-(m) E-(m)   (h = all SL2(Z)-classes, F = integral dual, E = U-dual)\n")
    for r in rows: fh.write("%2d %5d %5d %5d | %2d %5d %5d %5d\n" % (r[0], r[1], r[2], r[3], r[0], r[4], r[5], r[6]))
print("enumeration agrees with Theorem E for k, m <=", K)
