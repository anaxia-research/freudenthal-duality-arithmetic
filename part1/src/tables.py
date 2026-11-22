"""Plain-text tables of the counts from the closed formulas, cross-checked against the enumeration logs."""
import re, sys
from formula import Fplus
from asymptotics import Eminus_array
from series import arrays

def chi(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def Fminus(m):
    from math import comb
    v = 1; n = m; p = 2
    while n > 1:
        if n % p == 0:
            a = 0
            while n % p == 0: n //= p; a += 1
            v *= comb(a+2, 2)*p**a if p == 2 else comb(a+2, 2)*p**a - comb(a+1, 2)*p**(a-1)
        p += 1
    return v
Em = Eminus_array(40)
# enumeration logs
bps = {}
for line in open("../logs/bps40.log"):
    t = line.split()
    if t and t[0].isdigit():
        k, h, T, P, F, Fe, G = map(int, t[:7]); bps[k] = (h, T, F, Fe, G)
non = {}
for line in open("../logs/split12.log"):
    t = line.split()
    if t and t[0].isdigit():
        m, T, P, F, Fe = map(int, t[:5]); non[m] = (T, F, Fe)
for line in open("../logs/split_more.log"):
    t = line.split()
    if t and t[0].isdigit():
        m, T, F, Fe = map(int, t[:4]); non[m] = (T, F, Fe)
rows = []
for k in range(1, 25):
    F, G = Fplus(k)
    h, T, Fe_, Fe, Ge = bps[k]
    assert (F, F - G, G) == (Fe_, Fe, Ge), k
    rows.append(f"{k:3d} {4*k*k:6d} {h:4d} {T:6d} {F:5d} {F-G:5d} {G:4d}")
open("../logs/table_bps.txt", "w").write("  k  Delta    h     T+    F+    E+   G+   (positive definite face forms; full counts are 4x)\n" + "\n".join(rows) + "\n")
rows = []
for m in range(1, 33):
    F = Fminus(m); E = int(Em[m])
    if m in non:
        T, Fn, En = non[m]; assert (Fn, En) == (F, E), m
        Ts = str(T)
    else:
        Ts = "--"
    if m <= 24 or m in (25, 32):
        rows.append(f"{m:3d} {-m*m:6d} {Ts:>6} {F:6d} {E:4d}")
open("../logs/table_non.txt", "w").write("  m  Delta     T-     F-   E-\n" + "\n".join(rows) + "\n")
print("tables written to logs/; all enumerated values agree with the closed formulas")
