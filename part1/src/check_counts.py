"""Compare the enumeration logs (k <= 40, m <= 40) with the closed formulas of Theorems E and F."""
import sys
from formula import Fplus
from math import comb
def Fminus(m):
    v = 1; n = m; p = 2
    while n > 1:
        if n % p == 0:
            a = 0
            while n % p == 0: n //= p; a += 1
            v *= comb(a+2, 2)*p**a if p == 2 else comb(a+2, 2)*p**a - comb(a+1, 2)*p**(a-1)
        p += 1
    return v
from asymptotics import Eminus_array
logdir = sys.argv[1] if len(sys.argv) > 1 else "../logs"
Em = Eminus_array(40)
nb = 0
for line in open(f"{logdir}/bps40.log"):
    t = line.split()
    if t and t[0].isdigit():
        k, h, T, P, F, Fe, G = map(int, t[:7])
        Ff, Gf = Fplus(k)
        assert (F, Fe, G) == (Ff, Ff - Gf, Gf), k
        nb += 1
nn = 0
seen = set()
for fn in ("split12.log", "split_more.log"):
    for line in open(f"{logdir}/{fn}"):
        t = line.split()
        if t and t[0].isdigit():
            if fn == "split12.log": m, T, P, F, Fe = map(int, t[:5])
            else: m, T, F, Fe = map(int, t[:4])
            assert (F, Fe) == (Fminus(m), int(Em[m])), m
            seen.add(m)
assert seen == set(range(1, 41))
print(f"enumeration agrees with Theorem E for k = 1..{nb} and with Theorem F for m = 1..40")
