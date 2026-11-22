import sys
from math import gcd, isqrt
from cubic import *
from classes import classes_bps, classes_non, canon
def chi4(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def factor(n):
    out = {}; p = 2
    while p*p <= n:
        while n % p == 0: out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out
def nu(n): return sum(chi4(d) for d in range(1, n+1) if n % d == 0)
def dnum(n): return sum(1 for d in range(1, n+1) if n % d == 0)
def psi(f):
    r = f
    for p in factor(f): r = r * (p - chi4(p)) // p
    return r
def phi(f):
    r = f
    for p in factor(f): r = r * (p - 1) // p
    return r
def rho(f): return sum(1 for u in range(f) if (u*u + 1) % f == 0) if f > 1 else 1
def issq(n): return n >= 0 and isqrt(n)**2 == n
def Fp(k): return sum(psi(f)*nu(k//f**3) for f in range(1, k+1) if k % f**3 == 0)
def Ep(k): return nu(k)
def Fm(m):
    s = sum(phi(f)*dnum(m//f**3) for f in range(1, m+1) if m % f**3 == 0)
    s += sum(phi(f)*dnum(4*m//f**3) for f in range(2, 2*m+1, 2) if (4*m) % f**3 == 0 and (4*m//f**3) % 2 == 1)
    return s
def Em(m):
    s = sum(rho(f) for f in range(1, m+1) if m % f**3 == 0 and issq(m//f**3))
    s += sum(rho(f) for f in range(2, 2*m+1, 2) if (4*m) % f**3 == 0 and (4*m//f**3) % 2 == 1 and issq(4*m//f**3))
    return s
if __name__ == "__main__":
    K = int(sys.argv[1])
    print(" k  h   F  Fpred  E  Epred   | m  h   F  Fpred  E  Epred")
    for k in range(1, K+1):
        C = classes_bps(k); Fi = [f for f in C if dual(f) is not None]
        Ei = [f for f in Fi if canon(dual(f)) == f]
        Cn = classes_non(k); Fn = [f for f in Cn if dual(f) is not None]
        En = [f for f in Fn if canon(dual(f)) == f]
        print(f"{k:2d} {len(C):3d} {len(Fi):3d} {Fp(k):4d} {len(Ei):3d} {Ep(k):3d}    | {k:2d} {len(Cn):3d} {len(Fn):3d} {Fm(k):4d} {len(En):3d} {Em(k):3d}", flush=True)
