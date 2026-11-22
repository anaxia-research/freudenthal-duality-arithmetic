"""Forms produced from Bhargava data, compared with the classes found by Hessian reduction.
BPS (Delta = 4k^2): data (e, u, lam): e^3 | k, u in X_e, lam in Z[i] of norm k/e^3; f = Im(gamma xi^3)/k with gamma = k lam,
that is f = Re(alpha xi^3) and dual = Im(alpha xi^3) with alpha = -i lam as in the paper.
non-BPS (Delta = -m^2): data (e, t, A, B); f = (A xi_1^3 - B xi_2^3)/m on I_{e,t}."""
import sys
from math import gcd, isqrt
from fractions import Fraction as Fr
from cubic import *
from classes import classes_bps, classes_non, canon
from counts import Fp, Ep, Fm, Em

def cmul(z, w): return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])

def X(e):
    """Representatives u = (x, y) of (Z[i]/e)^x / (Z/e)^x, as Gaussian integers."""
    if e == 1: return [(1, 0)]
    seen = set(); reps = []
    for x in range(e):
        for y in range(e):
            n = (x*x + y*y) % e
            if gcd(n, e) != 1: continue
            key = min(((s*x) % e, (s*y) % e) for s in range(1, e) if gcd(s, e) == 1)
            if key in seen: continue
            seen.add(key); reps.append((x, y))
    return reps

def lattice_basis(gens):
    """Z-basis (Hermite form) of the lattice in Z^2 spanned by gens."""
    rows = [list(g) for g in gens]
    # column 1 (y) gcd first
    basis = []
    # simple approach: Euclid on y-coordinates then x
    import functools
    while True:
        rows = [r for r in rows if r != [0, 0]]
        ys = [r for r in rows if r[1] != 0]
        if len(ys) <= 1: break
        ys.sort(key=lambda r: abs(r[1]))
        p = ys[0]
        new = [p]
        for r in rows:
            if r is p: continue
            if r[1] != 0:
                q = r[1] // p[1]; new.append([r[0] - q*p[0], r[1] - q*p[1]])
            else: new.append(r)
        rows = new
    ys = [r for r in rows if r[1] != 0]; xs = [r for r in rows if r[1] == 0]
    g = 0
    for r in xs: g = gcd(g, r[0])
    v1 = ys[0]; v0 = [g, 0]
    return tuple(v0), tuple(v1)

def gen_bps(k):
    out = []
    for e in range(1, k + 1):
        if k % e**3: continue
        n = k // e**3
        alphas = [(x, y) for x in range(-isqrt(n), isqrt(n) + 1) for y in range(-isqrt(n), isqrt(n) + 1) if x*x + y*y == n and x > 0 and y >= 0]
        for u in X(e):
            l0, l1 = lattice_basis([u, (e, 0), (0, e)])
            if l0[0]*l1[1] - l1[0]*l0[1] < 0: l1 = (-l1[0], -l1[1])
            for al in alphas:
                gam = (k*al[0], k*al[1])
                def coef(p, q):
                    z = gam
                    for _ in range(p): z = cmul(z, l0)
                    for _ in range(q): z = cmul(z, l1)
                    return z
                zs = [coef(3, 0), coef(2, 1), coef(1, 2), coef(0, 3)]
                assert all(z[1] % k == 0 and z[0] % k == 0 for z in zs)
                f = tuple(z[1] // k for z in zs); ft = tuple(-z[0] // k for z in zs)
                out.append(((e, u, al), f, ft))
    return out

def phi(n): return sum(1 for t in range(1, n + 1) if gcd(t, n) == 1)

def gen_non(m):
    out = []
    for e in range(1, 2*m + 1):
        if m % e: continue
        units = [t for t in range(e) if gcd(t, e) == 1] if e > 1 else [0]
        pairs = []
        if m % e**3 == 0:
            n = m // e**3
            pairs += [(m*a, m*(n//a), 'a') for a in range(1, n + 1) if n % a == 0]
        if e % 2 == 0 and (4*m) % e**3 == 0 and (4*m // e**3) % 2 == 1:
            n = 4*m // e**3
            pairs += [(m*a//2, m*(n//a)//2, 'b') for a in range(1, n + 1) if n % a == 0]
        for t in units:
            l0, l1 = (0, e), (1, t)   # coordinate determinant -e < 0: positively oriented
            for A2, B2, case in pairs:
                A = Fr(A2, 1); B = Fr(B2, 1)  # delta^{-1} = (A, B) with A = m*a or m*a/2 (integers)
                def c(p, q):
                    x1 = l0[0]**p * l1[0]**q; x2 = l0[1]**p * l1[1]**q
                    return (A*x1 - B*x2) / m, (A*x1 + B*x2) / m
                cs = [c(3, 0), c(2, 1), c(1, 2), c(0, 3)]
                assert all(x.denominator == 1 and y.denominator == 1 for x, y in cs), (m, e, t, A, B)
                f = tuple(int(x) for x, y in cs); ft = tuple(int(y) for x, y in cs)
                out.append(((e, t, A2, B2, case), f, ft))
    return out

if __name__ == "__main__":
    K = int(sys.argv[1])
    for k in range(1, K + 1):
        G = gen_bps(k); C = classes_bps(k)
        Ci = {f for f in C if dual(f) is not None}
        cs = [canon(f) for _, f, _ in G]
        assert all(Delta(f) == 4*k*k for _, f, _ in G)
        assert all(dual(f) == ft for _, f, ft in G), k
        assert len(set(cs)) == len(cs) and set(cs) == Ci, (k, len(cs), len(set(cs)), len(Ci))
        for (e, u, al), f, ft in G:
            assert cont(hess(f)) == k // e
            assert (canon(ft) == canon(f)) == (e == 1)
        Gn = gen_non(k); Cn = classes_non(k)
        Cni = {f for f in Cn if dual(f) is not None}
        cn = [canon(f) for _, f, _ in Gn]
        assert all(Delta(f) == -k*k for _, f, _ in Gn), k
        assert all(dual(f) == ft for _, f, ft in Gn), k
        assert len(set(cn)) == len(cn) and set(cn) == Cni, (k, len(cn), len(set(cn)), len(Cni))
        for (e, t, A2, B2, case), f, ft in Gn:
            assert cont(hess(f)) == k // e
            g1 = cont([x + y for x, y in zip(f, ft)]); g2 = cont([x - y for x, y in zip(f, ft)])
            assert (g1, g2) == (2*A2 // k, 2*B2 // k), (k, e, t, A2, B2, g1, g2)
            U = canon(ft) == canon(f)
            assert U == ((t*t + 1) % e == 0 and A2 == B2), (k, e, t, A2, B2)
        print(f"k=m={k}: BPS {len(cs)} classes from data = enumeration; non-BPS {len(cn)} classes from data = enumeration", flush=True)
    print("generation check passed")
