"""Theorems C and D on every class with integral Freudenthal dual, k, m <= K (enumeration by Hessian reduction)."""
import sys
from math import gcd, isqrt
from cubic import *
from classes import classes_bps, classes_non, canon, roots_vectors, egcd
def qcanon(Q, m):
    best = None
    for (u0, v0) in roots_vectors(Q, m):
        g, x, y = egcd(u0, v0); M = ((u0, -y), (v0, x)); Q1 = qsub(Q, M)
        q1, R1 = Q1[1], Q1[2]
        t = (-(R1 // q1)) if q1 > 0 else (R1 // (-q1))
        Q2 = qsub(Q, mul(M, ((1, t), (0, 1))))
        best = Q2 if best is None or Q2 < best else best
    return best
def rational_root(f):
    a, b, c, d = f
    if a == 0 or d == 0: return True
    for p in range(1, abs(d) + 1):
        if d % p: continue
        for q in range(1, abs(a) + 1):
            if a % q: continue
            for s in (1, -1):
                if a*(s*p)**3 + 3*b*(s*p)**2*q + 3*c*(s*p)*q*q + d*q**3 == 0: return True
    return False
def cube_root_form(P):
    """P = A (p u + q v)^3 with A > 0 and (p, q) primitive: return (A, p, q)."""
    A = cont(P); x = [c // A for c in P]
    # x = (p^3, p^2 q, p q^2, q^3) up to a common sign
    for s in (1, -1):
        y = [s*c for c in x]
        p = round(abs(y[0]) ** (1/3)) * (1 if y[0] >= 0 else -1); q = round(abs(y[3]) ** (1/3)) * (1 if y[3] >= 0 else -1)
        if (p**3, p*p*q, p*q*q, q**3) == tuple(y): return s*A, p, q
    raise ValueError(P)
def g_of(f, ft):
    """the real U-duality g_f with g_f . f = ft when A = B (Theorem D); None if A != B"""
    A, p1, q1 = cube_root_form([x + y for x, y in zip(f, ft)])
    B, p2, q2 = cube_root_form([y - x for x, y in zip(f, ft)])
    if A < 0: A, p1, q1 = -A, -p1, -q1
    if B < 0: B, p2, q2 = -B, -p2, -q2
    if A != B: return None
    # g [l1 l2] = [l2, -l1]
    det = p1*q2 - q1*p2
    M = ((p2, -p1), (q2, -q1)); Ninv = ((q2, -p2), (-q1, p1))  # [l1 l2]^{-1} * det
    num = tuple(tuple(sum(M[i][k]*Ninv[k][j] for k in range(2)) for j in range(2)) for i in range(2))
    return num, det
K = int(sys.argv[1]); nb = nn = 0; ng = 0; proj = []
for k in range(1, K + 1):
    for f in classes_bps(k):
        ft = dual(f)
        if ft is None: continue
        U = canon(ft) == f
        p, q, r = hess(f)
        Dint = q % (2*k) == 0 and p % k == 0 and r % k == 0
        assert U == (cont(hess(f)) == k) == Dint, f
        if U: assert act(((q//(2*k), 2*p//(2*k)), (-2*r//(2*k), -q//(2*k))), f) == ft
        if cont(hess(f)) == 1: proj.append(4*k*k)
        nb += 1
    m = k
    for f in classes_non(m):
        ft = dual(f)
        if ft is None: continue
        U = canon(ft) == f
        Q = hess(f); h = (Q[0], -Q[1], Q[2])  # hat Q as a form (P, B, C) with B = bc - ad
        # proper equivalence of h with -h, using the substitution action on the coefficient triple
        hc = qcanon(Q, m); hnc = qcanon(tuple(-x for x in Q), m)
        g1 = cont([x + y for x, y in zip(f, ft)]); g2 = cont([x - y for x, y in zip(f, ft)])
        assert U == (hc == hnc and g1 == g2), f
        if U: assert rational_root(f)
        gd = g_of(f, ft)
        integral = gd is not None and all(x % gd[1] == 0 for row in gd[0] for x in row)
        assert integral == U, f
        if U:
            g = tuple(tuple(x // gd[1] for x in row) for row in gd[0])
            assert g[0][0]*g[1][1] - g[0][1]*g[1][0] == 1 and g[0][0] + g[1][1] == 0 and act(g, f) == ft, f
            ng += 1
        if cont(hess(f)) == 1: proj.append(-m*m)
        nn += 1
print("Theorem C checked on", nb, "BPS classes and Theorem D on", nn, "non-BPS classes with integral dual, k, m <=", K)
print("g_f integral exactly for the U-dual non-BPS classes; g_f . f = ft, g_f^2 = -1 checked on", ng, "classes")
print("projective classes with integral dual occur only for Delta in", sorted(set(proj)))
