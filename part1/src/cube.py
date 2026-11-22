"""Bhargava cubes as STU black hole charges: hyperdeterminant, triple product, quadratic forms."""
from math import gcd
import itertools, random

def idx(A, B, C): return 4*A + 2*B + C

def Det(a):
    a0,a1,a2,a3,a4,a5,a6,a7 = a
    return (a0*a0*a7*a7 + a1*a1*a6*a6 + a2*a2*a5*a5 + a3*a3*a4*a4
            - 2*(a0*a1*a6*a7 + a0*a2*a5*a7 + a0*a3*a4*a7 + a1*a2*a5*a6 + a1*a3*a4*a6 + a2*a3*a4*a5)
            + 4*(a0*a3*a5*a6 + a1*a2*a4*a7))

def T(a):
    a0,a1,a2,a3,a4,a5,a6,a7 = a
    return (a0*( a3*a4 + a2*a5 + a1*a6 - a0*a7) - 2*a1*a2*a4,
            a1*(-a3*a4 - a2*a5 + a1*a6 - a0*a7) + 2*a0*a3*a5,
            a2*(-a3*a4 + a2*a5 - a1*a6 - a0*a7) + 2*a0*a3*a6,
            a3*(-a3*a4 + a2*a5 + a1*a6 + a0*a7) - 2*a1*a2*a7,
            a4*( a3*a4 - a2*a5 - a1*a6 - a0*a7) + 2*a0*a5*a6,
            a5*( a3*a4 - a2*a5 + a1*a6 + a0*a7) - 2*a1*a4*a7,
            a6*( a3*a4 + a2*a5 - a1*a6 + a0*a7) - 2*a2*a4*a7,
            a7*(-a3*a4 - a2*a5 - a1*a6 + a0*a7) + 2*a3*a5*a6)

def slice_(a, i, s):
    """2x2 matrix of the face with index i fixed to s (rows/cols = remaining indices in order)."""
    m = [[0,0],[0,0]]
    for u in range(2):
        for v in range(2):
            ind = [0,0,0]; ind[i] = s
            rest = [j for j in range(3) if j != i]
            ind[rest[0]] = u; ind[rest[1]] = v
            m[u][v] = a[idx(*ind)]
    return m

def form(a, i):
    """Q_i(x,y) = -det(x M_i - y N_i) as (A,B,C) with Q = A x^2 + B x y + C y^2."""
    M = slice_(a, i, 0); N = slice_(a, i, 1)
    det = lambda m: m[0][0]*m[1][1] - m[0][1]*m[1][0]
    A = -det(M); C = -det(N)
    # -det(xM - yN) = -(x^2 det M - x y (M00 N11 + M11 N00 - M01 N10 - M10 N01) + y^2 det N)
    B = (M[0][0]*N[1][1] + M[1][1]*N[0][0] - M[0][1]*N[1][0] - M[1][0]*N[0][1])
    return (A, B, C)

def act(a, i, g):
    """Act by g in SL2(Z) on the i-th index."""
    out = [0]*8
    for A in range(2):
        for B in range(2):
            for C in range(2):
                ind = [A,B,C]
                s = 0
                for t in range(2):
                    j = list(ind); j[i] = t
                    s += g[ind[i]][t] * a[idx(*j)]
                out[idx(A,B,C)] = s
    return tuple(out)

if __name__ == "__main__":
    random.seed(1)
    bad = 0
    for _ in range(2000):
        a = tuple(random.randint(-6, 6) for _ in range(8))
        D = Det(a)
        t = T(a)
        assert Det(t) == D**3, (a, Det(t), D)
        tt = T(t)
        assert all(tt[j] == -D*D*a[j] for j in range(8)) or all(tt[j] == D*D*a[j] for j in range(8)), (a, tt, D)
        sgn = 1 if all(tt[j] == D*D*a[j] for j in range(8)) else -1
        for i in range(3):
            A, B, C = form(a, i)
            assert B*B - 4*A*C == D
        c = [-(t[j] + (D % 4)*a[j]) for j in range(8)]
        assert all(x % 2 == 0 for x in c)
        for i in range(3):
            g = random.choice([((1,1),(0,1)), ((0,-1),(1,0)), ((2,1),(1,1))])
            b = act(a, i, g)
            assert Det(b) == D
            tb = T(b); ta = act(T(a), i, g)
            if tb != ta: bad += 1
    print("Det(T)=Det^3, T(T(a)) = %+d Det^2 a, forms disc = Det, c integral; T equivariance failures:" % sgn, bad)
