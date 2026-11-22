"""T^3 model (from Part II): SL2(Z)-classes of classical integral binary cubic forms with square |Delta|, by reduction of the Hessian.
canon(f) is a canonical representative of the SL2(Z)-orbit of f."""
import sys
from math import isqrt, gcd
from t3_cubic import *

def egcd(a, b):
    if b == 0: return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)

def reduce_def(Q):
    """Positive definite Q: return (Qred, M) with qsub(Q, M) == Qred reduced, M in SL2(Z)."""
    M = ((1, 0), (0, 1))
    while True:
        P, q, R = Q
        t = -((q + P) // (2 * P))  # makes -P <= q + 2tP < P
        if t:
            Tt = ((1, t), (0, 1)); Q = qsub(Q, Tt); M = mul(M, Tt); P, q, R = Q
        if P > R:
            S = ((0, -1), (1, 0)); Q = qsub(Q, S); M = mul(M, S); continue
        break
    P, q, R = Q
    if q == -P:
        Tt = ((1, 1), (0, 1)); Q = qsub(Q, Tt); M = mul(M, Tt)
    P, q, R = Q
    if P == R and q < 0:
        S = ((0, -1), (1, 0)); Q = qsub(Q, S); M = mul(M, S)
    return Q, M

SMALL = [((a, b), (c, d)) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) for d in (-1, 0, 1) if a*d - b*c == 1]
def aut_def(Q): return [A for A in SMALL if qsub(Q, A) == Q]

def canon_bps(f):
    Q, M = reduce_def(hess(f))
    f1 = act(phi_inv(M), f)
    assert hess(f1) == Q
    return min(act(phi_inv(A), f1) for A in aut_def(Q))

def roots_vectors(Q, m):
    P, q, R = Q
    out = []
    if P == 0:
        out.append((1, 0)); u, v = -R, q
    else:
        for s in (1, -1):
            u, v = -q + s*m, 2*P
            g = gcd(u, v); out.append((u//g, v//g))
        return out
    g = gcd(u, v); out.append((u//g, v//g)); return out

def canon_non(f):
    D = Det(f); m = isqrt(D); assert m*m == D and m > 0
    Q = hess(f); cands = []
    for (u0, v0) in roots_vectors(Q, m):
        g, x, y = egcd(u0, v0)   # x u0 + y v0 = 1
        M = ((u0, -y), (v0, x))
        Q1 = qsub(Q, M); assert Q1[0] == 0, (Q, M, Q1)
        q1, R1 = Q1[1], Q1[2]
        t = -(R1 // abs(q1)) * (1 if q1 > 0 else -1)
        # want 0 <= R1 + t q1 < |q1|
        t = (-(R1 // q1)) if q1 > 0 else ((R1 // (-q1)))
        Tt = ((1, t), (0, 1)); M2 = mul(M, Tt); Q2 = qsub(Q, M2)
        assert Q2[0] == 0 and 0 <= Q2[2] < abs(Q2[1]), (Q, Q2)
        f2 = act(phi_inv(M2), f)
        cands.append(f2); cands.append(tuple(-x for x in f2))
    return min(cands)

def canon(f): return canon_bps(f) if Delta(f) > 0 else canon_non(f)

def classes_bps(k):
    """All SL2(Z)-classes with Delta = 4k^2 (Det = -4k^2)."""
    D = -4*k*k; found = set()
    Pmax = isqrt(-D // 3) + 1
    for P in range(1, Pmax + 1):
        for q in range(-P + 1, P + 1):
            if (q*q - D) % (4*P): continue
            R = (q*q - D) // (4*P)
            if R < P or (P == R and q < 0): continue
            Q = (P, q, R)
            amax = isqrt(P**3 // (k*k)) + 1
            for a in range(-amax, amax + 1):
                if a == 0:
                    b0 = isqrt(P)
                    if b0*b0 != P: continue
                    for b in (b0, -b0):
                        if q % b: continue
                        c = -q // b
                        if (c*c - R) % b: continue
                        d = (c*c - R) // b
                        f = (0, b, c, d)
                        if hess(f) == Q: found.add(canon_bps(f))
                    continue
                disc = 4*P**3 - 4*k*k*a*a
                if disc < 0: continue
                s = isqrt(disc)
                if s*s != disc: continue
                for num in (-q*a + s, -q*a - s):
                    if num % (2*P): continue
                    b = num // (2*P)
                    if (b*b - P) % a: continue
                    c = (b*b - P) // a
                    if (q + b*c) % a: continue
                    d = (q + b*c) // a
                    f = (a, b, c, d)
                    if hess(f) == Q: found.add(canon_bps(f))
    return found

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def classes_non(m):
    """All SL2(Z)-classes with Delta = -m^2 (Det = m^2)."""
    found = set()
    for q in (m, -m):
        for R in range(m):
            Q = (0, q, R)
            for a0 in divisors(m*m):
                for a in (a0, -a0):
                    if (R*a) % q: continue
                    b = -R*a // q
                    if (b*b) % a: continue
                    c = b*b // a
                    if (q + b*c) % a: continue
                    d = (q + b*c) // a
                    f = (a, b, c, d)
                    if hess(f) == Q: found.add(canon_non(f))
    return found

if __name__ == "__main__":
    import random
    random.seed(7)
    # canonical forms are SL2(Z)-invariant
    G = [((1,1),(0,1)), ((0,-1),(1,0)), ((2,1),(1,1)), ((1,0),(3,1)), ((1,-2),(1,-1))]
    n = 0
    for _ in range(3000):
        f = tuple(random.randint(-9, 9) for _ in range(4))
        D = Delta(f)
        if D == 0: continue
        if D < 0 and isqrt(-D)**2 != -D: continue
        g = ((1, 0), (0, 1))
        for _ in range(4): g = mul(g, random.choice(G))
        assert canon(f) == canon(act(g, f)), (f, g)
        n += 1
    print("canonical form invariance checked on", n, "forms")
