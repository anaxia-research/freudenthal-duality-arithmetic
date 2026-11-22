"""Bhargava-free canonical form for SL2(Z)^3-orbits of cubes with square positive discriminant."""
from cube import Det, T, form, act
from stu_canon import hnf_sl2, slice_rows, Tn
from math import gcd, isqrt
import itertools

def egcd(a, b):
    if b == 0: return (a, 1, 0)
    g, x, y = egcd(b, a % b); return (g, y, x - (a // b)*y)

def isotropic(Q, m):
    A, B, C = Q
    out = []
    if A != 0:
        for s in (1, -1):
            x, y = -B + s*m, 2*A
            g = gcd(abs(x), abs(y)); out.append((x//g, y//g))
    else:
        out.append((1, 0))
        g = gcd(abs(C), abs(B)); out.append((C//g, -B//g))
    return out

def reductions(a, j, m):
    res = []
    for v in isotropic(form(a, j), m):
        v1, v2 = v
        # Q_new(x, y) = Q(x g00 - y g10, -x g01 + y g11); C_new = Q(-g10, g11) = 0 for (g10, g11) = (-v1, v2)
        g0, p, q = egcd(v2, v1)            # p v2 + q v1 = g0 = +-1
        if g0 < 0: p, q = -p, -q
        assert p*v2 + q*v1 == 1
        g = ((p, q), (-v1, v2))
        b = act(a, j, g)
        A, B, C = form(b, j)
        assert C == 0 and abs(B) == m, (form(a, j), v, (A, B, C))
        n = A // B                           # A - n B in [0, |B|) if B > 0; handle sign
        b2 = act(b, j, Tn(n))
        A2, B2, C2 = form(b2, j)
        if not (0 <= A2 < abs(B2)):
            for nn in (n-1, n+1, -n, -n-1, -n+1):
                b3 = act(b, j, Tn(nn)); A3, B3, _ = form(b3, j)
                if 0 <= A3 < abs(B3): b2 = b3; break
        A2, B2, C2 = form(b2, j)
        assert C2 == 0 and 0 <= A2 < abs(B2)
        res.append(b2)
    return res

SMALL = [((p, q), (r, s)) for p, q, r, s in itertools.product(range(-3, 4), repeat=4) if p*s - q*r == 1]
def auts(a, j):
    Q = form(a, j)
    return [g for g in SMALL if form(act(a, j, g), j) == Q]

def canon_sq(a):
    D = Det(a); m = isqrt(D); assert m*m == D and m > 0
    best = None
    for b in reductions(a, 0, m):
        for u1 in auts(b, 0):
            c = act(b, 0, u1)
            for d in reductions(c, 1, m):
                for u2 in auts(d, 1):
                    e = act(d, 1, u2)
                    key = (form(e, 0), form(e, 1), hnf_sl2(slice_rows(e)))
                    if best is None or key < best: best = key
    return best

if __name__ == "__main__":
    import random
    random.seed(3)
    S_ = ((0, -1), (1, 0))
    gens = [S_, Tn(1), Tn(-1), ((2, 1), (1, 1))]
    n = 0
    while n < 300:
        a = tuple(random.randint(-3, 3) for _ in range(8))
        D = Det(a)
        if D <= 0 or isqrt(D)**2 != D: continue
        b = a
        for _ in range(12): b = act(b, random.randrange(3), random.choice(gens))
        assert canon_sq(a) == canon_sq(b), (a, b)
        n += 1
    print("square-discriminant canonical form invariant under 300 random moves: OK")
