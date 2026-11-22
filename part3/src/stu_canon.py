"""STU model (from Part I): canonical form for SL2(Z)^3-orbits of definite cubes (Det < 0):
reduce the face forms Q1, Q2 to reduced forms, then take the SL2(Z)-Hermite normal form of the
remaining 2x4 slice matrix, minimised over the (finite) automorphism groups of Q1 and Q2."""
from cube import Det, T, form, act
import itertools

S_ = ((0, -1), (1, 0)); Sinv = ((0, 1), (-1, 0))
def Tn(n): return ((1, n), (0, 1))

def reduce_index(a, j):
    """Reduce the j-th face form to |B| <= C <= A (B >= 0 on the boundary); with our action,
    T^n changes B by 2nC and fixes C, S swaps A and C."""
    for _ in range(100000):
        A, B, C = form(a, j)
        s = 1 if A > 0 else -1
        A, B, C = s*A, s*B, s*C
        if abs(B) > C:
            n = -((B + C) // (2*C))  # bring B into (-C, C]
            cand = [act(a, j, Tn(n)), act(a, j, Tn(-n))]
            a = min(cand, key=lambda b: abs(form(b, j)[1])); continue
        if C > A:
            a = act(a, j, S_); continue
        if B < 0 and (abs(B) == C or A == C):
            for g in (S_, Tn(1), Tn(-1)):
                b = act(a, j, g); A2, B2, C2 = [s*x for x in form(b, j)]
                if abs(B2) <= C2 <= A2 and B2 >= 0:
                    a = b; break
            else:
                raise RuntimeError("boundary")
            continue
        return a
    raise RuntimeError("no convergence")

SMALL = [((p, q), (r, s)) for p, q, r, s in itertools.product(range(-2, 3), repeat=4) if p*s - q*r == 1]
def auts(a, j):
    Q = form(a, j)
    return [g for g in SMALL if form(act(a, j, g), j) == Q]

def hnf_sl2(rows):
    """rows: two integer rows of length 4 (rank 2). Return (H, sign) complete SL2(Z) invariant."""
    M = [list(rows[0]), list(rows[1])]
    U = [[1, 0], [0, 1]]
    # column-by-column gcd on first nonzero column
    col = 0
    while col < 4 and M[0][col] == 0 and M[1][col] == 0: col += 1
    while M[1][col] != 0:
        q = M[0][col] // M[1][col]
        M[0] = [x - q*y for x, y in zip(M[0], M[1])]; U[0] = [x - q*y for x, y in zip(U[0], U[1])]
        M[0], M[1] = M[1], M[0]; U[0], U[1] = U[1], U[0]
    if M[0][col] < 0:
        M[0] = [-x for x in M[0]]; U[0] = [-x for x in U[0]]
    c2 = col + 1
    while c2 < 4 and M[1][c2] == 0: c2 += 1
    if M[1][c2] < 0:
        M[1] = [-x for x in M[1]]; U[1] = [-x for x in U[1]]
    q = M[0][c2] // M[1][c2]
    M[0] = [x - q*y for x, y in zip(M[0], M[1])]; U[0] = [x - q*y for x, y in zip(U[0], U[1])]
    d = U[0][0]*U[1][1] - U[0][1]*U[1][0]
    return (tuple(M[0]), tuple(M[1]), d)

def slice_rows(a):
    """index-3 (j=2) rows: C = 0 and C = 1, entries ordered by (A, B)."""
    return [tuple(a[4*A + 2*B + C] for A in range(2) for B in range(2)) for C in range(2)]

def canon(a):
    assert Det(a) < 0
    b = reduce_index(reduce_index(a, 0), 1)
    best = None
    for u1 in auts(b, 0):
        c = act(b, 0, u1)
        for u2 in auts(c, 1):
            e = act(c, 1, u2)
            key = (form(e, 0), form(e, 1), hnf_sl2(slice_rows(e)))
            if best is None or key < best: best = key
    return best

if __name__ == "__main__":
    import random
    random.seed(7)
    gens = [S_, Tn(1), Tn(-1), ((2, 1), (1, 1))]
    n = 0
    while n < 300:
        a = tuple(random.randint(-4, 4) for _ in range(8))
        if Det(a) >= 0: continue
        b = a
        for _ in range(12):
            b = act(b, random.randrange(3), random.choice(gens))
        assert canon(a) == canon(b), (a, b)
        n += 1
    print("canonical form invariant under 300 random SL2(Z)^3 moves: OK")
