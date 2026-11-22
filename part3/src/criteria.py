"""Theorems C and D and Corollary P on all orbits with integral dual, k, m <= K.
BPS (Delta = 4k^2): U-dual <=> cont(Q_S) = k <=> D_S = R_S/2k in SL2(Z), and then (D_S, 1) x = x~.
  Orbits with cont(Q_T) = k but cont(Q_S) < k: D_T on a single T index maps x to x~ (a U-duality of the STU model),
  but x~ is not in the SL2(Z)^2-orbit of x.
non-BPS (Delta = -m^2): U-dual <=> Q_S ~ -Q_S, Q_T ~ -Q_T and gcd(x + x~) = gcd(x - x~) <=> the matrices
  g_S, g_T built from the rank-one cubes x + x~ and x - x~ are integral; then (g_S, g_T) x = x~, g_S^2 = g_T^2 = -1.
Primitive Q_S and Q_T with integral dual occur only for Delta in {4, -1, -4}."""
import sys
from fractions import Fraction as Fr
from st2 import *

K = int(sys.argv[1])
neg = lambda Q: tuple(-t for t in Q)

def factor_vectors(R):
    """R = c lam (x) mu (x) mu (rank one, symmetric in the last two slots): primitive integral lam, mu."""
    x = cube(R)
    lam = mu = None
    for B in range(2):
        for C in range(2):
            v = (x[0*4 + 2*B + C], x[1*4 + 2*B + C])
            if v != (0, 0): lam = v; break
        if lam: break
    for A in range(2):
        for C in range(2):
            v = (x[4*A + 0 + C], x[4*A + 2 + C])
            if v != (0, 0): mu = v; break
        if mu: break
    g1 = cont(lam); g2 = cont(mu)
    return (lam[0]//g1, lam[1]//g1), (mu[0]//g2, mu[1]//g2)

def swap_matrix(c1, c2, p):
    """The rational g with g c1 = p c2 and g c2 = -c1/p (columns c1, c2); det g = 1, g^2 = -1."""
    M = ((Fr(c1[0]), Fr(c2[0])), (Fr(c1[1]), Fr(c2[1])))
    d = M[0][0]*M[1][1] - M[0][1]*M[1][0]
    Mi = ((M[1][1]/d, -M[0][1]/d), (-M[1][0]/d, M[0][0]/d))
    N = ((p*c2[0], -c1[0]/p), (p*c2[1], -c1[1]/p))
    return tuple(tuple(sum(N[i][k]*Mi[k][j] for k in range(2)) for j in range(2)) for i in range(2))

def integral(g): return all(t.denominator == 1 for row in g for t in row)
def toint(g): return tuple(tuple(int(t) for t in row) for row in g)

nb = nbU = nT = 0; nn = nnU = 0; prim = set()
for k in range(1, K + 1):
    for s in orbits(4*k*k):
        d = dual(s)
        if d is None: continue
        nb += 1
        U = canon(d) == s
        assert U == (cont(QS(s)) == k)
        R = RS(s)
        if U:
            nbU += 1
            D = tuple(tuple(t // (2*k) for t in row) for row in R)
            assert all(t % (2*k) == 0 for row in R for t in row) and det2(D) == 1
            assert act(D, I2, s) == d and mul(D, D) == MI
        else:
            assert any(t % (2*k) for row in R for t in row)
        if cont(QT(s)) == k and not U:
            nT += 1
            RTm = RT(s); DT = tuple(tuple(t // (2*k) for t in row) for row in RTm)
            assert all(t % (2*k) == 0 for row in RTm for t in row)
            assert act1(DT, s, 1) == cube(d) and act1(DT, s, 2) == cube(d)
            assert act(I2, DT, s) == tuple(-t for t in s)
        if cont(QS(s)) == 1 and cont(QT(s)) == 1: prim.add(4*k*k)
for m in range(1, K + 1):
    for s in orbits(-m*m):
        d = dual(s)
        if d is None: continue
        nn += 1
        U = canon(d) == s
        crit = (qcanon(QS(s)) == qcanon(neg(QS(s))) and qcanon(QT(s)) == qcanon(neg(QT(s)))
                and cont([p + q for p, q in zip(s, d)]) == cont([p - q for p, q in zip(s, d)]))
        assert U == crit, (m, s)
        P = tuple(p + q for p, q in zip(s, d)); N = tuple(p - q for p, q in zip(s, d))
        l1, m1 = factor_vectors(P); l2, m2 = factor_vectors(N)
        found = None
        Ac, Bc = cont(P), cont(N)
        for p in (Fr(Ac, Bc), Fr(-Ac, Bc), Fr(Bc, Ac), Fr(-Bc, Ac)):
            gS = swap_matrix(l1, l2, p); gT = swap_matrix(m1, m2, Fr(1))
            # check on rationals: (gS, gT) x = x~
            x = cube(s); y = [Fr(0)]*8
            for A in range(2):
                for B in range(2):
                    for C in range(2):
                        y[4*A + 2*B + C] = sum(gS[A][a]*gT[B][b]*gT[C][c]*x[4*a + 2*b + c]
                                               for a in range(2) for b in range(2) for c in range(2))
            if tuple(y) == tuple(Fr(t) for t in cube(d)): found = (gS, gT)
        assert found is not None, (m, s)
        gS, gT = found
        assert (integral(gS) and integral(gT)) == U, (m, s, gS, gT, U)
        if U:
            nnU += 1
            a, b = toint(gS), toint(gT)
            assert act(a, b, s) == d and mul(a, a) == MI and mul(b, b) == MI
        if cont(QS(s)) == 1 and cont(QT(s)) == 1: prim.add(-m*m)
print("BPS: %d orbits with integral dual (k <= %d); U-dual exactly when cont(Q_S) = k (%d orbits), and then "
      "D_S = R_S/2k lies in SL2(Z), D_S^2 = -1, (D_S, 1) x = x~" % (nb, K, nbU))
print("BPS: %d orbits have cont(Q_T) = k but are not U-dual; for them D_T on one T index maps x to x~ "
      "(a U-duality of the STU model) and (1, D_T) x = -x" % nT)
print("non-BPS: %d orbits with integral dual (m <= %d); U-dual exactly when Q_S ~ -Q_S, Q_T ~ -Q_T and "
      "gcd(x+x~) = gcd(x-x~) (%d orbits)" % (nn, K, nnU))
print("non-BPS: the real duality (g_S, g_T) with g_T exchanging the T-factors of x + x~ and x - x~ up to sign is found "
      "for every orbit; it is integral exactly for the U-dual orbits, and then g_S^2 = g_T^2 = -1 and (g_S, g_T) x = x~")
print("primitive Q_S and Q_T with integral dual occur only for Delta in", sorted(prim))
