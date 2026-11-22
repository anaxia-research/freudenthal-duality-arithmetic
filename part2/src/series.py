"""Dirichlet series of the counting functions (Theorem F), checked coefficientwise up to N."""
import sys
from math import gcd
from counts import Fp, Ep, Fm, Em, nu, dnum, psi, phi, rho, chi4, factor
N = int(sys.argv[1])
def conv(A, B):
    C = [0]*(N+1)
    for i in range(1, N+1):
        if A[i] == 0: continue
        for j in range(1, N//i + 1): C[i*j] += A[i]*B[j]
    return C
def euler(local):
    """multiplicative function from local(p, e)."""
    out = [0]*(N+1); out[1] = 1
    for n in range(2, N+1):
        r = 1
        for p, e in factor(n).items(): r *= local(p, e)
        out[n] = r
    return out
one = [0] + [1]*N
chi = [0] + [chi4(n) for n in range(1, N+1)]
zL = conv(one, chi)                                     # zeta(s) L(s,chi)
cubeP = [0]*(N+1)                                       # zeta(3s-1)/L(3s,chi) = sum psi(e) e^{-3s}
e = 1
while e**3 <= N: cubeP[e**3] = psi(e); e += 1
FP = conv(zL, cubeP)
ok1 = all(FP[n] == Fp(n) for n in range(1, N+1)) and all(zL[n] == Ep(n) for n in range(1, N+1))
# non-BPS: zeta(s)^2 zeta(3s-1)/zeta(3s) [1 + 2^{-s}(1-2^{-s})^2/(1-2^{-3s})]
d2 = conv(one, one)
cubePhi = [0]*(N+1); e = 1
while e**3 <= N: cubePhi[e**3] = phi(e); e += 1
base = conv(d2, cubePhi)
# correction series c(s) = 1 + 2^{-s}(1-2^{-s})^2 / (1-2^{-3s}): coefficients at powers of 2
corr = [0]*(N+1); corr[1] = 1
# 2^{-s}(1 - 2*2^{-s} + 2^{-2s}) * sum_j 2^{-3js}
j = 0
while 2**(3*j+1) <= N:
    for ex, cf in ((1, 1), (2, -2), (3, 1)):
        n = 2**(3*j + ex)
        if n <= N: corr[n] += cf
    j += 1
FM = conv(base, corr)
ok2 = all(FM[n] == Fm(n) for n in range(1, N+1))
# E_-: zeta(2s) zeta(3s) L(3s,chi)/zeta(6s) [1 + 2^{-s}(1-2^{-2s})/(1+2^{-3s})]
sq = [0]*(N+1); i = 1
while i*i <= N: sq[i*i] = 1; i += 1
rho3 = [0]*(N+1); e = 1
while e**3 <= N: rho3[e**3] = rho(e); e += 1
baseE = conv(sq, rho3)
corrE = [0]*(N+1); corrE[1] = 1
# 2^{-s}(1 - 2^{-2s}) * sum_j (-1)^j 2^{-3js}
j = 0
while 2**(3*j+1) <= N:
    for ex, cf in ((1, 1), (3, -1)):
        n = 2**(3*j + ex)
        if n <= N: corrE[n] += cf*(-1)**j
    j += 1
EM = conv(baseE, corrE)
ok3 = all(EM[n] == Em(n) for n in range(1, N+1))
# rho series: sum rho(n) n^{-s} = zeta(s) L(s,chi) / zeta(2s)
mu2 = [0]*(N+1)
for n in range(1, int(N**0.5)+1):
    fz = factor(n); mu2[n*n] = 0 if any(v > 1 for v in fz.values()) else (-1)**len(fz)
ok4 = all(conv(zL, mu2)[n] == rho(n) for n in range(1, N+1))
print("F_+ = zeta L(chi) zeta(3s-1)/L(3s,chi), E_+ = zeta L(chi):", ok1)
print("F_- = zeta^2 zeta(3s-1)/zeta(3s) (1 + 2^-s (1-2^-s)^2/(1-2^-3s)):", ok2)
print("E_- = zeta(2s) zeta(3s) L(3s,chi)/zeta(6s) (1 + 2^-s (1-2^-2s)/(1+2^-3s)):", ok3)
print("sum rho(n) n^-s = zeta L(chi)/zeta(2s):", ok4)
print("coefficients checked up to", N)
