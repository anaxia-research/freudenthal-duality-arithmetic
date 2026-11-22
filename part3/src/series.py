"""Dirichlet series of Theorem F, checked coefficientwise against the closed formulas of Theorem E up to N.
  sum F_+(k) k^-s = zeta(s) zeta(s-1) (1 + zeta(2s-1)/L(2s,chi))
  sum E_+(k) k^-s = zeta(s) L(s,chi) (1 + zeta(2s-1)/L(2s,chi))
  sum h_+(k) k^-s = zeta(s) zeta(s-1) (1 + zeta(s-1)/L(s,chi))
  sum F_-(m) m^-s = zeta(s) zeta(s-1) zeta(2s-1) (1 + 2^{1-s}) / (zeta(2s) (1 + 2^-s))
  sum E_-(m) m^-s = zeta(s) L(s,chi) zeta(2s) L(2s,chi) (1 + 2^-s) / (zeta(4s) (1 + 2^-2s))
  sum h_-(m) m^-s = zeta(s-1)^2
and the auxiliary series sum psi(e) e^-s = zeta(s-1)/L(s,chi), sum phi(e) e^-s = zeta(s-1)/zeta(s), sum rho(e) e^-s = zeta(s) L(s,chi)/zeta(2s)."""
import sys
from counts import Fp, Ep, hp, Fm, Em, hm, psi, phi, rho, chi
N = int(sys.argv[1])
def conv(A, B):
    C = [0]*(N + 1)
    for i in range(1, N + 1):
        if A[i] == 0: continue
        for j in range(1, N // i + 1): C[i*j] += A[i]*B[j]
    return C
def at_squares(fn, power=2):
    out = [0]*(N + 1); e = 1
    while e**power <= N: out[e**power] = fn(e); e += 1
    return out
def mobius(n):
    r = 1; p = 2
    while p*p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0: return 0
            r = -r
        p += 1
    return -r if n > 1 else r
one = [0] + [1]*N; idn = list(range(N + 1)); ch = [0] + [chi(n) for n in range(1, N + 1)]
mu = [0] + [mobius(n) for n in range(1, N + 1)]
muchi = [mu[n]*ch[n] for n in range(N + 1)]
sig = conv(one, idn); nuL = conv(one, ch)
ok_aux = (all(conv(idn, muchi)[n] == psi(n) for n in range(1, N + 1)) and
          all(conv(idn, mu)[n] == phi(n) for n in range(1, N + 1)) and
          all(conv(nuL, at_squares(mobius))[n] == rho(n) for n in range(1, N + 1)))
psq = at_squares(psi)
FP = [x + y for x, y in zip(sig, conv(sig, psq))]
EP = [x + y for x, y in zip(nuL, conv(nuL, psq))]
HP = [x + y for x, y in zip(sig, conv(sig, [0] + [psi(n) for n in range(1, N + 1)]))]
okP = all(FP[n] == Fp(n) and EP[n] == Ep(n) and HP[n] == hp(n) for n in range(1, N + 1))
# (1 + 2x)/(1 + x) = 1 + sum_{j>=1} (-1)^{j-1} x^j and (1 + x)/(1 + x^2), x = 2^-s
c2 = [0]*(N + 1); c2[1] = 1; j = 1
while 2**j <= N: c2[2**j] = (-1)**(j - 1); j += 1
c3 = [0]*(N + 1); j = 0
while 2**j <= N: c3[2**j] = (-1)**(j // 2); j += 1
FM = conv(conv(sig, at_squares(phi)), c2)
nusq = at_squares(lambda j: sum(chi(d) for d in range(1, j + 1) if j % d == 0))
EM = conv(conv(conv(nuL, nusq), at_squares(mobius, 4)), c3)
HM = conv(idn, idn)
okM = all(FM[n] == Fm(n) and EM[n] == Em(n) and HM[n] == hm(n) for n in range(1, N + 1))
print("auxiliary series for psi, phi, rho:", ok_aux)
print("F_+, E_+, h_+ series:", okP)
print("F_-, E_-, h_- series:", okM)
print("coefficients checked up to", N)
