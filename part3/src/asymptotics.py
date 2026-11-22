"""Partial sums of F_+, E_+, h_+ (k <= K) and F_-, E_-, h_- (m <= K) up to N, against the main terms of Theorem G.
For E_+, h_+ (double pole) and E_- (pole of zeta(2s) at s = 1/2) the next term is included; for F_+, F_- and h_-
the remainder has the size of the remainder in the sums of sigma(n) and d(n) and is shown normalized.
Writes ../logs/counts.dat."""
import sys, numpy as np, mpmath as mp
N = int(sys.argv[1]); mp.mp.dps = 30
chi = np.zeros(N + 1, dtype=np.int64); chi[1::4] = 1; chi[3::4] = -1
sig = np.zeros(N + 1, dtype=np.int64); nu = np.zeros(N + 1, dtype=np.int64); dd = np.zeros(N + 1, dtype=np.int64)
for d in range(1, N + 1):
    sig[d::d] += d; dd[d::d] += 1
    if chi[d]: nu[d::d] += chi[d]
# psi(n) = n prod (1 - chi(p)/p), phi(n), mu(n) by a prime sieve
isp = np.ones(N + 1, dtype=bool); isp[:2] = False
for p in range(2, int(N**0.5) + 1):
    if isp[p]: isp[p*p::p] = False
primes = np.nonzero(isp)[0]
psi = np.arange(N + 1, dtype=np.int64); phi = np.arange(N + 1, dtype=np.int64)
for p in primes:
    p = int(p); c = 0 if p == 2 else (1 if p % 4 == 1 else -1)
    psi[p::p] = psi[p::p] // p * (p - c); phi[p::p] = phi[p::p] // p * (p - 1)
r = int(N**0.5) + 2
mu = np.ones(r + 1, dtype=np.int64); mu[0] = 0
for p in range(2, r + 1):
    if all(p % q for q in range(2, int(p**0.5) + 1)):
        mu[p::p] *= -1; mu[p*p::p*p] = 0
# BPS
Fp = sig.copy(); Ep = nu.copy()
e = 1
while e*e <= N:
    c = e*e; L = N // c
    Fp[c::c] += psi[e] * sig[1:L + 1]; Ep[c::c] += psi[e] * nu[1:L + 1]
    e += 1
Fp *= 1; hp = sig.copy()
for d in range(1, N + 1):
    hp[d::d] += psi[d] * sig[1:N // d + 1]
# non-BPS
base = np.zeros(N + 1, dtype=np.int64); e = 1
while e*e <= N:
    c = e*e; base[c::c] += phi[e] * sig[1:N // c + 1]; e += 1
Fm = base.copy(); j = 1
while 2**j <= N:
    q = 2**j; Fm[q::q] += (-1)**(j - 1) * base[1:N // q + 1]; j += 1
w = np.zeros(r + 1, dtype=np.int64)          # w(j) = sum_{i^2 | j} mu(i) nu(j / i^2): coefficient of (nu_sq * mu_4) at j^2
for jj in range(1, r + 1):
    i = 1
    while i*i <= jj:
        if jj % (i*i) == 0: w[jj] += mu[i] * nu[jj // (i*i)]
        i += 1
E0 = np.zeros(N + 1, dtype=np.int64); jj = 1
while jj*jj <= N:
    c = jj*jj
    if w[jj]: E0[c::c] += w[jj] * nu[1:N // c + 1]
    jj += 1
Em = E0.copy(); j = 1
while 2**j <= N:
    q = 2**j; Em[q::q] += (-1)**(j // 2) * E0[1:N // q + 1]; j += 1
hm = np.arange(N + 1, dtype=np.int64) * dd
for arr in (Fp, Ep, hp, Fm, Em, hm): arr[0] = 0
S = {k: np.cumsum(v) for k, v in (("F+", Fp), ("E+", Ep), ("h+", hp), ("F-", Fm), ("E-", Em), ("h-", hm))}
# spot checks against the closed formulas
from counts import Fp as cFp, Ep as cEp, hp as chp, Fm as cFm, Em as cEm, hm as chm
for n in list(range(1, 300)) + [997, 1024, 3000, 7200, 9999]:
    if n > N: continue
    assert (Fp[n], Ep[n], hp[n], Fm[n], Em[n], hm[n]) == (cFp(n), cEp(n), chp(n), cFm(n), cEm(n), chm(n)), n
# constants
Lchi = lambda s: mp.dirichlet(s, [0, 1, 0, -1])
G = mp.catalan; g = mp.euler; z = mp.zeta
C2 = lambda s: (1 + 2**(1 - s))/(1 + 2**(-s)); C3 = lambda s: (1 + 2**(-s))/(1 + 2**(-2*s))
cF = z(2)*(1 + z(3)/Lchi(4))/2
B = lambda s: Lchi(s)/(Lchi(2*s)*s); B0 = B(1); B1 = mp.diff(B, 1)
dE = (B0/2, B1/2 + 3*g*B0/2 + Lchi(1))
H = lambda s: z(s)/(Lchi(s)*s); H0 = H(2); H1 = mp.diff(H, 2)
dh = (H0, H1 + 2*g*H0 + z(2)/2)                               # K^2 (dh0 log K + dh1)
cFm_ = z(2)*z(3)*C2(2)/(2*z(4))
cEm_ = C3(1)*Lchi(1)*z(2)*Lchi(2)/z(4)
sEm = C3(mp.mpf(1)/2)*z(mp.mpf(1)/2)*Lchi(mp.mpf(1)/2)*Lchi(1)/z(2)
print("F_+ : K^2 coefficient zeta(2)(1 + zeta(3)/L(4,chi))/2 =", mp.nstr(cF, 15))
print("E_+ : K (%s log K + %s), leading pi/(8G) =" % (mp.nstr(dE[0], 15), mp.nstr(dE[1], 12)), mp.nstr(mp.pi/(8*G), 15))
print("h_+ : K^2 (%s log K + %s), leading zeta(2)/(2G) =" % (mp.nstr(dh[0], 15), mp.nstr(dh[1], 12)), mp.nstr(z(2)/(2*G), 15))
print("F_- : M^2 coefficient 9 zeta(3)/pi^2 =", mp.nstr(cFm_, 15), mp.nstr(9*z(3)/mp.pi**2, 15))
print("E_- : M coefficient 9G/(2 pi) =", mp.nstr(cEm_, 15), mp.nstr(9*G/(2*mp.pi), 15), "; second term", mp.nstr(sEm, 12), "M^(1/2)")
print("h_- : (M^2/2)(log M + 2 gamma - 1/2)")
print("shares (X = |Delta|): U-dual BPS pairs ~ s_+ log X / sqrt X with s_+ = pi/(4 G zeta(2) (1 + zeta(3)/L(4,chi))) =",
      mp.nstr(mp.pi/(4*G*z(2)*(1 + z(3)/Lchi(4))), 15))
print("  U-dual non-BPS pairs ~ s_- / sqrt X with s_- = pi G/(2 zeta(3)) =", mp.nstr(mp.pi*G/(2*z(3)), 15))
print("  integral duals among all BPS orbits ~ 2G(1 + zeta(3)/L(4,chi))/log X =", mp.nstr(2*G*(1 + z(3)/Lchi(4)), 15), "/ log X")
print("  integral duals among all non-BPS orbits ~ 36 zeta(3)/(pi^2 log X) =", mp.nstr(36*z(3)/mp.pi**2, 15), "/ log X")
mainF = lambda K: cF*K**2
mainE = lambda K: K*(dE[0]*mp.log(K) + dE[1])
mainh = lambda K: K**2*(dh[0]*mp.log(K) + dh[1])
mainFm = lambda M: cFm_*M**2
mainEm = lambda M: cEm_*M + sEm*mp.sqrt(M)
mainhm = lambda M: M**2/2*(mp.log(M) + 2*g - mp.mpf(1)/2)
print("        K    S_F+/main  S_E+/main  S_h+/main  S_F-/main  S_E-/main  S_h-/main")
x = 10
while x <= N:
    rr = [float(S["F+"][x]/mainF(x)), float(S["E+"][x]/mainE(x)), float(S["h+"][x]/mainh(x)),
          float(S["F-"][x]/mainFm(x)), float(S["E-"][x]/mainEm(x)), float(S["h-"][x]/mainhm(x))]
    print(f"{x:9d}  " + "  ".join(f"{t:.7f}" for t in rr))
    x *= 10
print("        K  (S_F+ - main)/(K log K)  (S_E+ - main)/sqrt(K)  (S_h+ - main)/K^(3/2)  (S_F- - main)/(M log M)  (S_E- - main)/M^(1/4)  (S_h- - main)/M^(3/2)")
x = 10
while x <= N:
    lx = mp.log(x)
    rr = [(S["F+"][x] - mainF(x))/(x*lx), (S["E+"][x] - mainE(x))/mp.sqrt(x), (S["h+"][x] - mainh(x))/x**1.5,
          (S["F-"][x] - mainFm(x))/(x*lx), (S["E-"][x] - mainEm(x))/x**0.25, (S["h-"][x] - mainhm(x))/x**1.5]
    print(f"{x:9d}  " + "  ".join(f"{float(t):+.5f}" for t in rr))
    x *= 10
print("exact partial sums at K = M =", N, ":", ", ".join(f"{k} {int(S[k][N])}" for k in ("F+", "E+", "h+", "F-", "E-", "h-")))
print("main terms at", N, ":", ", ".join(f"{k} {mp.nstr(f(N), 12)}" for k, f in (("F+", mainF), ("E+", mainE), ("h+", mainh), ("F-", mainFm), ("E-", mainEm), ("h-", mainhm))))
sp_ = mp.pi/(4*G*z(2)*(1 + z(3)/Lchi(4))); sm_ = mp.pi*G/(2*z(3))
print("STU-U-dual share of the integral ST^2 BPS pairs: 2/(1 + zeta(3)/L(4,chi)) =", mp.nstr(2/(1 + z(3)/Lchi(4)), 15))
with open('../logs/counts.dat', 'w') as fh:
    # Figure 3: BPS: share * sqrt(X) / log X against log10 X (X = 4K^2), limit s_+;
    #           non-BPS: share * sqrt(X) against log10 X (X = M^2), limit s_-; pred = the same for the main terms
    fh.write("lgXB yB pB lgXN yN pN\n")
    for t in np.unique(np.round(np.logspace(1, np.log10(N), 100)).astype(int)):
        XB = 4.0*t*t; XN = 1.0*t*t
        sB = S["E+"][t]/S["F+"][t]; sN = S["E-"][t]/S["F-"][t]
        pB = float(mainE(t)/mainF(t)); pN = float(mainEm(t)/mainFm(t))
        fh.write(f"{np.log10(XB):.5f} {sB*np.sqrt(XB)/np.log(XB):.6f} {pB*np.sqrt(XB)/np.log(XB):.6f} {np.log10(XN):.5f} {sN*np.sqrt(XN):.6f} {pN*np.sqrt(XN):.6f}\n")
print("limits: s_+ =", mp.nstr(sp_, 10), " s_- =", mp.nstr(sm_, 10))
print("wrote ../logs/counts.dat")
