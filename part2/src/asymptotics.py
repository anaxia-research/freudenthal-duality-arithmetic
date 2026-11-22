"""Partial sums of F_+, E_+, F_-, E_- up to N against the main terms of Theorem G; writes ../logs/counts.dat."""
import sys, numpy as np, mpmath as mp
N = int(sys.argv[1]); mp.mp.dps = 30
n = np.arange(N + 1)
chi = np.zeros(N + 1, dtype=np.int64); chi[1::4] = 1; chi[3::4] = -1
nu = np.zeros(N + 1, dtype=np.int64); dd = np.zeros(N + 1, dtype=np.int64)
for d in range(1, N + 1):
    if chi[d]: nu[d::d] += chi[d]
    dd[d::d] += 1
def primes_of(e):
    out = []; p = 2
    while p*p <= e:
        if e % p == 0:
            out.append(p)
            while e % p == 0: e //= p
        p += 1
    if e > 1: out.append(e)
    return out
def chi4(p): return 0 if p % 2 == 0 else (1 if p % 4 == 1 else -1)
def psi(e):
    r = e
    for p in primes_of(e): r = r*(p - chi4(p))//p
    return r
def phi(e):
    r = e
    for p in primes_of(e): r = r*(p - 1)//p
    return r
def rho(e):
    if e == 1: return 1
    r = 1
    for p in primes_of(e):
        a = 0; t = e
        while t % p == 0: t //= p; a += 1
        r *= (1 if a == 1 else 0) if p == 2 else 1 + chi4(p)
    return r
Fp = np.zeros(N + 1, dtype=np.int64); Fm = np.zeros(N + 1, dtype=np.int64); Em = np.zeros(N + 1, dtype=np.int64)
e = 1
while e**3 <= N:
    c = e**3; L = N // c
    Fp[c::c] += psi(e) * nu[1:L + 1]; Fm[c::c] += phi(e) * dd[1:L + 1]
    s = 1
    while c*s*s <= N: Em[c*s*s] += rho(e); s += 1
    e += 1
# case (b): f = 2^j f' (f' odd), m = 2^{3j-2} f'^3 n with n odd
j = 1
while 2**(3*j - 2) <= N:
    fp = 1
    while 2**(3*j - 2) * fp**3 <= N:
        if fp % 2:
            c = 2**(3*j - 2) * fp**3; L = N // c
            ns = np.arange(1, L + 1); odd = ns % 2 == 1
            Fm[c*ns[odd]] += phi(2**j * fp) * dd[ns[odd]]
            if j == 1:
                s = 1
                while c*s*s <= N: Em[c*s*s] += rho(2*fp); s += 2
        fp += 1
    j += 1
Ep = nu
SFp, SEp, SFm, SEm = (np.cumsum(x) for x in (Fp, Ep, Fm, Em))
# constants
z2, z3 = mp.zeta(2), mp.zeta(3); L3 = (mp.zeta(3, 0.25) - mp.zeta(3, 0.75))/64
cFp = mp.pi/4 * z2/L3
Z = lambda s: mp.zeta(3*s - 1)/mp.zeta(3*s) * (1 + 2**(-s)*(1 - 2**(-s))**2/(1 - 2**(-3*s)))
c1 = Z(1); c0 = Z(1)*(2*mp.euler - 1) + mp.diff(Z, 1)
L32 = (mp.zeta(1.5, 0.25) - mp.zeta(1.5, 0.75))/8
cE = mp.zeta(1.5)*L32/z3 * (1 + 1/(2*mp.sqrt(2) + 1))
print("c_F+ =", mp.nstr(cFp, 15), " (4/3 =", mp.nstr(mp.mpf(4)/3, 15), ")")
print("c_E+ = pi/4 =", mp.nstr(mp.pi/4, 15), "; BPS share limit 3pi/16 =", mp.nstr(3*mp.pi/16, 15))
print("c1 = (8/7) zeta(2)/zeta(3) =", mp.nstr(c1, 15), mp.nstr(mp.mpf(8)/7*z2/z3, 15))
print("c0 =", mp.nstr(c0, 15))
print("c_E- =", mp.nstr(cE, 15))
print("   X=n        S_F+/main        S_E+/main        S_F-/main        S_E-/main")
rows = []
x = 10
while x <= N:
    r = (float(SFp[x]/(cFp*x)), float(SEp[x]/(mp.pi/4*x)), float(SFm[x]/(c1*x*mp.log(x) + c0*x)), float(SEm[x]/(cE*mp.sqrt(x))))
    print(f"{x:9d}  {r[0]:.9f}  {r[1]:.9f}  {r[2]:.9f}  {r[3]:.9f}")
    rows.append((x, int(SFp[x]), int(SEp[x]), int(SFm[x]), int(SEm[x])))
    x *= 10
print("exact partial sums at", N, ":", int(SFp[N]), int(SEp[N]), int(SFm[N]), int(SEm[N]))
# figure data: shares against log10 X with X = 4K^2 (BPS) and X = M^2 (non-BPS)
with open('../logs/counts.dat', 'w') as fh:
    fh.write("lgXB shareB leadB lgXN scaledN leadN\n")
    for t in np.unique(np.round(np.logspace(1, np.log10(N), 120)).astype(int)):
        XB = 4*t*t; XN = t*t
        sB = SEp[t]/SFp[t]; sN = SEm[t]/SFm[t]
        # scaled non-BPS share: X^{1/4} log X * share
        scaled = sN * XN**0.25 * np.log(XN)
        fh.write(f"{np.log10(XB):.5f} {sB:.6f} {float(3*mp.pi/16):.6f} {np.log10(XN):.5f} {scaled:.6f} {float(2*cE/c1):.6f}\n")
print("wrote ../logs/counts.dat")
