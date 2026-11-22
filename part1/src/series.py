"""Coefficient arrays of the Dirichlet series for F_+, E_+, F_-, E_- (fast sieve), consistency with the
finite formulas, and partial sums against the leading asymptotics."""
import numpy as np, math, sys
from formula import Fplus

def mobius_chi(N):
    mu = np.ones(N+1, dtype=np.int64); mu[0] = 0
    isp = np.ones(N+1, dtype=bool); isp[:2] = False
    for p in range(2, N+1):
        if isp[p]:
            isp[2*p::p] = False
            mu[p::p] *= -1
            mu[p*p::p*p] = 0
    return mu

def conv(a, b, N):
    c = np.zeros(N+1, dtype=object if a.dtype == object else np.int64)
    for i in range(1, N+1):
        if a[i]:
            c[i::i][:N//i] += a[i]*b[1:N//i+1]
    return c

def arrays(N):
    n = np.arange(N+1, dtype=np.int64)
    one = np.ones(N+1, dtype=np.int64); one[0] = 0
    chi = np.array([0 if k % 2 == 0 else (1 if k % 4 == 1 else -1) for k in range(N+1)], dtype=np.int64)
    mu = mobius_chi(N)
    Linv = mu*chi
    idn = n.copy()
    z1 = conv(one, idn, N)                    # zeta(s) zeta(s-1) = sigma
    z13 = conv(conv(z1, idn, N), idn, N)      # zeta(s) zeta(s-1)^3
    Fp = conv(conv(z13, Linv, N), Linv, N) + 3*z1
    z12 = conv(z1, idn, N)
    r = conv(one, chi, N)
    Ep = 3*conv(z12, Linv, N) + r
    # F_- : zeta(s-1)^3 / ((1-2^-s) zeta(s)) ;  1/((1-2^-s)zeta(s)) has coefficients mu(n) * (odd-part indicator?) -> compute
    w = np.zeros(N+1, dtype=np.int64)   # coefficients of prod_{p odd} (1 - p^-s) = mu(n) for odd n
    w[1::2] = mu[1::2]
    Fm = conv(conv(conv(w, idn, N), idn, N), idn, N)
    return Fp, Ep, Fm

if __name__ == "__main__":
    N = int(sys.argv[1])
    Fp, Ep, Fm = arrays(N)
    for k in range(1, 41):
        F, G = Fplus(k)
        assert Fp[k] == 4*F and Ep[k] == 4*(F - G), k
    print("Dirichlet-series coefficients agree with the finite formulas for k <= 40")
    G = 0.915965594177219015
    for K in (10**3, 10**4, N):
        X = 4*K*K
        sF = int(Fp[1:K+1].sum()); sE = int(Ep[1:K+1].sum())
        print(f"BPS  X=4K^2, K={K}: N_F/lead = {sF/(math.pi**2*X*math.log(X)**2/(384*G*G)):.4f},"
              f" N_E/lead = {sE/(math.pi**2*X*math.log(X)/(32*G)):.4f}")
        sFm = int(Fm[1:K+1].sum()); X = K*K
        print(f"nonBPS X=K^2, K={K}: N_F/lead = {sFm/(X*math.log(X)**2/(2*math.pi**2)):.4f}")
