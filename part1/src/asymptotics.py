"""Partial sums of F_+, E_+, F_-, E_- against (a) the leading constants and (b) the full main terms
Res_{s=s0} F(s) X^s / s computed by contour integration; proportions of U-dual pairs; table of counts."""
import numpy as np, math, sys, mpmath as mp
from series import arrays, mobius_chi, conv

mp.mp.dps = 30
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10**5

def rho_array(N):
    """rho(n) = #{t mod n : t^2 = -1}: multiplicative, rho(2)=1, rho(4^..)=0, rho(p^e)=1+chi(p) for odd p."""
    r = np.ones(N+1, dtype=np.int64); r[0] = 0
    spf = np.zeros(N+1, dtype=np.int64)
    for p in range(2, N+1):
        if spf[p] == 0:
            spf[p::p][spf[p::p] == 0] = p
    for n in range(2, N+1):
        p = spf[n]; m = n; e = 0
        while m % p == 0: m //= p; e += 1
        if p == 2: lp = 1 if e == 1 else 0
        else: lp = 2 if p % 4 == 1 else 0
        r[n] = r[m]*lp
    return r

def Eminus_array(N):
    rho = rho_array(N)
    r3 = conv(conv(rho, rho, N), rho, N)
    sq = np.zeros(N+1, dtype=np.int64)
    a = 1
    while a*a <= N: sq[a*a] = 1; a += 1
    Tm = conv(r3, sq, N)
    E = np.zeros(N+1, dtype=np.int64)
    for m in range(1, N+1):
        mm = m
        while mm % 2 == 0: mm //= 2
        E[m] = (4 if m % 2 == 0 else 1)*Tm[mm]
    return E

def L4(s): return mp.dirichlet(s, [0, 1, 0, -1])
def Fp_s(s): return mp.zeta(s)*mp.zeta(s-1)**3/L4(s)**2 + 3*mp.zeta(s)*mp.zeta(s-1)
def Ep_s(s): return 3*mp.zeta(s)*mp.zeta(s-1)**2/L4(s) + mp.zeta(s)*L4(s)
def Fm_s(s): return mp.zeta(s-1)**3/((1-mp.power(2, -s))*mp.zeta(s))
def Em_s(s): return (1+3*mp.power(2, -s))*mp.zeta(s)**3*L4(s)**3/((1+mp.power(2, -s))**2*mp.zeta(2*s)**2)

def main_term(F, s0, X, r=mp.mpf('0.35')):
    f = lambda t: F(s0 + r*mp.expj(t))*mp.power(X, s0 + r*mp.expj(t))/(s0 + r*mp.expj(t))*r*mp.expj(t)
    return mp.re(mp.quad(f, [0, mp.pi/2, mp.pi, 3*mp.pi/2, 2*mp.pi])/(2*mp.pi))

if __name__ == "__main__":
    Fp, Ep, Fm = arrays(N)
    Em = Eminus_array(N)
    # check E_- against the Dirichlet series (1+3*2^-s)/(1+2^-s)^2 * zeta^3 L^3 / zeta(2s)^2 for small n
    n = 2000
    one = np.zeros(n+1, dtype=np.int64); one[1:] = 1
    chi = np.array([0 if k % 2 == 0 else (1 if k % 4 == 1 else -1) for k in range(n+1)], dtype=np.int64); chi[0] = 0
    z3 = conv(conv(one, one, n), one, n); L3 = conv(conv(chi, chi, n), chi, n)
    mu = mobius_chi(n); mm = conv(mu, mu, n)
    inv2 = np.zeros(n+1, dtype=np.int64)
    a = 1
    while a*a <= n: inv2[a*a] = mm[a]; a += 1
    # (1+3x)/(1+x)^2 = sum c_j x^j
    c = [1]; j = 1
    while 2**j <= n:
        # coefficients of (1+3x) * sum_{i>=0} (i+1)(-x)^i
        c.append(3*((-1)**(j-1))*j + ((-1)**j)*(j+1)); j += 1
    two = np.zeros(n+1, dtype=np.int64)
    for j, cj in enumerate(c): two[2**j] = cj
    ser = conv(conv(conv(z3, L3, n), inv2, n), two, n)
    assert all(ser[k] == Em[k] for k in range(1, n+1)), "E_- Dirichlet series mismatch"
    print(f"E_- closed formula agrees with its Dirichlet series for m <= {n}")
    G = float(mp.catalan)
    rows = []
    print(" X              N_F+/lead  N_F+/main  N_E+/lead  N_E+/main | N_F-/lead  N_F-/main  N_E-/lead  N_E-/main | U-dual share BPS, non-BPS")
    Ks = [10**j for j in range(2, int(round(math.log10(N)))+1)]
    for K in Ks:
        sFp = int(Fp[1:K+1].sum()); sEp = int(Ep[1:K+1].sum())
        sFm = int(Fm[1:K+1].sum()); sEm = int(Em[1:K+1].sum())
        Xb = 4*K*K; Xn = K*K
        lFp = math.pi**2*Xb*math.log(Xb)**2/(384*G*G); lEp = math.pi**2*Xb*math.log(Xb)/(32*G)
        lFm = Xn*math.log(Xn)**2/(2*math.pi**2); lEm = 5/(64*math.pi)*math.sqrt(Xn)*math.log(Xn)**2
        mFp = float(main_term(Fp_s, 2, K)); mEp = float(main_term(Ep_s, 2, K))
        mFm = float(main_term(Fm_s, 2, K)); mEm = float(main_term(Em_s, 1, K))
        print(f" 10^{math.log10(K):.0f} (K) {sFp/lFp:10.4f} {sFp/mFp:10.5f} {sEp/lEp:10.4f} {sEp/mEp:10.5f} |"
              f" {sFm/lFm:10.4f} {sFm/mFm:10.5f} {sEm/lEm:10.4f} {sEm/mEm:10.5f} | {sEp/sFp:.4f} {sEm/sFm:.6f}", flush=True)
    with open("../logs/counts.dat", "w") as f:
        f.write("k Fp Ep Fm Em\n")
        for k in range(1, 61):
            f.write(f"{k} {Fp[k]//4} {Ep[k]//4} {Fm[k]} {Em[k]}\n")
    print("wrote counts.dat")
