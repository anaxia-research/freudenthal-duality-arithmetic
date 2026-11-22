"""Closed formulas for the number of BPS orbit classes (positively oriented triples) with integral F-dual."""
from math import gcd
from functools import lru_cache

def chi(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def primes_of(n):
    ps, p = [], 2
    while p*p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps
def psi(f):
    from fractions import Fraction as Fr
    v = Fr(f)
    for p in primes_of(f): v *= (1 - Fr(chi(p), p))
    assert v.denominator == 1
    return int(v)
def r(n): return sum(chi(d) for d in range(1, n+1) if n % d == 0)
def divisors(n): return [d for d in range(1, n+1) if n % d == 0]

def Fplus(k):
    tot = 0; Gp = 0
    from fractions import Fraction as Fr
    for f1 in divisors(k):
        for f2 in divisors(k):
            for f3 in divisors(k):
                if k % (f1*f2*f3): continue
                e = sum(1 for f in (f1, f2, f3) if f > 1)
                v = Fr(psi(f1)*psi(f2)*psi(f3)*r(k//(f1*f2*f3)), 2**min(e, 2))
                tot += v
                if e == 3: Gp += v
    assert tot.denominator == 1 and Gp.denominator == 1
    return int(tot), int(Gp)

def dirichlet_check(k):
    """4F+(k) = A(k) + 3 sigma(k),  A = coefficients of zeta(s) zeta(s-1)^3 / L(s,chi)^2."""
    # coefficients via convolution
    N = k
    one = [0] + [1]*N
    idn = [0] + list(range(1, N+1))
    def mu(n):
        ps = primes_of(n); m = 1
        for p in ps: m *= p
        return 0 if m != n else (-1)**len(ps)
    Linv = [0] + [mu(n)*chi(n) for n in range(1, N+1)]
    def conv(a, b):
        c = [0]*(N+1)
        for i in range(1, N+1):
            if a[i] == 0: continue
            for j in range(1, N//i + 1):
                c[i*j] += a[i]*b[j]
        return c
    A = conv(conv(conv(conv(conv(one, idn), idn), idn), Linv), Linv)
    sig = conv(one, idn)
    return A[k] + 3*sig[k]

if __name__ == "__main__":
    for k in range(1, 41):
        F, G = Fplus(k)
        assert 4*F == dirichlet_check(k)
        print(k, F, G)
