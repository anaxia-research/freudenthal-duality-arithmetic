"""Closed formulas of Theorem E (ST^2 model).
Delta = 4k^2 : F_+(k) = 2 sum_{e1 e2^2 n = k} psi(e1) h(e2) nu(n),  E_+(k) = 2 sum_{e2^2 n = k} h(e2) nu(n),
               h_+(k) = 2 sum_{c d = k} sigma(c) h(d)        (all orbits, both signs of Q_T)
Delta = -m^2 : F_-(m) = sum over (e1, e2, A, B) with A B e1 e2^2 = 4m, A = B mod 2, A, B odd only if e1, e2 even,
               of phi(e1) phi(e2);  E_-(m) = sum_{j^2 e1 e2^2 = m} rho(e1) rho(e2)
               + sum_{2 A^2 e1 e2^2 = m, A, e1, e2 odd} rho(e1) rho(e2);  h_-(m) = m d(m).
psi(e) = e prod_{p | e} (1 - chi(p)/p), h(e) = h(-4e^2) (h(1) = 1, h(e) = psi(e)/2 for e >= 2), nu(n) = sum_{d|n} chi(d),
rho(e) = #{t mod e : t^2 = -1}."""
from math import gcd, isqrt

def chi(n): return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
def divisors(n):
    out = []
    i = 1
    while i*i <= n:
        if n % i == 0:
            out.append(i)
            if i*i != n: out.append(n // i)
        i += 1
    return sorted(out)
def primes_of(n):
    ps = []; p = 2
    while p*p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps
def psi(e):
    from fractions import Fraction
    r = Fraction(e)
    for p in primes_of(e): r *= (1 - Fraction(chi(p), p))
    assert r.denominator == 1
    return int(r)
def hcl(e): return 1 if e == 1 else psi(e) // 2
def nu(n): return sum(chi(d) for d in divisors(n))
def phi(n): return sum(1 for t in range(1, n + 1) if gcd(t, n) == 1)
def ndiv(n): return len(divisors(n))
def sigma(n): return sum(divisors(n))
def rho(e): return sum(1 for t in range(e) if (t*t + 1) % e == 0) if e > 1 else 1

def Fp(k):
    s = 0
    for e1 in divisors(k):
        r = k // e1
        for e2 in divisors(r):
            if r % (e2*e2) == 0: s += psi(e1) * hcl(e2) * nu(r // (e2*e2))
    return 2*s
def Ep(k):
    return 2*sum(hcl(e2) * nu(k // (e2*e2)) for e2 in divisors(k) if k % (e2*e2) == 0)
def hp(k):
    return 2*sum(sigma(c) * hcl(k // c) for c in divisors(k))

def Fm(m):
    s = 0
    for e1 in divisors(4*m):
        r = 4*m // e1
        for e2 in divisors(r):
            if r % (e2*e2): continue
            n = r // (e2*e2)
            if n % 4 == 0: c = ndiv(n // 4)
            elif n % 2 == 1 and e1 % 2 == 0 and e2 % 2 == 0: c = ndiv(n)
            else: c = 0
            s += phi(e1) * phi(e2) * c
    return s
def Em(m):
    s = 0
    for j in divisors(m):
        if m % (j*j): continue
        r = m // (j*j)
        for e1 in divisors(r):
            q = r // e1
            for e2 in divisors(q):
                if q == e2*e2: s += rho(e1) * rho(e2)
    if m % 2 == 0:
        r = m // 2
        for A in divisors(r):
            if A % 2 == 0 or r % (A*A): continue
            q0 = r // (A*A)
            for e1 in divisors(q0):
                q = q0 // e1
                for e2 in divisors(q):
                    if q == e2*e2 and e1 % 2 and e2 % 2: s += rho(e1) * rho(e2)
    return s
def hm(m): return m * ndiv(m)
