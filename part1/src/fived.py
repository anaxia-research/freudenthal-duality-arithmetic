"""Five-dimensional STU charges B = (q1,q2,q3), cubic norm N = q1 q2 q3, Jordan dual B* = B^#/N^(1/3):
brute-force check of the integrality and U-duality criteria and of the orbit count j(n)."""
import itertools, math
from math import gcd
def icbrt(N):
    s = 1 if N > 0 else -1; a = round(abs(N)**(1/3))
    for c in (a-1, a, a+1):
        if c**3 == abs(N): return s*c
    return None
def sharp(B): q1, q2, q3 = B; return (q2*q3, q1*q3, q1*q2)
def jdual(B):
    n = icbrt(B[0]*B[1]*B[2])
    if n is None: return None
    S = sharp(B)
    if any(c % n for c in S): return None
    return tuple(c//n for c in S)
G5 = [(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)]
G5hat = [(p, e) for p in itertools.permutations(range(3)) for e in G5]
def act(g, B): p, e = g; return tuple(e[i]*B[p[i]] for i in range(3))
R = 40
bad = 0
for B in itertools.product(range(-R, R+1), repeat=3):
    N = B[0]*B[1]*B[2]
    if N == 0: continue
    n = icbrt(N)
    Bs = jdual(B)
    pred = n is not None and all((B[i]*B[j]) % n == 0 for i, j in ((0,1),(0,2),(1,2)))
    if (Bs is not None) != pred: bad += 1; continue
    if Bs is None: continue
    assert jdual(Bs) == B and Bs[0]*Bs[1]*Bs[2] == N
    u = any(act((tuple(range(3)), e), B) == Bs for e in G5)
    if u != (abs(B[0]) == abs(B[1]) == abs(B[2])) or u != (Bs == B): bad += 1
    uh = any(act(g, B) == Bs for g in G5hat)
    if uh != (abs(n) in map(abs, B)): bad += 1
print("criteria checked for |q_i| <=", R, "failures:", bad)
def j_brute(n):
    c = 0
    for a in range(1, n**3+1):
        if n**3 % a: continue
        for b in range(1, n**3//a+1):
            if (n**3//a) % b: continue
            d = n**3//(a*b)
            if (a*b) % n == 0 and (a*d) % n == 0 and (b*d) % n == 0: c += 1
    return c
def j_formula(n):
    v = 1; m = n; p = 2
    while m > 1:
        if m % p == 0:
            e = 0
            while m % p == 0: m //= p; e += 1
            v *= 3*e*e + 3*e + 1
        p += 1
    return v
for n in range(1, 31):
    assert j_brute(n) == j_formula(n), n
print("j(n) = prod (3e^2+3e+1) verified for n <= 30")
print("examples:", jdual((1, 5, 25)), jdual((2, 2, 16)), jdual((8, 8, 1)), jdual((12, 18, 1)))
# smallest n with an integral J-pair not U-dual even with triality
for n in range(1, 13):
    found = []
    for a in range(1, n**3+1):
        if n**3 % a: continue
        for b in range(1, n**3//a+1):
            if (n**3//a) % b: continue
            d = n**3//(a*b); B = (a, b, d)
            if jdual(B) and n not in B: found.append(B)
    if found: print("smallest n:", n, found); break
# leading constant c = H(1)/720, H(1) = prod_p (1 - 9/p^2 + 16/p^3 - 9/p^4 + 1/p^6)
import sympy
H = 1.0
for p in sympy.primerange(2, 2*10**6):
    x = 1.0/p; H *= 1 - 9*x*x + 16*x**3 - 9*x**4 + x**6
print("H(1) =", H, " c = H(1)/720 =", H/720)
