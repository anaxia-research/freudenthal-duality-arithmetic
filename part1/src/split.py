"""Non-BPS STU charges: Delta = -m^2 (Bhargava discriminant D = m^2 > 0). K = Q x Q, S = {(u,v): u = v mod m}.
Enumerate balanced triples of positively oriented S-lattices up to equivalence, test
(a) integrality of the F-dual:  eps P in S, eps = (1,-1), P = I1 I2 I3;
(b) U-equivalence of x and its F-dual: compare the canonical class of (eps Ibar1, eps Ibar2, eps Ibar3) with (I1, I2, I3)."""
from fractions import Fraction as Fr
from math import gcd
from lat import Lat, sublattices
import sys

def mul(u, v): return (u[0]*v[0], u[1]*v[1])
def conj(u): return (u[1], u[0])
def inv(u): return (1/Fr(u[0]), 1/Fr(u[1]))
def scale(k, L): return Lat([mul(k, v) for v in L.basis()])
def prod(L, M): return Lat([mul(u, v) for u in L.basis() for v in M.basis()])
def conjL(L): return Lat([conj(v) for v in L.basis()])
def colon(S, M):
    """(S : M) = {x : x m in S for m in M}, by duality (works with zero divisors)."""
    (a, b), (_, c) = S.basis()
    det = a*c
    Binv = ((c/det, -b/det), (Fr(0), a/det))     # inverse of [[a,b],[0,c]] (row convention z -> z Binv)
    gens = []
    for m in M.basis():
        # C = diag(m) Binv ; columns of C
        C = ((m[0]*Binv[0][0], m[0]*Binv[0][1]), (m[1]*Binv[1][0], m[1]*Binv[1][1]))
        gens.append((C[0][0], C[1][0])); gens.append((C[0][1], C[1][1]))
    return Lat(gens).dual()
def order(f): return Lat([(1, 1), (0, f)])
def conductor(L):
    O = colon(L, L)
    (a, b), (_, c) = O.basis()
    assert a == 1 and b == 1 or (a == 1 and b % c == 1 % c)
    return int(c)
def divisors(n): return [d for d in range(1, n+1) if n % d == 0]
EPS = (Fr(1), Fr(-1))

def It(f, t): return Lat([(1, t), (0, f)]) if f > 1 else Lat([(1, 0), (0, 1)])

def canon_class(L):
    """L invertible over S_f: return (f, t, lam) with N(lam) > 0 and lam L = I_t."""
    f = conductor(L)
    (a, b), (_, c) = L.basis()          # basis rows (a, b), (0, c)
    # p1(L) = a Z ; p2(L) = gcd(b, c) Z
    p2 = Fr(gcd((b*c.denominator*b.denominator).numerator, (c*c.denominator*b.denominator).numerator), c.denominator*b.denominator)
    lam = (1/a, 1/p2)
    Lp = scale(lam, L)
    (a1, b1), (_, c1) = Lp.basis()
    assert a1 == 1 and c1 == f, (L, Lp, f)
    t = int(b1) % f if f > 1 else 0
    return f, t, lam

def reps(m):
    out = []
    for f in divisors(m):
        ts = [t for t in range(f) if gcd(t, f) == 1] if f > 1 else [0]
        for t in ts:
            J = It(f, t)
            assert conductor(J) == f
            out.append((J, f, t))
    return out

def NS(L, m): return L.covol() / order(m).covol()

def triples(m):
    S = order(m)
    R = reps(m)
    res = []
    for (J1, f1, t1) in R:
        for (J2, f2, t2) in R:
            P12 = prod(J1, J2)
            L = colon(S, P12)
            n3 = 1 / (NS(J1, m)*NS(J2, m))
            d = n3 * S.covol() / L.covol()
            if d.denominator != 1: continue
            for I3 in sublattices(L, int(d)):
                tau = (Fr(m % 2 + m, 2), Fr(m % 2 - m, 2))
                if not all(I3.contains(mul(tau, v)) for v in I3.basis()): continue
                f3 = conductor(I3)
                P = prod(P12, I3)
                assert S.contains_lat(P)
                fint = S.contains_lat(scale(EPS, P))
                res.append(dict(f=(f1, f2, f3), t=(t1, t2), fint=fint, J1=J1, J2=J2, I3=I3))
    return res

def dual_triple_class(x):
    """Class of the F-dual triple (eps Ibar1, eps Ibar2, eps Ibar3) in canonical form (J1', J2', I3')."""
    A = [scale(EPS, conjL(L)) for L in (x['J1'], x['J2'], x['I3'])]
    f1, t1, l1 = canon_class(A[0]); f2, t2, l2 = canon_class(A[1])
    I3p = scale(inv(mul(l1, l2)), A[2])
    return (f1, t1), (f2, t2), I3p

if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    print(" m  T+  P+  F+  F+equiv")
    for m in range(1, M+1):
        r = triples(m)
        T = len(r); P = sum(1 for x in r if x['f'] == (m, m, m))
        Fi = [x for x in r if x['fint']]
        Fe = 0
        for x in Fi:
            (f1, t1), (f2, t2), I3p = dual_triple_class(x)
            if (f1, t1) == (x['f'][0], x['t'][0]) and (f2, t2) == (x['f'][1], x['t'][1]) and I3p == x['I3']:
                Fe += 1
        print(f"{m:2d} {T:4d} {P:4d} {len(Fi):4d} {Fe:6d}   patterns {sorted(set(x['f'] for x in Fi))[:10]}", flush=True)
