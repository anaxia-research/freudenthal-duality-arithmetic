"""BPS STU charges with Delta = 4k^2 (Bhargava discriminant D = -4k^2): enumerate balanced triples of
positively oriented S-lattices, S = Z[k i], up to equivalence (kappa1 kappa2 kappa3 = 1), and test
(a) integrality of the Freudenthal dual: I1 I2 I3 in k Z[i];
(b) U-equivalence of x and its F-dual: some multiplier ring equals Z[i]."""
from fractions import Fraction as Fr
from math import gcd, isqrt
from lat import Lat, sublattices
from gauss import Field
import sys

F = Field(-1)
I = (Fr(0), Fr(1))

def divisors(n): return [d for d in range(1, n+1) if n % d == 0]

def reduced_forms(D):
    """Reduced primitive positive definite forms of discriminant D < 0."""
    out = []
    A = 1
    while 3*A*A <= -D:
        for B in range(-A+1, A+1):
            if (B*B - D) % (4*A): continue
            C = (B*B - D) // (4*A)
            if C < A: continue
            if B < 0 and (A == C): continue
            if gcd(gcd(A, abs(B)), C) != 1: continue
            out.append((A, B, C))
        A += 1
    return out

def order(f): return Lat([(1, 0), (0, f)])

def reps(k):
    """Lattice classes with multiplier ring S_f = Z[f i], f | k: (J, f, N_S(J))."""
    out = []
    for f in divisors(k):
        for (A, B, C) in reduced_forms(-4*f*f):
            J = Lat([(A, 0), (Fr(-B, 2), f)])
            assert F.mult_ring(J) == order(f)
            NS = J.covol() / order(k).covol()   # Bhargava norm relative to S = S_k
            out.append((J, f, NS))
    return out

def canon(L, allow_i):
    if not allow_i: return L.key()
    return min(L.key(), F.scale(I, L).key())

def triples(k):
    S = order(k); kO = Lat([(k, 0), (0, k)])
    R = reps(k)
    res = []
    for (J1, f1, N1) in R:
        for (J2, f2, N2) in R:
            P12 = F.prod(J1, J2)
            L = F.colon(S, P12)
            n3 = 1 / (N1*N2)
            d = n3 * S.covol() / L.covol()
            if d.denominator != 1: continue
            allow = (f1 == 1 or f2 == 1)
            seen = set()
            for I3 in sublattices(L, int(d)):
                if not I3.contains_lat(F.scale((Fr(0), Fr(k)), I3)): continue   # S-module: k i I3 in I3
                key = canon(I3, allow)
                if key in seen: continue
                seen.add(key)
                O3 = F.mult_ring(I3)
                f3 = int(1 / O3.covol()) if O3.covol() <= 1 else None
                f3 = int(O3.basis()[1][1])
                P = F.prod(P12, I3)
                assert S.contains_lat(P)
                fint = kO.contains_lat(P)
                equiv = (1 in (f1, f2, f3))
                proj = (f1 == f2 == f3 == k)
                res.append(dict(f=(f1, f2, f3), fint=fint, equiv=equiv, proj=proj, J1=J1, J2=J2, I3=I3))
    return res

if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print(" k  h(S)  T+  P+  F+  F+equiv  F+pairs(inequiv)")
    for k in range(1, K+1):
        r = triples(k)
        h = len(reduced_forms(-4*k*k))
        T = len(r); P = sum(x['proj'] for x in r)
        Fi = [x for x in r if x['fint']]
        Fe = sum(x['equiv'] for x in Fi)
        print(f"{k:2d} {h:4d} {T:4d} {P:4d} {len(Fi):4d} {Fe:6d} {len(Fi)-Fe:6d}", "  f-patterns:", sorted(set(x['f'] for x in Fi))[:12])
        assert P == h*h
