"""Normal forms of Theorem B, generated from their data and compared with the enumeration (k, m <= K).
Delta = 4k^2, sign eps of Q_T: data (L_S, L_T, alpha) with L_S in X_{e1}, L_T in X_{e2}, alpha in Z[i],
    N(alpha) e1 e2^2 = k;  x + i x~ = alpha xi_S(w1) xi_T(w2) xi_T(w3)  (eps = +1)
                         or alpha xi_S(w1) conj(xi_T(w2)) conj(xi_T(w3))  (eps = -1).
Delta = -m^2: data (e1, t1, e2, t2, A, B), A, B >= 1, A = B mod 2, A, B odd only if e1, e2 even, A B e1 e2^2 = 4m;
    x = (A v1 v2 v3 - B l(w1) n(w2) n(w3)) / 2 with l = e1 u + t1 v, n = e2 u + t2 v."""
import sys
from math import gcd, isqrt
from st2 import *
from counts import Fp, Ep, Fm, Em, psi, hcl, nu

def cmul(z, w): return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])
def cj(z): return (z[0], -z[1])

def X(e):
    """Generators u of the lattices L_u = Z u + e Z[i] in X_e: representatives of (Z[i]/e)^x / (Z/e)^x."""
    if e == 1: return [(1, 0)]
    seen = set(); reps = []
    for x in range(e):
        for y in range(e):
            if gcd((x*x + y*y) % e, e) != 1: continue
            key = min(((s*x) % e, (s*y) % e) for s in range(1, e) if gcd(s, e) == 1)
            if key in seen: continue
            seen.add(key); reps.append((x, y))
    return reps

def hnf_basis(gens):
    """Positively oriented Z-basis (l0, l1) of the lattice in Z^2 = Z[i] spanned by gens."""
    rows = [list(g) for g in gens]
    while True:
        rows = [r for r in rows if r != [0, 0]]
        ys = sorted([r for r in rows if r[1] != 0], key=lambda r: abs(r[1]))
        if len(ys) <= 1: break
        p = ys[0]; new = [p]
        for r in rows:
            if r is p: continue
            if r[1] != 0:
                q = r[1] // p[1]; new.append([r[0] - q*p[0], r[1] - q*p[1]])
            else: new.append(r)
        rows = new
    ys = [r for r in rows if r[1] != 0]; xs = [r for r in rows if r[1] == 0]
    g = 0
    for r in xs: g = gcd(g, r[0])
    l0, l1 = (g, 0), tuple(ys[0])
    if l0[0]*l1[1] - l1[0]*l0[1] < 0: l1 = (-l1[0], -l1[1])
    return l0, l1

def lattice(e, u): return hnf_basis([u, (e, 0), (0, e)])

def bps_charge(alpha, LS, LT, eps):
    """(x, x~) with x + i x~ = alpha xi_S xi_T xi_T (eps = 1) or alpha xi_S conj(xi_T)^2 (eps = -1)."""
    zs = []
    for (A, B, C) in KEYS:
        w2 = LT[B] if eps > 0 else cj(LT[B]); w3 = LT[C] if eps > 0 else cj(LT[C])
        zs.append(cmul(cmul(cmul(alpha, LS[A]), w2), w3))
    return tuple(z[0] for z in zs), tuple(z[1] for z in zs)

def gauss(n):
    r = isqrt(n)
    return [(x, y) for x in range(-r, r + 1) for y in range(-r, r + 1) if x*x + y*y == n]

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]
def unit_lattice(kap, e, u):
    """kappa L_u as an element of X_e (returned as a basis)."""
    v = cmul(kap, u)
    return lattice(e, (v[0] % e, v[1] % e)) if e > 1 else lattice(1, (1, 0))

def gen_bps(k):
    out = []
    for eps in (1, -1):
        for e1 in range(1, k + 1):
            if k % e1: continue
            for e2 in range(1, k + 1):
                if k % (e1*e2*e2): continue
                n = k // (e1*e2*e2)
                for al in gauss(n):
                    for u1 in X(e1):
                        for u2 in X(e2):
                            x, xt = bps_charge(al, lattice(e1, u1), lattice(e2, u2), eps)
                            out.append(((eps, e1, e2, u1, u2, al), x, xt))
    return out

def rank1(A, lam, mu):
    return tuple(A*lam[i]*mu[j]*mu[k] for (i, j, k) in KEYS)

def units_mod(e): return [0] if e == 1 else [t for t in range(e) if gcd(t, e) == 1]
def inv_mod(t, e): return 0 if e == 1 else pow(t, -1, e)

def non_charge(e1, t1, e2, t2, A, B):
    P = rank1(A, (0, 1), (0, 1)); N = rank1(B, (e1, t1), (e2, t2))
    return tuple((p - q)//2 for p, q in zip(P, N)), tuple((p + q)//2 for p, q in zip(P, N))

def gen_non(m):
    out = []
    for e1 in range(1, 4*m + 1):
        if (4*m) % e1: continue
        for e2 in range(1, 4*m + 1):
            if (4*m) % (e1*e2*e2): continue
            n = 4*m // (e1*e2*e2)
            for A in range(1, n + 1):
                if n % A: continue
                B = n // A
                if (A - B) % 2 or (A % 2 and (e1 % 2 or e2 % 2)): continue
                for t1 in units_mod(e1):
                    for t2 in units_mod(e2):
                        x, xt = non_charge(e1, t1, e2, t2, A, B)
                        out.append(((e1, t1, e2, t2, A, B), x, xt))
    return out

if __name__ == "__main__":
    K = int(sys.argv[1])
    for k in range(1, K + 1):
        G = gen_bps(k)
        Ci = {s for s in orbits(4*k*k) if dual(s) is not None}
        fib = {}
        for (eps, e1, e2, u1, u2, al), x, xt in G:
            assert Delta(x) == 4*k*k and dual(x) == xt
            assert (QT(x)[0] > 0) == (eps > 0) and QS(x)[0] > 0
            assert cont(QS(x)) == k // e1 and cont(QT(x)) == k // e2
            c = canon(x); fib.setdefault(c, []).append((eps, e1, e2, u1, u2, al))
            assert (canon(xt) == c) == (e1 == 1)
            # the unit action (kappa1 L_S, kappa2 L_T, kappa1^-1 kappa2^-2 alpha) gives the same orbit
            for k1 in UNITS:
                for k2 in UNITS:
                    k1i = cj(k1); k2i2 = cj(cmul(k2, k2))
                    al2 = cmul(cmul(al, k1i), k2i2)
                    y, yt = bps_charge(al2, unit_lattice(k1, e1, u1), unit_lattice(k2, e2, u2), eps)
                    assert canon(y) == c
            # the dual has the datum (L_S, L_T, -i alpha)
            y, yt = bps_charge(cmul((0, -1), al), lattice(e1, u1), lattice(e2, u2), eps)
            assert y == xt
        assert set(fib) == Ci and len(fib) == Fp(k), (k, len(fib), len(Ci), Fp(k))
        for c, data in fib.items():
            e1s = {d[1] for d in data}; e2s = {d[2] for d in data}
            assert len(e1s) == 1 and len(e2s) == 1
            assert len(data) == (4 if data[0][2] == 1 else 8), (k, c, len(data))
        assert sum(1 for c, d in fib.items() if d[0][1] == 1) == Ep(k)
        Gn = gen_non(k)
        Cn = {s for s in orbits(-k*k) if dual(s) is not None}
        seen = {}
        for (e1, t1, e2, t2, A, B), x, xt in Gn:
            assert Delta(x) == -k*k and dual(x) == xt
            assert cont(QS(x)) == k // e1 and cont(QT(x)) == k // e2
            assert cont([p + q for p, q in zip(x, xt)]) == A and cont([p - q for p, q in zip(x, xt)]) == B
            c = canon(x); assert c not in seen, (k, (e1, t1, e2, t2, A, B), seen.get(c)); seen[c] = 1
            y, yt = non_charge(e1, (-inv_mod(t1, e1)) % e1 if e1 > 1 else 0, e2, (-inv_mod(t2, e2)) % e2 if e2 > 1 else 0, B, A)
            assert canon(y) == canon(xt)
            U = canon(xt) == c
            assert U == (A == B and (t1*t1 + 1) % e1 == 0 and (t2*t2 + 1) % e2 == 0)
        assert set(seen) == Cn and len(seen) == Fm(k)
        print(f"k = m = {k}: {len(fib)} orbits with Delta = 4k^2 and {len(seen)} with Delta = -m^2 from the data; "
              "pairwise inequivalent, exhaustive, contents, dual data and U-duality as stated", flush=True)
    print("generation check passed")
