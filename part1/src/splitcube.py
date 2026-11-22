"""Cross-check of the non-BPS enumeration with the Bhargava-free canonical form."""
from fractions import Fraction as Fr
from split import triples, dual_triple_class, mul
from cube import Det, T
from canonsq import canon_sq
import sys

def oriented_basis(L):
    (a, b), (_, c) = L.basis()
    return ((a, b), (Fr(0), -c))      # det < 0: positively oriented for our tau

def cube_from(m, J1, J2, I3):
    al, be, ga = oriented_basis(J1), oriented_basis(J2), oriented_basis(I3)
    a = [0]*8
    for A in range(2):
        for B in range(2):
            for C in range(2):
                x = mul(mul(al[A], be[B]), ga[C])
                v = (x[0] - x[1]) / m
                assert v.denominator == 1
                a[4*A + 2*B + C] = int(v)
    return tuple(a)

if __name__ == "__main__":
    M = int(sys.argv[1])
    for m in range(1, M+1):
        r = triples(m)
        keys = set(); nF = 0; nE = 0; nEideal = 0
        for x in r:
            a = cube_from(m, x['J1'], x['J2'], x['I3'])
            assert Det(a) == m*m, (m, Det(a))
            fint = all(t % m == 0 for t in T(a))
            assert fint == x['fint']
            c = canon_sq(a); keys.add(c)
            if fint:
                nF += 1
                ad = tuple(-t // m for t in T(a))
                same = canon_sq(ad) == c
                (f1, t1), (f2, t2), I3p = dual_triple_class(x)
                ideal_same = (f1, t1) == (x['f'][0], x['t'][0]) and (f2, t2) == (x['f'][1], x['t'][1]) and I3p == x['I3']
                assert same == ideal_same, (m, x['f'], same, ideal_same)
                nE += same
        print(f"m={m}: orbits {len(keys)} (triples {len(r)}), F-integral {nF}, F-dual U-equivalent {nE}", flush=True)
        assert len(keys) == len(r)
