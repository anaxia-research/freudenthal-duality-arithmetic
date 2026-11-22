"""Cubes from balanced triples (Bhargava's construction) and cross-check with the canonical form."""
from fractions import Fraction as Fr
from bps import triples, reduced_forms, F
from cube import Det, T, act, form
from canon import canon
from formula import Fplus
import sys

def cube_from(k, J1, J2, I3):
    al, be, ga = J1.basis(), J2.basis(), I3.basis()
    a = [0]*8
    for A in range(2):
        for B in range(2):
            for C in range(2):
                x = F.mul(F.mul(al[A], be[B]), ga[C])     # x = u + v i, tau = k i
                v = x[1] / k
                assert v.denominator == 1 and x[0].denominator == 1
                a[4*A + 2*B + C] = int(v)
    return tuple(a)

J = ((0, 1), (1, 0))
if __name__ == "__main__":
    K = int(sys.argv[1])
    for k in range(1, K+1):
        r = triples(k)
        keys = set(); nF = 0; nE = 0
        for x in r:
            a = cube_from(k, x['J1'], x['J2'], x['I3'])
            assert Det(a) == -4*k*k
            m = 2*k
            fint = all(t % m == 0 for t in T(a))
            assert fint == x['fint']
            variants = [a, act(act(a, 0, J), 1, J), act(act(a, 0, J), 2, J), act(act(a, 1, J), 2, J)]
            for b in variants:
                c = canon(b); keys.add(c)
                if fint:
                    nF += 1
                    bd = tuple(t // m for t in T(b))
                    same = canon(bd) == c
                    assert same == x['equiv'], (k, x['f'])
                    nE += same
        Fp, Gp = Fplus(k)
        print(f"k={k}: orbits {len(keys)} (=4*{len(r)}), F-integral {nF} (formula {4*Fp}), F-dual U-equivalent {nE} (formula {4*(Fp-Gp)})", flush=True)
        assert len(keys) == 4*len(r) and nF == 4*Fp and nE == 4*(Fp - Gp)
