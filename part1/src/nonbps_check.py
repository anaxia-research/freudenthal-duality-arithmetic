"""Non-BPS: check that x~ is U-equivalent to x iff (i) every face form Q_j is properly equivalent to -Q_j
and (ii) eps * conj(P) = P for P = I1 I2 I3; and that (i) for a face holds iff t^2 = -1 mod f."""
from split import triples, dual_triple_class, prod, conjL, scale, EPS
from splitcube import cube_from
from cube import form, T, act
from canonsq import canon_sq, isotropic
from math import gcd, isqrt
import sys

def form_class(Q, m):
    """canonical SL2(Z) class of a form of discriminant m^2: min over isotropic lines of (A mod |B|, B, 0)."""
    A, B, C = Q
    best = None
    for v in isotropic(Q, m):
        v1, v2 = v
        # change variables so that v becomes (0,1): Q'(x,y) = Q(x*w + y*v) with (w, v) unimodular
        from canonsq import egcd
        g0, p, q = egcd(v2, v1)
        if g0 < 0: p, q = -p, -q
        w = (p, -q)          # p*v2 - (-q)*v1 = 1  -> det [[w],[v]] = w1 v2 - w2 v1 = 1
        assert w[0]*v2 - w[1]*v1 == 1
        Ev = lambda x, y: A*x*x + B*x*y + C*y*y
        A2 = Ev(*w); C2 = Ev(*v); B2 = Ev(w[0]+v[0], w[1]+v[1]) - A2 - C2
        assert C2 == 0 and abs(B2) == m
        key = (A2 % abs(B2), B2, 0)
        best = key if best is None or key < best else best
    return best

M = int(sys.argv[1])
bad = 0; tot = 0
for m in range(1, M+1):
    for x in triples(m):
        a = cube_from(m, x['J1'], x['J2'], x['I3'])
        cond_i = all(form_class(form(a, j), m) == form_class(tuple(-c for c in form(a, j)), m) for j in range(3))
        # t^2 = -1 mod f for each face: J1, J2 known classes; I3 class via canonical class
        from split import canon_class
        f3, t3, _ = canon_class(x['I3'])
        ts = [(x['f'][0], x['t'][0]), (x['f'][1], x['t'][1]), (f3, t3)]
        cond_t = all((t*t + 1) % f == 0 for f, t in ts)
        assert cond_i == cond_t, (m, ts)
        if not x['fint']: continue
        tot += 1
        P = prod(prod(x['J1'], x['J2']), x['I3'])
        cond_ii = (scale(EPS, conjL(P)) == P)
        (f1, t1), (f2, t2), I3p = dual_triple_class(x)
        same = (f1, t1) == (x['f'][0], x['t'][0]) and (f2, t2) == (x['f'][1], x['t'][1]) and I3p == x['I3']
        if same != (cond_i and cond_ii): bad += 1
print(f"checked all classes with Delta = -m^2, m <= {M}: (i) <=> t^2=-1 mod f always; F-integral classes {tot}; criterion failures {bad}")
