"""Checks of statements used in the text of Part III that are not covered by the other scripts:
the explicit triple product, covariance of R_S, the involution iota, the two real types, the T^3 embedding,
the rank-one factorizations, stabilizers, real dualities and Freudenthal rotations (30 digits),
the elements of SL2(Z)^2 taking x to its dual, and the worked examples."""
import random
import sympy as sp
from mpmath import mp, mpf, mpc, sqrt, cos, sin, cosh, sinh, matrix, inverse, fabs, exp, pi
from st2 import *
from cube import act as cact, Det as cDet
from stu_canon import canon as stu_canon_def
from stu_canonsq import canon_sq as stu_canon_sq
from generate import bps_charge, non_charge, lattice
mp.dps = 30
random.seed(11)

# ---------- 1. the triple product written out ----------
a, b, c, e, f, g = sp.symbols('a b c e f g')
s = (a, b, c, e, f, g)
t = [sp.expand(z) for z in T(s)]
expl = [-a**2*g + 2*a*b*f + a*c*e - 2*b**2*e, -a*b*g + 2*a*c*f - b*c*e, a*c*g - 2*b**2*g + 2*b*c*f - c**2*e,
        -a*e*g + 2*a*f**2 - 2*b*e*f + c*e**2, a*f*g - 2*b*e*g + c*e*f, a*g**2 - 2*b*f*g - c*e*g + 2*c*f**2]
print("triple product T(x) as written in Section 2:", all(sp.expand(p - q) == 0 for p, q in zip(t, expl)))
D = sp.expand(Delta(s))
print("Delta as written:", sp.expand(D - (-a**2*g**2 + 4*a*b*f*g + 2*a*c*e*g - 4*a*c*f**2 - 4*b**2*e*g + 4*b*c*e*f - c**2*e**2)) == 0)
# Pluecker description: A x B for the rows A = (a, b, c), B = (e, f, g)
A = sp.Matrix([a, b, c]); B = sp.Matrix([e, f, g]); P, q, R = QT(s)
print("A x B = (-R, q, -P) for Q_T = (P, q, R):", list(A.cross(B)) == [sp.expand(-R), sp.expand(q), sp.expand(-P)] or
      sp.simplify(A.cross(B) - sp.Matrix([-R, q, -P])) == sp.zeros(3, 1))
# Q_S is the form b^2 - ac restricted to the plane: Q_S(u, v) = beta(uA + vB) with beta(a, b, c) = b^2 - ac
u, v = sp.symbols('u v')
W = u*A + v*B
print("Q_S(u, v) = (b^2 - ac)(uA + vB):", sp.expand(W[1]**2 - W[0]*W[2] - (QS(s)[0]*u**2 + QS(s)[1]*u*v + QS(s)[2]*v**2)) == 0)

# ---------- 2. covariance of R_S, the involution iota ----------
E = ((1, 0), (0, -1))
def iota(s): return (s[0], -s[1], s[2], s[3], -s[4], s[5])
gens = [((1, 1), (0, 1)), ((0, -1), (1, 0)), ((2, 1), (1, 1)), ((1, 0), (3, 1)), MI]
okR = okI = okTypes = okT3 = True; nb = 0
for _ in range(3000):
    s0 = tuple(random.randint(-6, 6) for _ in range(6))
    g1 = I2; g2 = I2
    for _ in range(3): g1 = mul(g1, random.choice(gens)); g2 = mul(g2, random.choice(gens))
    y = act(g1, g2, s0)
    # R_S(g.x) = g1 R_S(x) g1^{-1}  (det g1 = det g2 = 1)
    g1i = ((g1[1][1], -g1[0][1]), (-g1[1][0], g1[0][0]))
    okR &= RS(y) == mul(mul(g1, RS(s0)), g1i)
    # iota is the action of (1, E), commutes with T and duality, and normalizes the group
    okI &= iota(s0) == uncube(cact(cact(cube(s0), 1, E), 2, E))
    okI &= T(iota(s0)) == iota(T(s0)) and Delta(iota(s0)) == Delta(s0) and QS(iota(s0)) == QS(s0)
    okI &= QT(iota(s0)) == tuple(-z for z in qact(E, QT(s0)))
    okI &= iota(y) == act(g1, mul(mul(E, g2), E), iota(s0))
    d = dual(s0); okI &= (d is None and dual(iota(s0)) is None) or (d is not None and dual(iota(s0)) == iota(d))
    Dl = Delta(s0)
    if Dl > 0:
        nb += 1
        PS, qS, RS_ = QS(s0); PT, qT, RT_ = QT(s0)
        okTypes &= PS > 0 and RS_ > 0 and (PT > 0) == (RT_ > 0) and PT != 0
        okTypes &= (QT(iota(s0))[0] > 0) != (PT > 0)
for _ in range(3000):
    f4 = tuple(random.randint(-6, 6) for _ in range(4))
    a4, b4, c4, d4 = f4
    x = (a4, b4, c4, b4, c4, d4)
    Df = -(a4*a4*d4*d4 - 6*a4*b4*c4*d4 + 4*a4*c4**3 + 4*b4**3*d4 - 3*b4*b4*c4*c4)
    okT3 &= Delta(x) == Df and QS(x) == QT(x) == (b4*b4 - a4*c4, b4*c4 - a4*d4, c4*c4 - b4*d4)
    if Df > 0: okT3 &= QT(x)[0] > 0
print("R_S(g.x) = g1 R_S(x) g1^{-1} for (g1, g2) in SL2(Z)^2 (3000 random charges):", okR)
print("iota = (1, E): commutes with T and with duality, fixes Q_S and Delta, Q_T(iota x)(w) = -Q_T(x)(wE),"
      " iota (g1, g2) iota = (g1, E g2 E):", okI)
print("Delta > 0: Q_S positive definite, Q_T definite, iota changes the sign of Q_T (%d charges):" % nb, okTypes)
print("T^3 charges (a, b, c, d) -> (a, b, c, b, c, d): same Delta, Q_S = Q_T = Hessian, Q_T > 0 when Delta > 0:", okT3)

# ---------- 3. factorizations, stabilizers, real dualities and rotations, 30 digits ----------
def num_cube(s): return [mpf(z) for z in cube(s)]
def c_act(x, slot, gm): return list(cact(tuple(x), slot, gm))
def maxdiff(x, y): return max(fabs(p - q) for p, q in zip(x, y))
def factor_rank1(Phi):
    """Phi_{ABC} = lam_A nu_B nu_C (symmetric in B, C)."""
    BC = max(((B, C) for B in range(2) for C in range(2)), key=lambda p: abs(Phi[2*p[0] + p[1]]) + abs(Phi[4 + 2*p[0] + p[1]]))
    lam = [Phi[4*A + 2*BC[0] + BC[1]] for A in range(2)]
    A0 = max(range(2), key=lambda A: sum(abs(Phi[4*A + j]) for j in range(4)))
    C0 = max(range(2), key=lambda C: abs(Phi[4*A0 + C]) + abs(Phi[4*A0 + 2 + C]))
    nu = [Phi[4*A0 + 2*B + C0] for B in range(2)]
    kap = Phi[4*A0 + 2*BC[0] + BC[1]] / (lam[A0]*nu[BC[0]]*nu[BC[1]])
    lam = [kap*z for z in lam]
    err = max(abs(Phi[4*A + 2*B + C] - lam[A]*nu[B]*nu[C]) for A in range(2) for B in range(2) for C in range(2))
    return lam, nu, err
def rho(vec, z):
    """Real matrix M with w M (vec) = z w (vec) and w M (conj vec) = conj(z) w (conj vec): lam(wM) = z lam(w)."""
    Pm = matrix([[vec[0], vec[0].conjugate()], [vec[1], vec[1].conjugate()]])
    M = Pm * matrix([[z, 0], [0, z.conjugate()]]) * inverse(Pm)
    return ((M[0, 0].real, M[0, 1].real), (M[1, 0].real, M[1, 1].real)), max(abs(M[i, j].imag) for i in range(2) for j in range(2))
def realdiag(vecs, zs):
    Pm = matrix([[vecs[0][0], vecs[1][0]], [vecs[0][1], vecs[1][1]]])
    M = Pm * matrix([[zs[0], 0], [0, zs[1]]]) * inverse(Pm)
    return ((M[0, 0], M[0, 1]), (M[1, 0], M[1, 1]))
def m2(M): return ((mpf(M[0][0]), mpf(M[0][1])), (mpf(M[1][0]), mpf(M[1][1])))
def mdet(M): return M[0][0]*M[1][1] - M[0][1]*M[1][0]
def mmul(X, Y): return tuple(tuple(sum(X[i][k]*Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))
worst = mpf(0); nB = nN = 0
while nB < 200 or nN < 200:
    s0 = tuple(random.randint(-9, 9) for _ in range(6)); Dl = Delta(s0)
    if Dl == 0: continue
    if Dl > 0 and nB >= 200: continue
    if Dl < 0 and nN >= 200: continue
    sg = 1 if Dl > 0 else -1; rt = sqrt(mpf(abs(Dl)))
    x = num_cube(s0); xt = [sg*z/rt for z in num_cube(T(s0))]
    DS = tuple(tuple(sg*mpf(z)/rt for z in row) for row in RS(s0))
    DT = tuple(tuple(sg*mpf(z)/rt for z in row) for row in RT(s0))
    errs = [maxdiff(c_act(x, 0, DS), xt), fabs(mdet(DS) - sg)]
    errs.append(maxdiff(c_act(c_act(x, 1, DT), 2, DT), [-sg*z for z in x]))
    th = mpf(random.random()*6)
    if Dl > 0:
        nB += 1
        DS2 = mmul(DS, DS); errs.append(max(fabs(DS2[0][0] + 1), fabs(DS2[1][1] + 1), fabs(DS2[0][1]), fabs(DS2[1][0])))
        rot = ((cos(th) + sin(th)*DS[0][0], sin(th)*DS[0][1]), (sin(th)*DS[1][0], cos(th) + sin(th)*DS[1][1]))
        errs.append(maxdiff(c_act(x, 0, rot), [cos(th)*p + sin(th)*q for p, q in zip(x, xt)]))
        Phi = [mpc(p, q) for p, q in zip(x, xt)]
        lam, nu, err = factor_rank1(Phi); errs.append(err)
        # lam(w D_S) = -i lam(w), so D_S = rho_S(-i)
        errs.append(max(abs(DS[A][0]*lam[0] + DS[A][1]*lam[1] + 1j*lam[A]) for A in range(2)))
        # Q_S is a positive multiple of |lam|^2 and Q_T a multiple of |nu|^2 with the sign of Q_T
        QSs = QS(s0); QTs = QT(s0)
        rS = [QSs[0]/abs(lam[0])**2 if lam[0] != 0 else None]
        for w in [(1, 0), (0, 1), (1, 1), (1, -2)]:
            lw = lam[0]*w[0] + lam[1]*w[1]; nw = nu[0]*w[0] + nu[1]*w[1]
            qs = QSs[0]*w[0]**2 + QSs[1]*w[0]*w[1] + QSs[2]*w[1]**2; qt = QTs[0]*w[0]**2 + QTs[1]*w[0]*w[1] + QTs[2]*w[1]**2
            rS.append(qs/abs(lw)**2)
            rT = qt/abs(nw)**2
            errs.append(0 if (rT > 0) == (QTs[0] > 0) else 1)
        errs.append(max(fabs(z - rS[1]) for z in rS[1:]) / fabs(rS[1]))
        # stabilizer (rho_S(s^-2), rho_T(s)) and real dualities (rho_S(-i s^-2), rho_T(s))
        for _ in range(3):
            z = exp(1j*mpf(random.random()*6))
            g1, i1 = rho(lam, z**-2); g2, i2 = rho(nu, z); errs += [i1, i2]
            y = c_act(c_act(c_act(x, 0, g1), 1, g2), 2, g2)
            errs.append(maxdiff(y, x))
            g1d, i3 = rho(lam, -1j*z**-2); errs.append(i3)
            y = c_act(c_act(c_act(x, 0, g1d), 1, g2), 2, g2)
            errs.append(maxdiff(y, xt))
    else:
        nN += 1
        DS2 = mmul(DS, DS); errs.append(max(fabs(DS2[0][0] - 1), fabs(DS2[1][1] - 1), fabs(DS2[0][1]), fabs(DS2[1][0])))
        rot = ((cosh(th) + sinh(th)*DS[0][0], sinh(th)*DS[0][1]), (sinh(th)*DS[1][0], cosh(th) + sinh(th)*DS[1][1]))
        errs.append(maxdiff(c_act(x, 0, rot), [cosh(th)*p + sinh(th)*q for p, q in zip(x, xt)]))
        l1, m1, e1 = factor_rank1([p + q for p, q in zip(x, xt)])
        l2, m2_, e2 = factor_rank1([p - q for p, q in zip(x, xt)])
        errs += [e1, e2]
        # torus stabilizer: lam1 -> s^-2 lam1, lam2 -> s^2 lam2, mu1 -> s mu1, mu2 -> s^-1 mu2
        for _ in range(3):
            z = mpf(random.random()*3 + 0.2)
            g1 = realdiag([l1, l2], [z**-2, z**2]); g2 = realdiag([m1, m2_], [z, 1/z])
            errs += [fabs(mdet(g1) - 1), fabs(mdet(g2) - 1)]
            errs.append(maxdiff(c_act(c_act(c_act(x, 0, g1), 1, g2), 2, g2), x))
        # real dualities: swaps lam1 <-> lam2 and mu1 <-> mu2; each has square -1 and maps x to x~
        for _ in range(3):
            z = mpf(random.random()*3 + 0.2) * random.choice([1, -1])
            # mu1(w g2) = z mu2(w), mu2(w g2) = -mu1(w)/z ; lam determined by (g1, g2)(x + x~) = x~ - x
            Pm = matrix([[m1[0], m2_[0]], [m1[1], m2_[1]]]); N = matrix([[z*m2_[0], -m1[0]/z], [z*m2_[1], -m1[1]/z]])
            G2 = N * inverse(Pm); g2 = ((G2[0, 0], G2[0, 1]), (G2[1, 0], G2[1, 1]))
            # (x + x~)(w g) = lam1(w g1) z^2 mu2^2 must equal -(x - x~) = -lam2 mu2^2: lam1(w g1) = -lam2(w)/z^2
            # (x - x~)(w g) = lam2(w g1) mu1^2 / z^2 must equal (x + x~) = lam1 mu1^2: lam2(w g1) = z^2 lam1(w)
            Pl = matrix([[l1[0], l2[0]], [l1[1], l2[1]]]); Nl = matrix([[-l2[0]/z**2, z**2*l1[0]], [-l2[1]/z**2, z**2*l1[1]]])
            G1 = Nl * inverse(Pl); g1 = ((G1[0, 0], G1[0, 1]), (G1[1, 0], G1[1, 1]))
            errs += [fabs(mdet(g1) - 1), fabs(mdet(g2) - 1)]
            for gg in (g1, g2):
                sq = mmul(gg, gg); errs.append(max(fabs(sq[0][0] + 1), fabs(sq[1][1] + 1), fabs(sq[0][1]), fabs(sq[1][0])))
            errs.append(maxdiff(c_act(c_act(c_act(x, 0, g1), 1, g2), 2, g2), xt))
    worst = max(worst, max(mpf(abs(z)) for z in errs))
print("rank-one factorizations, D_S = rho_S(-i), stabilizers, real dualities (square -1 in the non-BPS case),"
      " (1, D_T) x = -sgn(Delta) x and Freudenthal rotations at %d BPS and %d non-BPS charges; largest error %s"
      % (nB, nN, mp.nstr(worst, 3)))

# ---------- 4. all elements of SL2(Z)^2 mapping x to x~, BPS, k <= KB ----------
KB = 40
tot = {0: 0, 2: 0, 4: 0}; stab = {2: 0, 4: 0}; okG = True
for k in range(1, KB + 1):
    for s0 in orbits(4*k*k):
        d = dual(s0)
        if d is None: continue
        A1 = [tr(M) for M in aut(QS(s0))]; A2 = [tr(M) for M in aut(QT(s0))]
        assert len(A1) in (2, 4, 6) and len(A2) in (2, 4, 6)
        st = [(p, q) for p in A1 for q in A2 if act(p, q, s0) == s0]
        to = [(p, q) for p in A1 for q in A2 if act(p, q, s0) == d]
        stab[len(st)] += 1; tot[len(to)] += 1
        cs, ct = cont(QS(s0)) == k, cont(QT(s0)) == k
        okG &= len(st) == (4 if ct else 2)
        okG &= len(to) == (0 if not cs else (4 if ct else 2))
        if cs:
            DSi = tuple(tuple(z // (2*k) for z in row) for row in RS(s0))
            DTi = tuple(tuple(z // (2*k) for z in row) for row in RT(s0)) if ct else None
            want = {(DSi, I2), (DSi, MI)} | ({(tuple(tuple(-z for z in r) for r in DSi), DTi),
                                              (tuple(tuple(-z for z in r) for r in DSi), tuple(tuple(-z for z in r) for r in DTi))} if ct else set())
            okG &= set(to) == want
print("BPS, k <= %d: the stabilizer in SL2(Z)^2 is {(1, +-1)}, or {(1, +-1), (-1, +-D_T)} when cont(Q_T) = k;"
      " the elements taking x to x~ are (D_S, +-1), and also (-D_S, +-D_T) when cont(Q_T) = k:" % KB, okG)
print("   orbits by number of elements taking x to x~:", tot, " by stabilizer order:", stab)

# ---------- 5. examples ----------
def stu_U(s0):
    x8 = cube(s0); d = dual(s0)
    if d is None: return None
    d8 = cube(d)
    if cDet(x8) < 0: return stu_canon_def(x8) == stu_canon_def(d8)
    return stu_canon_sq(x8) == stu_canon_sq(d8)
def U(s0): d = dual(s0); return d is not None and canon(d) == canon(s0)
x = (1, 0, -1, 0, -2, 0)
print("Delta = 16 example:", x, "Delta", Delta(x), "dual", dual(x), "Q_S", QS(x), "Q_T", QT(x),
      "ST^2-U-dual", U(x), "STU-U-dual", stu_U(x), "canonical forms", canon(x), canon(dual(x)))
print("   x + i x~ = (u1 + 2i v1)(u2 + i v2)^2:", bps_charge((1, 0), ((1, 0), (0, 2)), ((1, 0), (0, 1)), 1) == (x, dual(x)),
      "; the dual lies in the orbit of Re((2u1 + i v1)(u2 + i v2)^2) =", bps_charge((1, 0), ((2, 0), (0, 1)), ((1, 0), (0, 1)), 1)[0],
      canon(bps_charge((1, 0), ((2, 0), (0, 1)), ((1, 0), (0, 1)), 1)[0]) == canon(dual(x)))
y1 = (1, 0, -1, 0, -1, 0)
print("Delta = 4 example:", y1, "Delta", Delta(y1), "dual", dual(y1), "= Im((u1 + i v1)(u2 + i v2)^2):",
      bps_charge((1, 0), ((1, 0), (0, 1)), ((1, 0), (0, 1)), 1) == (y1, dual(y1)), "U-dual", U(y1), "D_S", tuple(tuple(z // 2 for z in r_) for r_ in RS(y1)))
for dat in [(5, 2, 1, 0, 2, 2), (2, 1, 2, 1, 1, 1)]:
    y, yd = non_charge(*dat); e1, t1, e2, t2_, A_, B_ = dat
    gS = ((-t1, e1), (-(1 + t1*t1)//e1, t1)); gT = ((-t2_, e2), (-(1 + t2_*t2_)//e2, t2_))
    print("   datum", dat, ": g_S =", gS, "g_T =", gT, "(g_S, g_T) x = x~:", act(gS, gT, y) == yd, "(g_S, -g_T) x = x~:", act(gS, MI if gT == I2 else tuple(tuple(-z for z in r_) for r_ in gT), y) == yd)
okF = True
for r in range(1, 13):
    for t2 in range(1, 13):
        y = (1, 0, -r*r, 0, -t2, 0)
        okF &= Delta(y) == 4*r*r*t2*t2 and (dual(y) is not None) == (t2 % r == 0) and QT(y) == (t2, 0, r*r*t2)
        if t2 % r == 0: okF &= bps_charge((1, 0), ((1, 0), (0, t2 // r)), ((1, 0), (0, r)), 1) == (y, dual(y))
        if dual(y) is not None:
            okF &= U(y) == (r == t2) and stu_U(y) == (r == 1 or r == t2)
print("D6-D2 family (1, 0, -r^2, 0, -t, 0), r, t <= 12: Delta = 4r^2t^2, Q_T = t(u^2 + r^2 v^2), dual integral iff r | t, and then"
      " x + i x~ = (u1 + i(t/r)v1)(u2 + i r v2)^2;"
      " ST^2-U-dual iff r = t, STU-U-dual iff r = 1 or r = t:", okF)
okD = True
for p0 in range(-12, 13):
    for q0 in range(-12, 13):
        if p0*q0 == 0: continue
        y = (p0, 0, 0, 0, 0, q0)
        okD &= Delta(y) == -p0*p0*q0*q0 and dual(y) is not None and U(y) == (abs(p0) == abs(q0))
print("D0-D6 charges (p0, 0, 0, 0, 0, q0), |p0|, |q0| <= 12: Delta = -p0^2 q0^2, dual integral, U-dual iff |p0| = |q0|:", okD)
for dat in [(1, 0, 1, 0, 2, 4), (1, 0, 1, 0, 4, 2), (3, 1, 1, 0, 2, 2), (5, 2, 1, 0, 2, 2), (2, 1, 2, 1, 1, 1)]:
    y, yd = non_charge(*dat)
    print("   non-BPS datum (e1, t1, e2, t2, A, B) =", dat, "-> x =", y, "Delta", Delta(y), "dual", dual(y), "= data dual",
          yd == dual(y), "U-dual", U(y), "gcd(x+x~), gcd(x-x~) =", cont([p + q for p, q in zip(y, yd)]), cont([p - q for p, q in zip(y, yd)]),
          "canon", canon(y))
print("   D0-D6 charge (-2, 0, 0, 0, 0, 1) lies in the orbit of the datum (1, 0, 1, 0, 2, 4) or its dual:",
      canon((-2, 0, 0, 0, 0, 1)) in (canon(non_charge(1, 0, 1, 0, 2, 4)[0]), canon(non_charge(1, 0, 1, 0, 4, 2)[0])))

# ---------- 6. comparison with the STU and T^3 models ----------
from t3_classes import classes_bps as t3_bps, classes_non as t3_non, canon as t3_canon
from t3_cubic import dual as t3_dual
KC = 30
okS = True; nS = [0, 0, 0]
for k in range(1, KC + 1):
    for s0 in orbits(4*k*k):
        if dual(s0) is None: continue
        want = cont(QS(s0)) == k or cont(QT(s0)) == k
        okS &= stu_U(s0) == want; nS[0] += 1; nS[1] += want; nS[2] += U(s0)
    for s0 in orbits(-k*k):
        if dual(s0) is None: continue
        okS &= stu_U(s0) == U(s0)
print("ST^2 orbits with integral dual, k, m <= %d: BPS pairs are STU-U-dual iff cont(Q_S) = k or cont(Q_T) = k;"
      " non-BPS pairs are STU-U-dual iff ST^2-U-dual:" % KC, okS, " (BPS: %d orbits, %d STU-U-dual, %d ST^2-U-dual)" % tuple(nS))
okT = True; nT = 0
for k in range(1, KC + 1):
    for f4 in list(t3_bps(k)) + list(t3_non(k)):
        d4 = t3_dual(f4)
        if d4 is None: continue
        a4, b4, c4, dd4 = f4; x = (a4, b4, c4, b4, c4, dd4)
        assert dual(x) == (d4[0], d4[1], d4[2], d4[1], d4[2], d4[3])
        u3 = t3_canon(d4) == t3_canon(f4)
        okT &= U(x) == u3 and stu_U(x) == u3; nT += 1
print("T^3 orbits with integral dual, k, m <= %d (%d orbits): T^3-, ST^2- and STU-U-duality agree:" % (KC, nT), okT)
