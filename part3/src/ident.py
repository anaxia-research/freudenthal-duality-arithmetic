"""Polynomial identities for ST^2 charges, in Z[a,b,c,e,f,g] (and generic matrices where they occur)."""
import sympy as sp
from st2 import cube, uncube, QS, QT, T, RS, RT, KEYS
from cube import Det as cDet, act as cact
a, b, c, e, f, g, u1, v1, u2, v2 = sp.symbols('a b c e f g u1 v1 u2 v2')
s = (a, b, c, e, f, g)
D = sp.expand(-cDet(cube(s)))
t = [sp.expand(x) for x in T(s)]
def ev(y):
    """y(w1, w2, w2) with w1 = (u1, v1), w2 = (u2, v2)."""
    x = cube(y); w1 = (u1, v1); w2 = (u2, v2)
    return sp.expand(sum(x[4*A + 2*B + C]*w1[A]*w2[B]*w2[C] for A in range(2) for B in range(2) for C in range(2)))
def qev(Q, w): return sp.expand(Q[0]*w[0]**2 + Q[1]*w[0]*w[1] + Q[2]*w[1]**2)
def matQ(M): return [[M[i][j] for j in range(2)] for i in range(2)]
# definitions of the face forms
x8 = cube(s)
MS = sp.Matrix(2, 2, lambda B, C: u1*x8[2*B + C] + v1*x8[4 + 2*B + C])
MT = sp.Matrix(2, 2, lambda A, C: u2*x8[4*A + C] + v2*x8[4*A + 2 + C])
print("Q_S(w) = -det x(w,.,.), Q_T(w) = -det x(.,w,.):", sp.expand(qev(QS(s), (u1, v1)) + MS.det()) == 0,
      sp.expand(qev(QT(s), (u2, v2)) + MT.det()) == 0)
dS = sp.expand(QS(s)[1]**2 - 4*QS(s)[0]*QS(s)[2]); dT = sp.expand(QT(s)[1]**2 - 4*QT(s)[0]*QT(s)[2])
print("disc Q_S = disc Q_T = -Delta:", sp.expand(dS + D) == 0, sp.expand(dT + D) == 0)
# covariance under (g1, g2)
G1 = sp.Matrix(2, 2, sp.symbols('p0:4')); G2 = sp.Matrix(2, 2, sp.symbols('r0:4'))
g1 = ((G1[0,0], G1[0,1]), (G1[1,0], G1[1,1])); g2 = ((G2[0,0], G2[0,1]), (G2[1,0], G2[1,1]))
y = cact(cact(cact(cube(s), 0, g1), 1, g2), 2, g2); ys = uncube(tuple(sp.expand(z) for z in y))
w1g = (G1[0,0]*u1 + G1[1,0]*v1, G1[0,1]*u1 + G1[1,1]*v1); w2g = (G2[0,0]*u2 + G2[1,0]*v2, G2[0,1]*u2 + G2[1,1]*v2)
print("(g.x)(w1,w2,w2) = x(w1 g1, w2 g2, w2 g2):", sp.expand(ev(ys) - ev(s).subs({u1: w1g[0], v1: w1g[1], u2: w2g[0], v2: w2g[1]}, simultaneous=True)) == 0)
print("Delta(g.x) = det(g1)^2 det(g2)^4 Delta(x):", sp.expand(-cDet(cube(ys)) - G1.det()**2*G2.det()**4*D) == 0)
print("Q_S(g.x)(w) = det(g2)^2 Q_S(x)(w g1):", sp.expand(qev(QS(ys), (u1, v1)) - G2.det()**2*qev(QS(s), w1g)) == 0)
print("Q_T(g.x)(w) = det(g1) det(g2) Q_T(x)(w g2):", sp.expand(qev(QT(ys), (u2, v2)) - G1.det()*G2.det()*qev(QT(s), w2g)) == 0)
ty = [sp.expand(z) for z in T(ys)]
ygT = uncube(tuple(sp.expand(z) for z in cact(cact(cact(cube(tuple(t)), 0, g1), 1, g2), 2, g2)))
print("T(g.x) = det(g1) det(g2)^2 g.T(x):", all(sp.expand(p - G1.det()*G2.det()**2*q) == 0 for p, q in zip(ty, ygT)))
# triple product identities
print("Delta(T(x)) = Delta^3:", sp.expand(-cDet(cube(tuple(t))) - D**3) == 0)
tt = T(tuple(t))
print("T(T(x)) = -Delta^2 x:", all(sp.expand(p + D**2*q) == 0 for p, q in zip(tt, s)))
print("Q_S(T(x)) = Delta Q_S(x), Q_T(T(x)) = Delta Q_T(x):",
      all(sp.expand(p - D*q) == 0 for p, q in zip(QS(tuple(t)), QS(s))), all(sp.expand(p - D*q) == 0 for p, q in zip(QT(tuple(t)), QT(s))))
print("syzygy T(x)(w1,w2,w2)^2 + Delta x(w1,w2,w2)^2 = 4 Q_S(w1) Q_T(w2)^2:",
      sp.expand(ev(tuple(t))**2 + D*ev(s)**2 - 4*qev(QS(s), (u1, v1))*qev(QT(s), (u2, v2))**2) == 0)
# the duality matrices
R = RS(s); p, q, r = QS(s)
H = sp.Matrix([[p, q/sp.Integer(2)], [q/sp.Integer(2), r]]); Jm = sp.Matrix([[0, 1], [-1, 0]])
print("R_S = -2 H_S J:", sp.simplify(sp.Matrix(R) + 2*H*Jm) == sp.zeros(2))
RR = sp.Matrix(R)*sp.Matrix(R)
print("R_S^2 = -Delta * 1, det R_S = Delta:", sp.expand(RR[0,0] + D) == 0 and sp.expand(RR[0,1]) == 0 and sp.expand(RR[1,0]) == 0 and sp.expand(RR[1,1] + D) == 0,
      sp.expand(sp.Matrix(R).det() - D) == 0)
print("T(x) = R_S acting on the S index:", all(sp.expand(p1 - p2) == 0 for p1, p2 in zip(uncube(cact(cube(s), 0, R)), t)))
Rt = RT(s)
print("T(x) = R_T acting on one T index (either one):", all(sp.expand(p1 - p2) == 0 for p1, p2 in zip(cact(cube(s), 1, Rt), cube(tuple(t)))),
      all(sp.expand(p1 - p2) == 0 for p1, p2 in zip(cact(cube(s), 2, Rt), cube(tuple(t)))))
both = cact(cact(cube(s), 1, Rt), 2, Rt)
print("(1, R_T) x = -Delta x:", all(sp.expand(p1 + D*p2) == 0 for p1, p2 in zip(both, cube(s))))
print("Delta =", sp.factor(D))
