# Checks of statements used in the text of Part II.
import sympy as sp
from mpmath import mp, zeta, mpf, sqrt, log, euler, diff
from cubic import Det, T, act, hess
a,b,c,d=sp.symbols('a b c d')
H=sp.Matrix([[b**2-a*c,(b*c-a*d)/2],[(b*c-a*d)/2,c**2-b*d]])
J=sp.Matrix([[0,1],[-1,0]])
R=sp.Matrix([[b*c-a*d,2*(a*c-b**2)],[2*(c**2-b*d),a*d-b*c]])
print("R = -2HJ:", sp.simplify(R+2*H*J)==sp.zeros(2))
print("det(2HJ) = Delta:", sp.expand((2*H*J).det()+Det((a,b,c,d)))==0)
# g^T J g = det(g) J
x,y,z,w=sp.symbols('x y z w'); g=sp.Matrix([[x,y],[z,w]])
print("g^T J g = det g J:", sp.simplify(g.T*J*g-g.det()*J)==sp.zeros(2))
# traces of D rho, D rho^2 using D^2=-1
D=sp.Matrix([[0,1],[-1,0]]); r=-sp.Rational(1,2)*sp.eye(2)+sp.sqrt(3)/2*D
print("tr D rho, tr D rho^2:", sp.simplify((D*r).trace()), sp.simplify((D*r*r).trace()))
# examples
def show(f):
    t=T(f); De=-Det(f); print(f, "Delta",De,"T",t,"hess",hess(f))
for f in [(1,0,-1,0),(1,0,-4,0),(8,0,-2,0),(1,0,0,1),(2,0,0,1),(-4,-2,-1,0),(-125,-50,-20,-7)]:
    show(f)
mp.dps=30
c1=mpf(8)/7*zeta(2)/zeta(3)
lc=3*diff(zeta,2)/zeta(2)-3*diff(zeta,3)/zeta(3)+log(2)/14
c0=(2*euler-1)*c1+c1*lc
print("c1",c1,"c0",c0)
# single-slot lemma F_f(w1 R, w2, w3) = F_T(w1,w2,w3), and R_T = Delta R_f
def cubeof(f):
    a,b,c,d=f; x={}
    for A in (0,1):
        for B in (0,1):
            for C in (0,1):
                x[(A,B,C)]=[a,b,c,d][A+B+C]
    return x
f=(a,b,c,d); Tf=T(f); X=cubeof(f); XT=cubeof(Tf)
ok=True
for A in (0,1):
    for B in (0,1):
        for C in (0,1):
            # F_f(e_A R, e_B, e_C) = sum_A' R[A,A'] x[A',B,C]
            lhs=sum(R[A,Ap]*X[(Ap,B,C)] for Ap in (0,1))
            ok = ok and sp.expand(lhs-XT[(A,B,C)])==0
print("single slot:",ok)
RT=sp.Matrix([[Tf[1]*Tf[2]-Tf[0]*Tf[3],2*(Tf[0]*Tf[2]-Tf[1]**2)],[2*(Tf[2]**2-Tf[1]*Tf[3]),Tf[0]*Tf[3]-Tf[1]*Tf[2]]])
print("R_T = Delta R:", sp.simplify(RT+Det(f)*R)==sp.zeros(2))
# rotation formulas of Theorem 3.4, numerically with 30 digits at random charges and angles
import random
from mpmath import mpf, cos, sin, cosh, sinh, sqrt as msqrt, matrix
random.seed(1)
def actR(g, f):
    # (g.f)(w) = f(w g) for row vectors w; coefficients (a,b,c,d) of a u^3+3b u^2v+3c uv^2+d v^3
    (p, q), (r, s) = g  # w g = (p u + r v, q u + s v)
    a, b, c, d = f
    # trilinear: x_{ABC} -> sum g_{A A'} g_{B B'} g_{C C'} x_{A'B'C'}
    X = lambda i, j, k: [a, b, c, d][i + j + k]
    G = ((p, q), (r, s))
    Y = lambda i, j, k: sum(G[i][I]*G[j][J]*G[k][K]*X(I, J, K) for I in (0, 1) for J in (0, 1) for K in (0, 1))
    return (Y(0, 0, 0), Y(0, 0, 1), Y(0, 1, 1), Y(1, 1, 1))
worst = 0
for trial in range(200):
    f = tuple(random.randint(-9, 9) for _ in range(4))
    D = -Det(f)
    if D == 0: continue
    t = T(f); s = msqrt(abs(D)); ft = tuple((1 if D > 0 else -1)*x/s for x in t)
    a_, b_, c_, d_ = f
    Dm = ((mpf(a_*d_ - b_*c_)/s, mpf(2*(b_*b_ - a_*c_))/s), (mpf(2*(b_*d_ - c_*c_))/s, mpf(b_*c_ - a_*d_)/s))
    th = mpf(random.random()*6)
    if D > 0:
        M = ((cos(th) + sin(th)*Dm[0][0], sin(th)*Dm[0][1]), (sin(th)*Dm[1][0], cos(th) + sin(th)*Dm[1][1]))
        want = tuple(cos(3*th)*x - sin(3*th)*y for x, y in zip(f, ft))
    else:
        M = ((cosh(th) + sinh(th)*Dm[0][0], sinh(th)*Dm[0][1]), (sinh(th)*Dm[1][0], cosh(th) + sinh(th)*Dm[1][1]))
        want = tuple(cosh(3*th)*x + sinh(3*th)*y for x, y in zip(f, ft))
    got = actR(M, f)
    worst = max(worst, max(abs(x - y)/(1 + abs(y)) for x, y in zip(got, want)))
    assert max(abs(x - y)/(1 + abs(y)) for x, y in zip(got, ft if False else want)) < mpf(10)**-20
print("rotation formulas: max relative error", worst)
