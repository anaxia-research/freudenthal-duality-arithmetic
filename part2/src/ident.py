"""Polynomial identities for classical integral binary cubic forms (T^3 charges)."""
import sympy as sp
from cubic import Det, Delta, hess, T, act
a, b, c, d, u, v = sp.symbols('a b c d u v')
f = (a, b, c, d); D = sp.expand(Delta(f)); t = [sp.expand(x) for x in T(f)]
ev = lambda F: F[0]*u**3 + 3*F[1]*u**2*v + 3*F[2]*u*v**2 + F[3]*v**3
p, q, r = hess(f)
Qh = p*u**2 - q*u*v + r*v**2     # Hessian form hat Q = -H/36
H = sp.diff(ev(f), u, 2)*sp.diff(ev(f), v, 2) - sp.diff(ev(f), u, v)**2
print("hat Q = -H/36:", sp.expand(Qh + H/36) == 0)
print("Det(T) = Det^3:", sp.expand(Det(t) - Det(f)**3) == 0)
print("T(T(f)) = -Delta^2 f:", all(sp.expand(x + D**2*y) == 0 for x, y in zip(T(t), f)))
pt, qt, rt = hess(t)
print("hess(T(f)) = Delta hess(f):", all(sp.expand(x - D*y) == 0 for x, y in zip((pt, qt, rt), (p, q, r))))
print("syzygy T(f)^2 + Delta f^2 = 4 hatQ^3:", sp.expand(ev(t)**2 + D*ev(f)**2 - 4*Qh**3) == 0)
R = ((-q, -2*p), (2*r, q))
print("R f = -Delta T(f) (diagonal action):", all(sp.expand(x + D*y) == 0 for x, y in zip(act(R, f), t)))
print("R^2 = Det * 1:", sp.expand(R[0][0]**2 + R[0][1]*R[1][0] - Det(f)) == 0 and sp.expand(R[0][0]*R[0][1] + R[0][1]*R[1][1]) == 0)
g = sp.Matrix(2, 2, sp.symbols('g0:4'))
gf = act(((g[0,0], g[0,1]), (g[1,0], g[1,1])), f)
print("(g f)(u,v) = f(g^T (u,v)):", sp.expand(ev(gf) - ev(f).subs({u: g[0,0]*u + g[1,0]*v, v: g[0,1]*u + g[1,1]*v}, simultaneous=True)) == 0)
pg, qg, rg = hess(gf)
Qg = pg*u**2 - qg*u*v + rg*v**2
print("hatQ(g f) = det(g)^2 hatQ(g^T (u,v)):", sp.expand(Qg - g.det()**2*Qh.subs({u: g[0,0]*u + g[1,0]*v, v: g[0,1]*u + g[1,1]*v}, simultaneous=True)) == 0)
print("Delta(g f) = det(g)^6 Delta(f):", sp.expand(Delta(gf) - g.det()**6*D) == 0)
