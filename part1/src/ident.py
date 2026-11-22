"""Symbolic identities for the STU Freudenthal triple system (sympy)."""
import sympy as sp
from cube import Det, T, form
a = sp.symbols('a0:8')
D = sp.expand(Det(a)); t = [sp.expand(x) for x in T(a)]
print("Det(T) = Det^3:", sp.expand(Det(t) - D**3) == 0)
print("T(T(a)) = -Det^2 a:", all(sp.expand(x + D**2*y) == 0 for x, y in zip(T(t), a)))
for i in range(3):
    q = [sp.expand(v) for v in form(t, i)]; q0 = form(a, i)
    lam = [sp.simplify(q[j] / (D*q0[j])) for j in range(3)]
    print(f"Q_{i+1}(T(a)) / (Det * Q_{i+1}(a)) =", lam)
# symplectic form and quartic norm: Delta = -Det; check <T(a), a> relation
