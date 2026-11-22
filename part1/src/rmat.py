import sympy as sp
from cube import Det, T, form, act
a = sp.symbols('a0:8')
t = [sp.expand(x) for x in T(a)]
for j in range(3):
    p, q, r = [sp.expand(v) for v in form(a, j)]
    found = None
    for R in ([[q, -2*p], [2*r, -q]], [[-q, 2*p], [-2*r, q]], [[q, 2*r], [-2*p, -q]], [[-q, -2*r], [2*p, q]],
              [[q, -2*r], [2*p, -q]], [[-q, 2*r], [-2*p, q]], [[q, 2*p], [-2*r, -q]], [[-q, -2*p], [2*r, q]]):
        b = act(a, j, R)
        if all(sp.expand(x - y) == 0 for x, y in zip(b, t)):
            found = R; break
    print("index", j+1, "T(a) = R_j a with R_j =", found)
