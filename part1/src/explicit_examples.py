"""Independent check of every explicit cube quoted in the examples section."""
from math import gcd, isqrt
from functools import reduce
from cube import Det, T, form, act
def g(v): return reduce(gcd, [abs(t) for t in v])
def cont(Q): return reduce(gcd, [abs(t) for t in Q])
def dual(x):
    D = -Det(x)
    s = isqrt(abs(D)); assert s*s == abs(D)
    t = T(x)
    assert all(c % s == 0 for c in t)
    return tuple((1 if D > 0 else -1)*c//s for c in t)
def Rmat(Q):
    p, q, r = Q
    return ((-q, -2*p), (2*r, q))
def show(name, x):
    D = -Det(x); xd = dual(x)
    print(name, "x =", x, "Delta =", D, "gcd", g(x), "gcdT", g(T(x)))
    print("   xt =", xd, "Delta(xt) =", -Det(xd), "gcd", g(xd), "gcdT", g(T(xd)))
    print("   forms x ", [form(x, j) for j in range(3)], "contents", [cont(form(x, j)) for j in range(3)])
    print("   forms xt", [form(xd, j) for j in range(3)])
    print("   gcd(x+xt) =", g([a+b for a, b in zip(x, xd)]), "gcd(x-xt) =", g([a-b for a, b in zip(x, xd)]))
    return xd
# Theorem A on a random example
x = (2, -1, 3, 0, 1, 4, -2, 5)
for j in range(3):
    assert act(x, j, Rmat(form(x, j))) == T(x)
print("Theorem A holds for", x)
# Delta = 4
x4 = (1, 0, 0, -1, 0, -1, -1, 0)
print("Det", Det(x4))
xd = show("Delta=4", x4)
for j in range(3):
    Q = form(x4, j); R = Rmat(Q); Phi = tuple(tuple(c//2 for c in row) for row in R)
    print("   Phi_%d =" % (j+1), Phi, "Phi x =", act(x4, j, Phi) == xd)
print()
x256a = (-4, -6, -2, 2, -2, 2, 0, 2); show("Delta=256 a", x256a)
x256b = (0, 1, 4, 0, 4, 0, 0, -4); show("Delta=256 b", x256b)
x576 = (-4, -9, -2, 3, -2, 3, 0, 3); show("Delta=576", x576)
for m in (1, 2, 3):
    show("D0D6 m=%d" % m, (1, 0, 0, 0, 0, 0, 0, m))
J = ((0, 1), (-1, 0)); Ji = ((0, -1), (1, 0))
for m, x in ((2, (-1, 1, 0, 0, 0, 0, -1, -1)), (5, (-2, 1, 0, 0, 0, 0, -1, -2))):
    xd = show("nonBPS U-dual m=%d" % m, x)
    y = act(act(act(x, 0, Ji), 1, Ji), 2, J)
    print("   (J^-1,J^-1,J) x == xt:", y == xd)
print("\nD6-D2 family x = e000 + q3 e011 + q2 e101 + q1 e110")
import itertools
from math import isqrt
cnt = 0
for q1, q2, q3 in itertools.product(range(-12, 13), repeat=3):
    if q1*q2*q3 >= 0: continue
    x = (1, 0, 0, q3, 0, q2, q1, 0)
    assert -Det(x) == -4*q1*q2*q3
    assert T(x) == (0, 2*q2*q3, 2*q1*q3, 0, 2*q1*q2, 0, 0, 2*q1*q2*q3)
    assert [form(x, j) for j in range(3)] == [(-q3, 0, q1*q2), (-q2, 0, q1*q3), (-q1, 0, q2*q3)]
    K2 = -q1*q2*q3; k = isqrt(K2)
    integral = k*k == K2 and all(c % (2*k) == 0 for c in T(x))
    pred = k*k == K2 and (q1*q2) % k == 0 and (q1*q3) % k == 0 and (q2*q3) % k == 0
    assert integral == pred
    if integral:
        cts = [reduce(gcd, [abs(c) for c in form(x, j)]) for j in range(3)]
        udual = k in cts
        assert udual == (k in (abs(q1), abs(q2), abs(q3)))
        cnt += 1
print("D6-D2 family claims verified for |q_i| <= 12;", cnt, "integral cases")
show("D6-D2 q=(-4,-4,-4)", (1, 0, 0, -4, 0, -4, -4, 0))
x16 = (1, 0, 0, -1, 0, -1, -4, 0)
print("q=(-4,-1,-1): Delta", -Det(x16), "T", T(x16))
