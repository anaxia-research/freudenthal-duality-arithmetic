"""T^3 model: charges are classical integral binary cubic forms f = (a,b,c,d) <-> a u^3 + 3b u^2 v + 3c u v^2 + d v^3,
i.e. triply symmetric cubes.  U-duality group SL2(Z) acting diagonally on the cube."""
from math import gcd, isqrt
from cube import Det as cDet, T as cT, form as cform, act as cact

def cube(f):
    a, b, c, d = f
    return (a, b, b, c, b, c, c, d)

def uncube(x):
    if all(isinstance(t, int) for t in x):
        assert x[1] == x[2] == x[4] and x[3] == x[5] == x[6], x
    else:
        import sympy
        assert all(sympy.expand(x[i] - x[j]) == 0 for i, j in ((1, 2), (1, 4), (3, 5), (3, 6)))
    return (x[0], x[1], x[3], x[7])

def Det(f):
    a, b, c, d = f
    return a*a*d*d - 6*a*b*c*d + 4*a*c**3 + 4*b**3*d - 3*b*b*c*c

def Delta(f): return -Det(f)

def hess(f):
    a, b, c, d = f
    return (b*b - a*c, a*d - b*c, c*c - b*d)

def T(f): return uncube(cT(cube(f)))

def act(g, f):
    """Diagonal action of g on the cube; equals f(g^T (u,v))."""
    x = cube(f)
    for i in range(3): x = cact(x, i, g)
    return uncube(x)

def qsub(Q, M):
    """Quadratic form Q(M (u,v)) for an integer matrix M (substitution)."""
    P, q, R = Q
    (m00, m01), (m10, m11) = M
    # u' = m00 u + m01 v, v' = m10 u + m11 v
    A = P*m00*m00 + q*m00*m10 + R*m10*m10
    B = 2*P*m00*m01 + q*(m00*m11 + m01*m10) + 2*R*m10*m11
    C = P*m01*m01 + q*m01*m11 + R*m11*m11
    return (A, B, C)

def tr(M): return ((M[0][0], M[1][0]), (M[0][1], M[1][1]))
def mul(A, B): return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def cont(v):
    g = 0
    for t in v: g = gcd(g, t)
    return g

def dual(f):
    """Freudenthal dual sgn(Delta) T / sqrt|Delta| if integral, else None."""
    D = Delta(f)
    n = isqrt(abs(D))
    if n*n != abs(D): return None
    t = T(f)
    if any(x % n for x in t): return None
    s = 1 if D > 0 else -1
    return tuple(s*x//n for x in t)

J = ((1, 0), (0, -1))
def phi(g):
    """hess(act(g, f)) == qsub(hess(f), phi(g))."""
    return mul(mul(J, tr(g)), J)
def phi_inv(M): return mul(mul(J, tr(M)), J)
