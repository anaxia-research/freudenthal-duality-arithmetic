"""ST^2 model: a charge is a 2x2x2 cube symmetric in its last two indices,
    s = (a, b, c, e, f, g) <-> x000 = a, x001 = x010 = b, x011 = c, x100 = e, x101 = x110 = f, x111 = g,
with x(w1, w2, w3) = sum x_ABC w1_A w2_B w3_C.  Charges (p0, pS, pT, q0, qS, qT) = (a, -e, -b, g, c, f).
The U-duality group SL2(Z) x SL2(Z) acts by (g1, g2) . x (w1, w2, w3) = x(w1 g1, w2 g2, w3 g2).
Quadratic forms are triples (P, q, R) <-> P u^2 + q uv + R v^2."""
from math import gcd, isqrt
from cube import Det as cDet, T as cT, act as cact

KEYS = [(0, 0, 0), (0, 0, 1), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 1)]

def cube(s):
    a, b, c, e, f, g = s
    return (a, b, b, c, e, f, f, g)

def uncube(x):
    assert x[1] == x[2] and x[5] == x[6], x
    return (x[0], x[1], x[3], x[4], x[5], x[7])

def entry(s, A, B, C):
    return cube(s)[4*A + 2*B + C]

def Delta(s): return -cDet(cube(s))

def QS(s):
    """Q_S(w) = -det [x(w, e_B, e_C)]_{B,C}."""
    a, b, c, e, f, g = s
    return (b*b - a*c, 2*b*f - a*g - c*e, f*f - e*g)

def QT(s):
    """Q_T(w) = -det [x(e_A, w, e_C)]_{A,C}."""
    a, b, c, e, f, g = s
    return (b*e - a*f, c*e - a*g, c*f - b*g)

def T(s): return uncube(cT(cube(s)))

def act(g1, g2, s):
    x = cact(cube(s), 0, g1)
    x = cact(x, 1, g2); x = cact(x, 2, g2)
    return uncube(x)

def act1(g, s, slot):
    """g acting on one slot only (slot 0, 1 or 2); returns an 8-entry cube."""
    return cact(cube(s), slot, g)

def cont(v):
    r = 0
    for t in v: r = gcd(r, t)
    return abs(r)

def sqrtabs(D):
    n = isqrt(abs(D)); return n if n*n == abs(D) else None

def dual(s):
    """Freudenthal dual sgn(Delta) T(x)/sqrt|Delta| if integral, else None."""
    D = Delta(s); n = sqrtabs(D)
    if not n: return None
    t = T(s)
    if any(v % n for v in t): return None
    sg = 1 if D > 0 else -1
    return tuple(sg*v//n for v in t)

def mul(A, B): return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def det2(A): return A[0][0]*A[1][1] - A[0][1]*A[1][0]
def tr(M): return ((M[0][0], M[1][0]), (M[0][1], M[1][1]))
I2 = ((1, 0), (0, 1)); S2 = ((0, -1), (1, 0)); MI = ((-1, 0), (0, -1)); J = ((0, 1), (-1, 0))

def RS(s):
    """T(x) = R_S acting on the first index; R_S = -2 H_S J with H_S the matrix of Q_S."""
    p, q, r = QS(s)
    return ((q, -2*p), (2*r, -q))

def RT(s):
    """T(x) = R_T acting on one of the two T indices."""
    p, q, r = QT(s)
    return ((q, -2*p), (2*r, -q))

# ---------- binary quadratic forms ----------
def qsub(Q, M):
    """Q(M (u,v)^T), i.e. Q(w M^T) for the row vector w."""
    P, q, R = Q; (m00, m01), (m10, m11) = M
    return (P*m00*m00 + q*m00*m10 + R*m10*m10,
            2*P*m00*m01 + q*(m00*m11 + m01*m10) + 2*R*m10*m11,
            P*m01*m01 + q*m01*m11 + R*m11*m11)

def qact(g, Q):
    """The form w -> Q(w g)."""
    return qsub(Q, tr(g))

def disc(Q): return Q[1]*Q[1] - 4*Q[0]*Q[2]

def egcd(a, b):
    if b == 0: return (a, 1, 0) if a >= 0 else (-a, -1, 0)
    g, x, y = egcd(b, a % b)
    return (g, y, x - (a // b) * y)

def reduce_def(Q):
    """Positive definite Q: (Qred, M) with qsub(Q, M) == Qred reduced, M in SL2(Z)."""
    M = I2
    while True:
        P, q, R = Q
        t = -((q + P) // (2 * P))
        if t:
            Tt = ((1, t), (0, 1)); Q = qsub(Q, Tt); M = mul(M, Tt); P, q, R = Q
        if P > R:
            Q = qsub(Q, S2); M = mul(M, S2); continue
        break
    P, q, R = Q
    if q == -P:
        Tt = ((1, 1), (0, 1)); Q = qsub(Q, Tt); M = mul(M, Tt)
    P, q, R = Q
    if P == R and q < 0:
        Q = qsub(Q, S2); M = mul(M, S2)
    return Q, M

SMALL = [((a, b), (c, d)) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) for d in (-1, 0, 1) if a*d - b*c == 1]
def aut(Q): return [A for A in SMALL if qsub(Q, A) == Q]

def normalizers(Q):
    """All M in SL2(Z) with qsub(Q, M) in normal shape (a finite set).  Definite Q: qsub(Q, M) = +-(reduced form).
    Square discriminant: qsub(Q, M) = (0, q', r) with 0 <= r < |q'|."""
    P, q, R = Q; D = disc(Q)
    if D < 0:
        sg = 1 if P > 0 else -1
        Qr, M = reduce_def(tuple(sg*t for t in Q))
        return [mul(M, A) for A in aut(Qr)]
    m = isqrt(D); assert m*m == D and m > 0
    out = []
    roots = [(1, 0), (-R, q)] if P == 0 else [(-q + m, 2*P), (-q - m, 2*P)]
    for (u0, v0) in roots:
        gg = gcd(u0, v0); u0 //= gg; v0 //= gg
        _, xx, yy = egcd(u0, v0)
        M = ((u0, -yy), (v0, xx))
        Q1 = qsub(Q, M); assert Q1[0] == 0
        q1, R1 = Q1[1], Q1[2]
        t = -(R1 // q1) if q1 > 0 else (R1 // (-q1))
        M2 = mul(M, ((1, t), (0, 1))); Q2 = qsub(Q, M2)
        assert Q2[0] == 0 and 0 <= Q2[2] < abs(Q2[1])
        out += [M2, mul(M2, MI)]
    return out

def qcanon(Q):
    """Canonical representative of the SL2(Z)-class of a definite or square-discriminant form."""
    return min(qsub(Q, M) for M in normalizers(Q))

def canon(s):
    """Canonical representative of the SL2(Z)^2-orbit of s (Delta != 0, -Delta not a nonsquare > 0):
    Q_S is brought to normal shape by the first factor and Q_T by the second; the stabilizers are finite."""
    best = None
    for M1 in normalizers(QS(s)):
        y = act(tr(M1), I2, s)
        for M2 in normalizers(QT(y)):
            z = act(I2, tr(M2), y)
            if best is None or z < best: best = z
    return best

# ---------- enumeration of all orbits with a given Delta ----------
def forms_with_disc(D):
    """Forms of discriminant D (D < 0, or D a positive square) meeting every SL2(Z)-class; definite: both signs."""
    out = []
    if D < 0:
        Pmax = isqrt(-D // 3) + 1
        for P in range(1, Pmax + 1):
            for q in range(-P + 1, P + 1):
                if (q*q - D) % (4*P): continue
                R = (q*q - D) // (4*P)
                if R < P or (P == R and q < 0): continue
                out.append((P, q, R)); out.append((-P, -q, -R))
    else:
        m = isqrt(D); assert m*m == D
        for q in (m, -m):
            for R in range(m): out.append((0, q, R))
    return out

def plane_basis(n):
    """Basis (A0, B0) of {v in Z^3 : v . n = 0} with A0 x B0 = n, for primitive n."""
    n1, n2, n3 = n
    if n1 == 0 and n2 == 0:
        A0, B0 = (1, 0, 0), (0, 1, 0)
    else:
        g, a, b = egcd(n1, n2)
        A0 = (n2 // g, -n1 // g, 0)
        B0 = (-n3 * a, -n3 * b, g)
    c = (A0[1]*B0[2]-A0[2]*B0[1], A0[2]*B0[0]-A0[0]*B0[2], A0[0]*B0[1]-A0[1]*B0[0])
    if c == tuple(n): return A0, B0
    assert c == tuple(-t for t in n), (n, c)
    return B0, A0

def plucker(Q):
    """The rows A = (a,b,c), B = (e,f,g) of a charge with Q_T = Q satisfy A x B = plucker(Q)."""
    return (-Q[2], Q[1], -Q[0])

def orbits(Dq):
    """Canonical representatives of all SL2(Z)^2-orbits of ST^2 charges with quartic norm Dq.
    For each normalized form Q of discriminant -Dq the charges with Q_T = Q are the bases C (A0; B0) of the
    sublattices of index cont(Q) of the plane lattice, and the first SL2(Z) acts on C from the left."""
    found = set()
    for Q in forms_with_disc(-Dq):
        nvec = plucker(Q)
        c = cont(nvec); n0 = tuple(t // c for t in nvec)
        A0, B0 = plane_basis(n0)
        for d1 in range(1, c + 1):
            if c % d1: continue
            d2 = c // d1
            for j in range(d2):
                A = tuple(d1*A0[i] + j*B0[i] for i in range(3))
                B = tuple(d2*B0[i] for i in range(3))
                s = A + B
                assert QT(s) == Q, (Q, s, QT(s))
                found.add(canon(s))
    return found

if __name__ == "__main__":
    import random
    random.seed(5)
    G = [((1,1),(0,1)), ((0,-1),(1,0)), ((2,1),(1,1)), ((1,0),(3,1)), ((1,-2),(1,-1)), MI]
    n = 0; nc = 0
    for _ in range(4000):
        s = tuple(random.randint(-7, 7) for _ in range(6)); D = Delta(s)
        g = I2; h = I2
        for _ in range(4): g = mul(g, random.choice(G)); h = mul(h, random.choice(G))
        y = act(g, h, s)
        assert Delta(y) == D and QS(y) == qact(g, QS(s)) and QT(y) == qact(h, QT(s))
        assert QS(act(I2, h, s)) == QS(s) and QT(act(g, I2, s)) == QT(s)
        assert disc(QS(s)) == -D and disc(QT(s)) == -D
        n += 1
        if D == 0 or (D < 0 and sqrtabs(D) is None): continue
        assert canon(s) == canon(y), (s, g, h)
        nc += 1
    print("covariance of Q_S, Q_T and disc = -Delta checked on", n, "charges; canonical form invariant on", nc)
