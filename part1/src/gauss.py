"""Arithmetic in K = Q(sqrt D) with elements x + y*w, w = sqrt(D) stored as pairs of Fractions.
For the BPS case D = -4k^2 we use w = i (so sqrt D = 2k i); generally we store (x, y) meaning x + y*r, r^2 = R."""
from fractions import Fraction as Fr
from lat import Lat

class Field:
    def __init__(self, R):
        self.R = R           # r^2 = R; elements x + y r
    def mul(self, u, v):
        return (u[0]*v[0] + self.R*u[1]*v[1], u[0]*v[1] + u[1]*v[0])
    def conj(self, u): return (u[0], -u[1])
    def norm(self, u): return u[0]*u[0] - self.R*u[1]*u[1]
    def inv(self, u):
        n = self.norm(u); return (u[0]/n, -u[1]/n)
    def div(self, u, v): return self.mul(u, self.inv(v))
    def scale(self, k, L):
        return Lat([self.mul(k, v) for v in L.basis()])
    def prod(self, L, M):
        return Lat([self.mul(u, v) for u in L.basis() for v in M.basis()])
    def colon(self, S, M):
        """(S : M) = {x : x M subset S} = intersection over basis m of m^{-1} S."""
        out = None
        for m in M.basis():
            L = self.scale(self.inv(m), S)
            out = L if out is None else (out & L)
        return out
    def mult_ring(self, I):
        return self.colon(I, I)
