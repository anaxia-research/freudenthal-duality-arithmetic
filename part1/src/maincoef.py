"""Coefficients of the main terms: Res F(s) X^s/s = X^{s0} * P(log X), P of degree w-1."""
import mpmath as mp
from asymptotics import Fp_s, Ep_s, Fm_s, Em_s, main_term
mp.mp.dps = 30
for name, F, s0, w in (("F+", Fp_s, 2, 3), ("E+", Ep_s, 2, 2), ("F-", Fm_s, 2, 3), ("E-", Em_s, 1, 3)):
    Xs = [mp.mpf(10)**j for j in (3, 6, 9)][:w]
    vals = [main_term(F, s0, X)/X**s0 for X in Xs]
    A = mp.matrix([[mp.log(X)**i for i in range(w)] for X in Xs])
    c = mp.lu_solve(A, mp.matrix(vals))
    print(name, "P(L) =", " + ".join(f"{mp.nstr(c[i], 10)} L^{i}" for i in range(w)))
G = mp.catalan
print("check: F+ lead pi^2/(24 G^2) =", mp.nstr(mp.pi**2/(24*G**2), 10), "; E+ lead pi^2/(4G) =", mp.nstr(mp.pi**2/(4*G), 10),
      "; F- lead 2/pi^2 =", mp.nstr(2/mp.pi**2, 10), "; E- lead 5/(16 pi) =", mp.nstr(5/(16*mp.pi), 10))
