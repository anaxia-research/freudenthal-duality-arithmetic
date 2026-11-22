"""Second terms predicted by the poles at s=2/3 (F_+, F_-) and s=1/3 (E_-), compared with the exact partial sums at 10^7."""
import mpmath as mp
mp.mp.dps = 30
L = lambda s: (mp.zeta(s, 0.25) - mp.zeta(s, 0.75)) / 4**s
cF = lambda y: 1 + y*(1 - y)**2/(1 - y**3)
cE = lambda y: 1 + y*(1 - y**2)/(1 + y**3)
N = 10**7
exact = {'F+': 13288838, 'E+': 7854006, 'F-': 236350648, 'E-': 7194}
z2, z3 = mp.zeta(2), mp.zeta(3)
c1 = mp.mpf(8)/7*z2/z3
Zf = lambda s: mp.zeta(3*s - 1)/mp.zeta(3*s)*cF(2**(-s))
c0 = (2*mp.euler - 1)*Zf(1) + mp.diff(Zf, 1)
c2 = mp.zeta(1.5)*L(1.5)/z3*cE(2**-0.5)
# residues of F(s) x^s / s
sF = mp.zeta(mp.mpf(2)/3)*L(mp.mpf(2)/3)/(2*L(2))            # F_+ at s=2/3
sFm = mp.zeta(mp.mpf(2)/3)**2*cF(2**(-mp.mpf(2)/3))/(2*z2)     # F_- at s=2/3
sEm = mp.zeta(mp.mpf(2)/3)*mp.pi/4*cE(2**(-mp.mpf(1)/3))/z2       # E_- at s=1/3
x = mp.mpf(N)
pred = {
 'F+': (mp.mpf(4)/3*x, mp.mpf(4)/3*x + sF*x**(mp.mpf(2)/3)),
 'E+': (mp.pi/4*x, mp.pi/4*x),
 'F-': (c1*x*mp.log(x) + c0*x, c1*x*mp.log(x) + c0*x + sFm*x**(mp.mpf(2)/3)),
 'E-': (c2*mp.sqrt(x), c2*mp.sqrt(x) + sEm*x**(mp.mpf(1)/3)),
}
print("second-term constants: F+ %s  F- %s  E- %s" % (mp.nstr(sF, 10), mp.nstr(sFm, 10), mp.nstr(sEm, 10)))
for k in ('F+', 'E+', 'F-', 'E-'):
    a, b = pred[k]
    print(f"{k}: exact {exact[k]}  main {mp.nstr(a, 12)}  main+second {mp.nstr(b, 12)}")
