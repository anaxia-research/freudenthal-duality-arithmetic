"""Data for the share figure: actual proportions of U-dual pairs, the ratio of the full main terms
(residue polynomials, coefficients from maincoef.py) and the leading-order predictions of Corollary H."""
import numpy as np, math, sys
from series import arrays
from asymptotics import Eminus_array
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10**5
Fp, Ep, Fm = arrays(N); Em = Eminus_array(N)
cFp, cEp, cFm, cEm = (np.cumsum(a.astype(np.float64)) for a in (Fp, Ep, Fm, Em))
# main-term polynomials in L = log K (BPS, K = sqrt(X)/2) and L = log M (non-BPS, M = sqrt(X))
PFp = (2.949719306, 0.4740306734, 0.4901513261); PEp = (-0.01237773447, 2.693770504)
PFm = (0.3005327114, 0.6365235581, 0.2026423673); PEm = (0.6855584093, 0.7548976649, 0.09947183943)
ev = lambda P, L: sum(c*L**i for i, c in enumerate(P))
G = 0.915965594177219015
with open("../logs/share.dat", "w") as f:
    f.write("lgXB shareB mainB leadB lgXN sqrtXshareN mainN leadN\n")
    for t in np.linspace(1.0, math.log10(N), 60):
        K = int(round(10**t)); L = math.log(K)
        Xb = 4*K*K; Xn = K*K
        sB = cEp[K]/cFp[K]; mB = ev(PEp, L)/ev(PFp, L); lB = 12*G/math.log(Xb)
        sN = cEm[K]/cFm[K]; mN = (K*ev(PEm, L))/(K*K*ev(PFm, L)); lN = 5*math.pi/(32*K)
        f.write(f"{math.log10(Xb):.4f} {sB:.5f} {mB:.5f} {lB:.5f} {math.log10(Xn):.4f} {sN*K:.5f} {mN*K:.5f} {lN*K:.5f}\n")
print("wrote share.dat")
