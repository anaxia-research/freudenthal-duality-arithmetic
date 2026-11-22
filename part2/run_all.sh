#!/bin/sh
# Reproduces every computation of Part II. Run from this directory:  sh run_all.sh
# Requires Python 3 with numpy, sympy and mpmath.
set -e
cd src
python3 ident.py              > ../logs/ident.log         # covariants, Det(T) = Det^3, T(T) = -Delta^2 f, syzygy, R f = -Delta T
python3 paper_checks.py       > ../logs/paper_checks.log  # R = -2HJ, single slot, R_T = Delta R, rotations, examples, c0
python3 tables.py 200         > ../logs/tables200.log     # Hessian-reduction enumeration = Theorem E, k, m <= 200
python3 tables.py 40          > ../logs/tables.log        # logs/table.txt (Table 2 of the paper)
python3 generate.py 40        > ../logs/generate.log      # normal forms from (e, L, alpha) and (e, t, A, B) = enumeration
python3 criteria.py 200       > ../logs/criteria.log      # Theorems C and D, g_f, primitive Hessians
python3 box.py 30 40          > ../logs/box.log           # completeness: all forms with entries in [-30, 30]
python3 series.py 3000        > ../logs/series.log        # Dirichlet series of Theorem F
python3 asymptotics.py 10000000 > ../logs/asymptotics.log # partial sums to 10^7; logs/counts.dat (Figure 3)
python3 secondary.py          > ../logs/secondary.log     # second terms of Remark 9.5
echo "all checks passed"
