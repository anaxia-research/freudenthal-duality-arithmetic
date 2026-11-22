#!/bin/sh
# Reproduces every computation of Part III. Run from this directory:  sh run_all.sh
# Requires Python 3 with numpy, sympy and mpmath.
set -e
cd src
python3 st2.py                > ../logs/st2.log           # covariance of Q_S, Q_T; canonical form invariant under the group
python3 ident.py              > ../logs/ident.log         # face forms, covariance, Delta(T) = Delta^3, T(T) = -Delta^2 x, syzygy, R_S, R_T
python3 paper_checks.py       > ../logs/paper_checks.log  # explicit T, R_S covariance, iota, real types, rotations, stabilizers, examples
python3 tables.py 200         > ../logs/tables200.log     # enumeration of all orbits = Theorem E, k, m <= 200
python3 tables.py 40          > ../logs/tables.log        # logs/table.txt (Table 3 of the paper)
python3 generate.py 40        > ../logs/generate.log      # normal forms from (alpha, L_S, L_T) and (e1, t1, e2, t2, A, B) = enumeration
python3 criteria.py 200       > ../logs/criteria.log      # Theorems C and D, g_S and g_T, primitive face forms
python3 box.py 4 40           > ../logs/box.log           # completeness: all charges with entries in [-4, 4]
python3 series.py 3000        > ../logs/series.log        # Dirichlet series of Theorem F
python3 asymptotics.py 10000000 > ../logs/asymptotics.log # partial sums to 10^7, second terms; logs/counts.dat (Figure 3)
echo "all checks passed"
