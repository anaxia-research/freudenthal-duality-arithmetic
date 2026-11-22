#!/bin/sh
# Reproduces every computation of Part I. Run from this directory:  sh run_all.sh
# Requires Python 3 with numpy, sympy and mpmath. Total time: under an hour on a laptop,
# most of it in the two box checks at the end.
set -e
cd src
python3 ident.py                 > ../logs/ident.log          # Det(T) = Det^3, T(T) = -Det^2 x, Q_j(T) = Delta Q_j
python3 rmat.py                  > ../logs/rmat.log           # T(x) = R_j x on each index (Theorem A)
python3 bps.py 40                > ../logs/bps40.log          # BPS balanced triples, k <= 40
python3 split.py 12              > ../logs/split12.log        # non-BPS balanced triples, m <= 12
python3 splitrun.py 13 15 16 25 17 18 20 24 32 14 19 21 22 23 26 27 28 29 30 31 33 34 35 36 37 38 39 40 > ../logs/split_more.log
python3 check_counts.py          > ../logs/check_counts.log   # enumeration = Theorems E and F, k, m <= 40
python3 formula.py               > ../logs/formula.log        # finite formula = Dirichlet series, k <= 40
python3 tables.py                > ../logs/tables.log         # logs/table_bps.txt, logs/table_non.txt
python3 tricube.py 24            > ../logs/tricube24.log      # Bhargava-free canonical forms, BPS
python3 splitcube.py 24          > ../logs/splitcube24.log    # Bhargava-free canonical forms, non-BPS
python3 criteria.py 16 24        > ../logs/criteria.log       # explicit criteria of Theorems C and D
python3 triality_check.py 16     > ../logs/triality.log       # triality criterion
python3 nonbps_check.py 12       > ../logs/nonbps_check.log
python3 series.py 100000         > ../logs/series.log
python3 asymptotics.py 100000    > ../logs/asymptotics.log    # partial sums against the main terms; logs/counts.dat
python3 maincoef.py              > ../logs/maincoef.log       # residue polynomials
python3 figdata.py 100000        > ../logs/figdata.log        # logs/share.dat
python3 d0d6.py                  > ../logs/d0d6.log
python3 examples.py              > ../logs/examples.log
python3 examples2.py             > ../logs/examples2.log
python3 explicit_examples.py     > ../logs/explicit_examples.log
python3 fived.py                 > ../logs/fived.log          # five-dimensional Jordan dual
python3 boxcompare.py 2 8        > ../logs/boxcompare.log     # completeness in the box [-2,2]^8, BPS
python3 boxcompare_non.py 2 1 8  > ../logs/boxcompare_non.log # completeness in the box [-2,2]^8, non-BPS
echo "all checks passed"
