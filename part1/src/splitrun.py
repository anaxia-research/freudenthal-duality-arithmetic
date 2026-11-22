import sys
from split import triples, dual_triple_class
for m in [int(x) for x in sys.argv[1:]]:
    r = triples(m)
    Fi = [x for x in r if x['fint']]
    Fe = 0
    for x in Fi:
        (f1, t1), (f2, t2), I3p = dual_triple_class(x)
        Fe += ((f1, t1) == (x['f'][0], x['t'][0]) and (f2, t2) == (x['f'][1], x['t'][1]) and I3p == x['I3'])
    print(m, len(r), len(Fi), Fe, flush=True)
