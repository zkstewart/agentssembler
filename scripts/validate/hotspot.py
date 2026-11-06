#!/usr/bin/env python3
'''Contact profile of named regions: pairs per 50 kb bin along the region, partners by scaffold, and partner concentration. Usage: hotspot.py pairs.tsv.gz scaffold:start-end ...'''
import sys, gzip, collections
pa = sys.argv[1]; R = []
for r in sys.argv[2:]:
    c, x = r.split(":"); a, b = x.split("-"); R.append((c, int(a), int(b)))
P = {r: collections.Counter() for r in R}; Bn = {r: collections.Counter() for r in R}; T = {r: collections.Counter() for r in R}
for l in gzip.open(pa, "rt"):
    a, x, b, y = l.rstrip("\n").split("\t"); x = int(x); y = int(y)
    for r in R:
        for p, q, u, v in ((a, x, b, y), (b, y, a, x)):
            if p == r[0] and r[1] <= q < r[2]:
                P[r][u] += 1; Bn[r][q // 50000] += 1; T[r][(u, v // 500000)] += 1
for r in R:
    print(r, "total", sum(P[r].values()))
    print(" along region:", " ".join(f"{k*0.05:.2f}:{v}" for k, v in sorted(Bn[r].items())))
    print(" by scaffold:", " ".join(f"{k}:{v}" for k, v in P[r].most_common(8)))
    print(" top 500 kb partners:", " ".join(f"{k[0]}@{k[1]*0.5:.1f}:{v}" for k, v in T[r].most_common(8)))
