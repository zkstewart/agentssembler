#!/usr/bin/env python3
'''Each breakpoint of a block: reads (primary, MAPQ >= 20) spanning with >= 5 kb each side, the fewest in +-25 kb and at the point, plus the nearest gap. Pre-fill coordinates.
Use the AGP the PAFs were mapped against. Each PAF counts only for the scaffolds of its own haplotype when given as hN=path; a plain path counts for all.
Usage: bp_reads.py scaffolds.agp [h1=]a.paf.gz[,[h2=]b.paf.gz...] scaffold:start:end ...'''
import sys, gzip, collections, bisect
agp, pafs = sys.argv[1:3]; B = [x.split(":") for x in sys.argv[3:]]
G = collections.defaultdict(list); U = collections.defaultdict(list)
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] != "W": G[s[0]].append(int(s[1]))
    else: U[s[0]].append((int(s[1]), int(s[2]), s[5]))
P = {}
for b in B:
    for x in (int(b[1]), int(b[2])):
        for d in range(-25000, 25001, 5000): P.setdefault(b[0], {})[x + d] = 0
for f, h in ((p.split("=", 1)[1], "_" + p.split("=", 1)[0]) if "=" in p else (p, "") for p in pafs.split(",")):
    for l in gzip.open(f, "rt"):
        s = l.split("\t", 12)
        if not s[5].endswith(h) or s[5] not in P or int(s[11]) < 20 or "tp:A:P" not in l: continue
        a, c = int(s[7]) + 5000, int(s[8]) - 5000
        if c <= a: continue
        for q in P[s[5]]:
            if a <= q <= c: P[s[5]][q] += 1
print("block\tedge\tat_point\tmin_within_25kb\tnearest_gap_kb\tunitig")
for b in B:
    for x in (int(b[1]), int(b[2])):
        vals = [P[b[0]][x + d] for d in range(-25000, 25001, 5000)]
        g = min((abs(y - x) for y in G[b[0]]), default=None)
        u = next((r for r in U[b[0]] if r[0] <= x <= r[1]), None)
        print(f"{b[0]}:{int(b[1])/1e6:.2f}-{int(b[2])/1e6:.2f}\t{x/1e6:.2f}\t{P[b[0]][x]}\t{min(vals)}\t{g/1e3 if g is not None else 'NA'}\t{u[2] if u else 'gap'}")
