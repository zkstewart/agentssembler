#!/usr/bin/env python3
'''Per gap edge: reads reaching the edge (within slack, >=5 kb anchored, MAPQ>=20), how many have >=2 kb of read left beyond the
edge, and where those reads' other alignments go. Run on a haplotype's filtered PAF. Usage: gap_edges.py in.agp ont.paf.gz out.tsv'''
import sys, gzip, collections
agp, paf, out = sys.argv[1:4]
G = collections.defaultdict(list); rows = []
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] in "NU": G[s[0]].append((int(s[1]) - 1, int(s[2])))
slack = 1000; edge = collections.defaultdict(lambda: [0, 0]); A = collections.defaultdict(list)
for l in gzip.open(paf, "rt"):
    s = l.split("\t", 12)
    A[s[0]].append((s[5], int(s[7]), int(s[8]), s[4], int(s[2]), int(s[3]), int(s[1]), int(s[11])))
go = collections.defaultdict(collections.Counter)
for q, xs in A.items():
    for x in xs:
        t, ts, te, st, qs, qe, ql, mq = x
        if mq < 20 or t not in G: continue
        for g0, g1 in G[t]:
            for side in ("L", "R"):
                if side == "L" and g0 - slack <= te <= g0 and ts <= g0 - 5000: left = (ql - qe) if st == "+" else qs
                elif side == "R" and g1 <= ts <= g1 + slack and te >= g1 + 5000: left = qs if st == "+" else (ql - qe)
                else: continue
                k = (t, g0, side); edge[k][0] += 1
                if left >= 2000:
                    edge[k][1] += 1
                    for y in xs:
                        if y is not x and y[7] >= 1: go[k][y[0] if y[0] != t else ("same:" + str(round((y[1] - g0) / 1e6, 1)))] += 1
with open(out, "w") as o:
    o.write("#scaffold\tgap_start\tside\treach\tbeyond\twhere\n")
    for t in sorted(G):
        for g0, g1 in G[t]:
            for side in "LR":
                k = (t, g0, side); e = edge[k]
                o.write(f"{t}\t{g0}\t{side}\t{e[0]}\t{e[1]}\t{' '.join(f'{a}:{b}' for a, b in go[k].most_common(5))}\n")
v = [edge[k] for k in edge]
print("edges with reads", len(v), "reach sum", sum(x[0] for x in v), "beyond sum", sum(x[1] for x in v))
