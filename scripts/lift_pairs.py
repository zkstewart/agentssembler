#!/usr/bin/env python3
'''
Lift chain-level Hi-C pairs, in which shared chains have one copy per haplotype, to final scaffold coordinates. A pair on a chain copy
goes to that copy's own scaffold, so a unitig that occurs in both haplotypes is not ambiguous. Positions outside the AGP spans (trimmed
overlaps) are dropped. Output columns: scaffold1 pos1 scaffold1 pos1 ... as scaffoldA posA scaffoldB posB, plus a second file with only pairs where
neither end lies on a chain copy or a unitig placed in two scaffolds: those measure phasing, since shared pairs were given to both copies by design.

Usage: lift_pairs.py pairs.tsv.gz chains.agp final.agp groups.tsv names.tsv out_all.tsv.gz out_specific.tsv.gz
'''
import sys, gzip, bisect, collections
pa, cagp, fagp, gr, nm, oa, os_ = sys.argv[1:8]
names = {}
for l in open(nm):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); names[s[0]] = s[1]
group = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); group[s[0]] = s[2]
C = collections.defaultdict(list)
for l in open(cagp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": C[s[0]].append((int(s[1]), int(s[2]), s[5], int(s[6]), int(s[7]), s[8]))
F = collections.defaultdict(lambda: collections.defaultdict(list)); where = collections.defaultdict(set)
for l in open(fagp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W":
        F[s[5]][s[0]].append((int(s[6]), int(s[7]), int(s[1]), s[8])); where[s[5]].add(s[0])
def to_unitig(c, p):
    for a, b, u, us, ue, o in C[c]:
        if a <= p <= b: return (u, us + p - a if o == "+" else ue - (p - a))
    return None
def lift(c, p):
    r = to_unitig(c, p)
    if r is None: return None
    u, x = r; g = group.get(c)
    sc = names.get(g)
    if sc is None:                         # an arm or other chain with no group: it has one scaffold or none
        if len(where[u]) != 1: return None
        sc = next(iter(where[u]))
    for us, ue, st, o in F[u].get(sc, ()):
        if us <= x <= ue: return (sc, st + x - us if o == "+" else st + ue - x, len(where[u]) > 1)
    return None
n = collections.Counter()
with gzip.open(oa, "wt", compresslevel=1) as A, gzip.open(os_, "wt", compresslevel=1) as S:
    for l in gzip.open(pa, "rt"):
        s = l.rstrip("\n").split("\t")
        x = lift(s[0], int(s[1])); y = lift(s[2], int(s[3])); n["in"] += 1
        if x is None or y is None: n["dropped"] += 1; continue
        A.write(f"{x[0]}\t{x[1]}\t{y[0]}\t{y[1]}\n"); n["lifted"] += 1
        if not (x[2] or y[2] or ".g" in s[0] or ".g" in s[2]): S.write(f"{x[0]}\t{x[1]}\t{y[0]}\t{y[1]}\n"); n["specific"] += 1
print(dict(n))
