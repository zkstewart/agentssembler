#!/usr/bin/env python3
'''Depth outside graph overlaps for named chains. Usage: hub_depth.py chains.agp gfa paf chain ...'''
import sys, gzip, collections
agp, gfa, paf = sys.argv[1:4]; Q = sys.argv[4:]
fl = {"+": "-", "-": "+"}; parts = collections.defaultdict(list); L = {}
for l in open(agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4] == "W": parts[s[0]].append((s[5], s[8])); L[s[0]] = int(s[2])
ov = collections.Counter()
for l in open(gfa):
    if l[0] == "L":
        s = l.rstrip("\n").split("\t"); o = int(s[5].rstrip("M")) if s[5] != "*" else 0
        ov[(s[1], s[2])] = max(ov[(s[1], s[2])], o); k = (s[3], fl[s[4]]); ov[k] = max(ov[k], o)
iv = {}
for c in Q:
    u, o = parts[c][-1]; r = ov[(u, o)]; u, o = parts[c][0]; a = ov[(u, fl[o])]; iv[c] = (a, L[c] - r)
cov = collections.Counter(); q = collections.Counter()
for l in gzip.open(paf, "rt"):
    s = l.split("\t", 12); c = s[5]
    if c in iv and "tp:A:P" in l:
        lo, hi = max(iv[c][0], int(s[7])), min(iv[c][1], int(s[8]))
        if hi > lo:
            cov[c] += hi - lo
            if int(s[11]) >= 20: q[c] += hi - lo
for c in Q:
    w = iv[c][1] - iv[c][0]
    print(c, L[c], "interior", w, "depth", round(cov[c] / w, 1) if w > 0 else "NA", "q20", round(q[c] / w, 1) if w > 0 else "NA")
