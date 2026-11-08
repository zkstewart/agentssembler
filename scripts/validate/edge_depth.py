#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Depth against distance from every gap edge, any MAPQ, primary alignments, on one haplotype's mapping.
Usage: edge_depth.py in.agp ont.paf.gz'''
import sys, gzip, collections
agp, paf = sys.argv[1:3]
G = collections.defaultdict(list)
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] in "NU": G[s[0]].append((int(s[1]) - 1, int(s[2])))
bins = [(0, 200), (200, 1000), (1000, 3000), (3000, 6000), (6000, 12000), (12000, 25000), (25000, 50000), (50000, 100000)]
D = {}
for t, v in G.items():
    for g0, g1 in v:
        D[(t, g0, "L")] = [g0 - b[1] for b in bins[::-1]] and [0] * len(bins); D[(t, g0, "R")] = [0] * len(bins)
idx = collections.defaultdict(list)
for t, v in G.items():
    for g0, g1 in v: idx[t].append((g0, g1))
ends_near = collections.Counter(); covered = 0
for l in gzip.open(paf, "rt"):
    s = l.split("\t", 12); t = s[5]
    if t not in idx: continue
    x, y = int(s[7]), int(s[8])
    for g0, g1 in idx[t]:
        if x < g0 and y > g0 - 100000:
            for i, (a, b) in enumerate(bins):
                lo, hi = g0 - b, g0 - a; o = min(y, hi) - max(x, lo)
                if o > 0: D[(t, g0, "L")][i] += o
        if y > g1 and x < g1 + 100000:
            for i, (a, b) in enumerate(bins):
                lo, hi = g1 + a, g1 + b; o = min(y, hi) - max(x, lo)
                if o > 0: D[(t, g0, "R")][i] += o
for k, v in D.items():
    print(k[0], k[1], k[2], *[f"{c/(b[1]-b[0]):.1f}" for c, b in zip(v, bins)], sep="\t")
