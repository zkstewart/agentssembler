#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Where does a small chain belong in its group by Hi-C alone? For each query chain prints contact density with each ordered
chain (pairs per Mb of target) with order rank, and the end-to-end profile: pairs per Mb to the flanking 100 kb of each target end.
Usage: place_test.py order.tsv pairs.tsv.gz chains.fasta.fai query1 query2 ...'''
import sys, gzip, collections
od, pa, fai = sys.argv[1:4]; Q = set(sys.argv[4:])
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
O = collections.defaultdict(list); gof = {}
for l in open(od):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); O[s[0]].append(s[2]); gof[s[2]] = s[0]
C = {q: collections.Counter() for q in Q}
for l in gzip.open(pa, "rt"):
    s = l.split("\t")
    a, b = s[0], s[2]
    if a in C: C[a][b] += 1
    if b in C: C[b][a] += 1
for q in sorted(Q):
    g = gof[q]; order = O[g]; i = order.index(q)
    print(f"\n{q} {g} rank {i+1}/{len(order)} len {L[q]/1e3:.0f} kb; total to group {sum(C[q][c] for c in order if c != q)}; to other {sum(v for c, v in C[q].items() if c not in order)}")
    rows = [(C[q][c] / (L[c] / 1e6), c, k) for k, c in enumerate(order) if c != q]
    top = sorted(rows, reverse=True)[:6]
    print("  top density: " + " ".join(f"{c}@{k+1}:{d:.0f}" for d, c, k in top))
    print("  neighbours: " + " ".join(f"{order[k]}@{k+1}:{C[q][order[k]]/(L[order[k]]/1e6):.0f}" for k in range(max(0, i-2), min(len(order), i+3)) if k != i))
