#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Test every unitig junction of every chain with Hi-C: a true join is crossed by pairs as often as the rest
of the chain, and both sides talk to the same other chains.
  cross_ratio    pairs crossing the junction with both ends within --dist, relative to the chain's median
                 at other positions (insulation).
  partner_cos    cosine between left- and right-flank contact profiles over other chains (--flank bp), each
                 side normalised; a chimera has two different chromosomes' profiles.
Sides shorter than --minSide are not tested. Output sorted by cross_ratio.
Usage: chain_breaks.py pairs.tsv.gz chains.agp out.tsv [--dist 500000 --flank 1000000 --minSide 100000]
'''
import gzip, argparse, collections, bisect, math
import numpy as np
p = argparse.ArgumentParser(); p.add_argument("pairs"); p.add_argument("agp"); p.add_argument("out")
p.add_argument("--dist", type=int, default=500000); p.add_argument("--flank", type=int, default=1000000); p.add_argument("--minSide", type=int, default=100000)
a = p.parse_args()
J = collections.defaultdict(list); L = collections.Counter()
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip().split("\t")
    L[s[0]] = max(L[s[0]], int(s[2]))
    if s[4] == "W": J[s[0]].append((int(s[1]) - 1, int(s[2]), s[5] + s[8]))
W = collections.defaultdict(list); T = collections.defaultdict(list)
with gzip.open(a.pairs, "rt") as fh:
    for l in fh:
        x, u, y, v = l.split("\t"); u = int(u); v = int(v)
        if x == y: W[x].append((min(u, v), max(u, v)))
        else: T[x].append((u, y)); T[y].append((v, x))
def crossing(arr, pos):
    '''pairs with u<pos<v, pos-u<dist, v-pos<dist'''
    if len(arr) == 0: return 0
    m = (arr[:, 0] < pos) & (arr[:, 1] >= pos) & (arr[:, 0] >= pos - a.dist) & (arr[:, 1] < pos + a.dist)
    return int(m.sum())
out = []
for c, us in J.items():
    if len(us) < 2: continue
    arr = np.asarray(W[c]) if W[c] else np.zeros((0, 2))
    grid = np.arange(a.dist // 2, max(L[c] - a.dist // 2, a.dist // 2 + 1), 100000)
    bg = np.array([crossing(arr, g) for g in grid]) if len(grid) else np.array([0])
    med = max(float(np.median(bg)), 1.0)
    tr = T[c]
    for (s1, e1, u1), (s2, e2, u2) in zip(us, us[1:]):
        pos = e1
        if pos < a.minSide or L[c] - pos < a.minSide: continue
        lc = collections.Counter(); rc = collections.Counter()
        for q, y in tr:
            if pos - a.flank <= q < pos: lc[y] += 1
            elif pos <= q < pos + a.flank: rc[y] += 1
        ks = set(lc) | set(rc)
        dot = sum(lc[k] * rc[k] for k in ks); na = math.sqrt(sum(v * v for v in lc.values())) or 1; nb = math.sqrt(sum(v * v for v in rc.values())) or 1
        cr = crossing(arr, pos)
        out.append((cr / med, c, pos, u1, u2, cr, med, dot / na / nb, sum(lc.values()), sum(rc.values())))
out.sort()
with open(a.out, "w") as o:
    o.write("#cross_ratio\tchain\tpos\tleft_unitig\tright_unitig\tcross\tmedian\tpartner_cos\tleft_trans\tright_trans\n")
    for r in out: o.write(f"{r[0]:.3f}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.0f}\t{r[7]:.3f}\t{r[8]}\t{r[9]}\n")
print(len(out), "junctions tested")
