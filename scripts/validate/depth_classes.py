#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Depth of the final scaffolds in 50 kb windows, split by whether the window lies in sequence placed in both haplotypes, from each
haplotype's masked mapping, plus the share of each chromosome that is shared. Coordinates are the pre-fill ones, as for the PAFs.
Usage: depth_classes.py scaffolds.agp target.masked.tsv a.paf.gz [b.paf.gz ...] out.tsv'''
import sys, gzip, collections
import numpy as np
agp, mask, *pafs, out = sys.argv[1:]; W = 50000
L = {}; N = collections.defaultdict(list)
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); L[s[0]] = max(L.get(s[0], 0), int(s[2]))
    if s[4] in "NU": N[s[0]].append((int(s[1]) - 1, int(s[2])))
sh = {c: np.zeros(n // W + 2) for c, n in L.items()}; bp = collections.Counter()
for l in open(mask):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); c, a, b = s[0], int(s[1]) - 1, int(s[2]); bp[c] += b - a
    for w in range(a // W, b // W + 1):
        sh[c][w] += max(0, min(b, (w + 1) * W) - max(a, w * W))
D = {c: np.zeros(n // W + 2) for c, n in L.items()}; Q = {c: np.zeros(n // W + 2) for c, n in L.items()}
for p in pafs:
    for l in gzip.open(p, "rt"):
        s = l.split("\t", 12)
        if s[5] not in D or "tp:A:P" not in l: continue
        a, b = int(s[7]), int(s[8]); q = int(s[11]) >= 20
        for w in range(a // W, b // W + 1):
            o = max(0, min(b, (w + 1) * W) - max(a, w * W)); D[s[5]][w] += o
            if q: Q[s[5]][w] += o
rows = []
for c in sorted(L):
    for w in range(L[c] // W):
        a = w * W; b = a + W
        if any(x < b and y > a for x, y in N[c]): continue
        f = sh[c][w] / W; rows.append((c, a, D[c][w] / W, Q[c][w] / W, f))
with open(out, "w") as o:
    o.write("#scaffold\tstart\tdepth\tdepth_q20\tshared_fraction\n")
    for r in rows: o.write(f"{r[0]}\t{r[1]}\t{r[2]:.1f}\t{r[3]:.1f}\t{r[4]:.2f}\n")
def summ(name, v):
    v = np.array(v); print(f"{name}: {len(v)} windows, median {np.median(v):.1f}, q1 {np.quantile(v,.25):.1f}, q3 {np.quantile(v,.75):.1f}, <10x {np.mean(v<10):.4f}, <20x {np.mean(v<20):.4f}, >1.5x median {np.mean(v>1.5*np.median(v)):.4f}, CV {v.std()/v.mean():.3f}")
summ("specific (<10% shared)", [r[2] for r in rows if r[4] < .1])
summ("shared (>90% shared)", [r[2] for r in rows if r[4] > .9])
summ("mixed", [r[2] for r in rows if .1 <= r[4] <= .9])
print("share of each scaffold placed in both haplotypes:")
for c in sorted(L, key=lambda x: (int(x[3:].split('_')[0]), x)): print(f"  {c}\t{L[c]/1e6:.2f} Mb\t{bp[c]/1e6:.2f} Mb\t{bp[c]/L[c]:.3f}")
