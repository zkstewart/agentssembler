#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Depth of the final scaffolds in 50 kb windows, split by whether the window lies in sequence placed in both haplotypes, from each
haplotype's masked mapping, plus the share of each chromosome that is shared. Coordinates are the pre-fill ones, as for the PAFs.
With --chainsAgp (the chain AGP), also writes depth.2x_chains.tsv beside out.tsv: haplotype-specific windows (< 10% shared) above 1.5x the
specific median, each given the chain that covers most of it. Two or more on one chain is the stage 15 alarm: a chain placed on one haplotype
that holds both collapsed (a hub core, or a shared chain placed once), which the median depth that classed it could not show.
Usage: depth_classes.py scaffolds.agp target.masked.tsv a.paf.gz [b.paf.gz ...] out.tsv [--chainsAgp chains.agp]'''
import sys, gzip, collections, os, re
import numpy as np
argv = sys.argv[1:]; chainsAgp = None
if "--chainsAgp" in argv: i = argv.index("--chainsAgp"); chainsAgp = argv[i + 1]; del argv[i:i + 2]
agp, mask, *pafs, out = argv; W = 50000
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
if chainsAgp:
    U = {}
    for l in open(chainsAgp):
        s = l.split("\t")
        if l[0] != "#" and s[4] not in "NU": U[s[5]] = s[0]
    C = collections.defaultdict(list)
    for l in open(agp):
        s = l.rstrip("\n").split("\t")
        if l[0] != "#" and s[4] not in "NU": C[s[0]].append((int(s[1]) - 1, int(s[2]), U.get(s[5], s[5])))
    sp = [r for r in rows if r[4] < .1]; m = float(np.median([r[2] for r in sp])) if sp else 0
    with open(os.path.join(os.path.dirname(out) or ".", "depth.2x_chains.tsv"), "w") as o:
        o.write(f"#scaffold\tstart\tdepth\tx_specific_median\tchain\n")
        for c, a, d, q, f in sp:
            if d <= 1.5 * m: continue
            ov = collections.Counter()
            for x, y, ch in C[c]:
                if min(a + W, y) > max(a, x): ov[ch] += min(a + W, y) - max(a, x)
            o.write(f"{c}\t{a}\t{d:.1f}\t{d / m:.2f}\t{ov.most_common(1)[0][0] if ov else '.'}\n")
summ("shared (>90% shared)", [r[2] for r in rows if r[4] > .9])
summ("mixed", [r[2] for r in rows if .1 <= r[4] <= .9])
print("share of each scaffold placed in both haplotypes:")
# natural sort: scaffold names need not be chrN_hM
for c in sorted(L, key=lambda x: [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', x)]):  print(f"  {c}\t{L[c]/1e6:.2f} Mb\t{bp[c]/1e6:.2f} Mb\t{bp[c]/L[c]:.3f}")
