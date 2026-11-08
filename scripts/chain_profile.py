#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Hi-C partners along a chain, in windows, to detect chimeric chains: for each window, pairs to each other
group (or to the chain's own remainder). Also lists each unitig's window.
Usage: chain_profile.py chain pairs.tsv.gz groups.tsv chains.agp [window]
'''
import sys, gzip, collections
ch, pairs, groups, agp = sys.argv[1:5]; W = int(sys.argv[5]) if len(sys.argv) > 5 else 500000
g = {}
for l in open(groups):
    if l[0] != "#":
        c, n, x = l.rstrip().split("\t"); g[c] = x
C = collections.defaultdict(collections.Counter)
with gzip.open(pairs, "rt") as fh:
    for l in fh:
        x, u, y, v = l.split("\t")
        if x == ch and y != ch: C[int(u) // W][g.get(y, "?")] += 1
        elif y == ch and x != ch: C[int(v) // W][g.get(x, "?")] += 1
        elif x == ch and y == ch and abs(int(u) - int(v)) > 200000:
            C[int(u) // W]["self"] += 1; C[int(v) // W]["self"] += 1
for w in sorted(C):
    t = sum(C[w].values())
    print(f"{w*W/1e6:7.1f} {t:6d} " + " ".join(f"{k}:{n}" for k, n in C[w].most_common(5)))
print("unitig joins:")
for l in open(agp):
    if l.startswith(ch + "\t"):
        s = l.rstrip().split("\t")
        if s[4] == "W": print(f"  {int(s[1])/1e6:7.2f}-{int(s[2])/1e6:7.2f} {s[5]}{s[8]}")
