#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Where do the two flanks of a breakpoint make Hi-C contact? For a flank (--skip to --win bp from the breakpoint), counts
pairs with one end in it, by destination: same flank, the opposite flank, and elsewhere by scaffold.
A true junction: flanks well covered and in contact with each other. A chimeric one: flanks well covered, little contact
across, and each flank's long-range contacts going to a different place. Flanks with few pairs in total say nothing.
Also reports read depth around the breakpoint (all alignments, and MAPQ >= 20).

Usage: sv_flanks.py agp ont.paf.gz pairs.tsv.gz scaffold:pos [scaffold:pos ...]
'''
import argparse, collections, gzip
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("agp"); p.add_argument("paf"); p.add_argument("pairs"); p.add_argument("bps", nargs="+")
p.add_argument("--skip", type=int, default=100000); p.add_argument("--win", type=int, default=1000000)
a = p.parse_args()
U = collections.defaultdict(list); slen = {}
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); slen[s[0]] = int(s[2])
    if s[4] == "W": U[s[5]].append((int(s[6]), int(s[7]), s[0], int(s[1]), s[8]))
def lift(u, x):
    for u0, u1, c, cs, o in U.get(u, ()):
        if u0 <= x <= u1: return c, (cs + x - u0 if o == "+" else cs + u1 - x)
    return None
B = [(b.split(":")[0], int(b.split(":")[1])) for b in a.bps]
byc = collections.defaultdict(list)
for i, (c, x) in enumerate(B): byc[c].append(i)
def side(i, y):
    x = B[i][1]; d = y - x
    if -a.win <= d <= -a.skip: return "L"
    if a.skip <= d <= a.win: return "R"
    return None
C = [collections.Counter() for _ in B]
with gzip.open(a.pairs, "rt") as fh:
    for l in fh:
        s = l.split("\t")
        if s[0] not in U or s[2] not in U: continue
        P = lift(s[0], int(s[1])); Q = lift(s[2], int(s[3]))
        if not P or not Q: continue
        for X, Y in ((P, Q), (Q, P)):
            for i in byc.get(X[0], ()):
                f = side(i, X[1])
                if not f: continue
                if Y[0] == X[0]:
                    g = side(i, Y[1])
                    if g == f: C[i][f + "self"] += 1
                    elif g: C[i]["cross"] += 1
                    else: C[i][f + ":" + ("same_left" if Y[1] < B[i][1] else "same_right")] += 1
                else: C[i][f + "@" + Y[0]] += 1
dep = [[0, 0, 0] for _ in B]
with gzip.open(a.paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 12)
        if s[5] not in byc or "tp:A:P" not in s[12]: continue
        t0, t1, q = int(s[7]), int(s[8]), int(s[11])
        for i in byc[s[5]]:
            x = B[i][1]
            if t0 <= x <= t1:
                dep[i][0] += 1
                if q >= 20: dep[i][1] += 1
                if t0 <= x - 5000 and t1 >= x + 5000 and q >= 20: dep[i][2] += 1
print("#bp\tdepth_any\tdepth_q20\tspan5k_q20\tL_self\tR_self\tcross(both counted)\tL_top_other\tR_top_other")
for i, (c, x) in enumerate(B):
    n = C[i]
    top = lambda f: ",".join(f"{k[2:]}:{v}" for k, v in sorted(((k, v) for k, v in n.items() if k.startswith(f + "@")), key=lambda t: -t[1])[:3])
    print(f"{c}:{x}\t{dep[i][0]}\t{dep[i][1]}\t{dep[i][2]}\t{n['Lself']//2}\t{n['Rself']//2}\t{n['cross']//2}\t{top('L')}\t{top('R')}")
