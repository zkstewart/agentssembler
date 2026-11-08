#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Hi-C test of an inverted block [x, y] on a scaffold. Each flank is compared with the near and far ends of the block:
if the layout is right the near end has more contacts; if the block is really the other way round, the far end does.
Output: left flank vs block near / far, right flank vs block near / far.
Usage: inversion_test.py agp pairs.tsv.gz scaffold:x:y ... [--win 500000]
'''
import argparse, collections, gzip
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("agp"); p.add_argument("pairs"); p.add_argument("blocks", nargs="+"); p.add_argument("--win", type=int, default=500000)
a = p.parse_args()
W = a.win; U = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": U[s[5]].append((int(s[6]), int(s[7]), s[0], int(s[1]), s[8]))
def lift(u, x):
    for u0, u1, c, cs, o in U.get(u, ()):
        if u0 <= x <= u1: return c, (cs + x - u0 if o == "+" else cs + u1 - x)
B = [(b.split(":")[0], int(b.split(":")[1]), int(b.split(":")[2])) for b in a.blocks]
want = {c for c, *_ in B}; P = collections.defaultdict(list)
with gzip.open(a.pairs, "rt") as fh:
    for l in fh:
        s = l.split("\t")
        if s[0] not in U or s[2] not in U: continue
        A = lift(s[0], int(s[1])); C = lift(s[2], int(s[3]))
        if A and C and A[0] == C[0] and A[0] in want: P[A[0]].append((A[1], C[1]))
def n(c, r1, r2):
    return sum(1 for u, v in P[c] if (r1[0] <= u < r1[1] and r2[0] <= v < r2[1]) or (r2[0] <= u < r2[1] and r1[0] <= v < r1[1]))
print("#block\tleft_vs_near\tleft_vs_far\tright_vs_near\tright_vs_far")
for c, x, y in B:
    L, R = (x - W, x), (y, y + W); ne, fa = (x, x + W), (y - W, y)
    print(f"{c}:{x/1e6:.2f}-{y/1e6:.2f}\t{n(c,L,ne)}\t{n(c,L,fa)}\t{n(c,R,fa)}\t{n(c,R,ne)}")
