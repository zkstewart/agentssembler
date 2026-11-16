#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Hi-C tie-break for a pairing TTT cannot decide: which haplotype anchor does each exit chain contact more?
Counts pairs between each candidate chain and the part of each anchor chain within --window bp of the anchor's
tangle-facing end (anchor given as chainNNNN+ = 3' end faces the tangle, chainNNNN- = 5' end).
Usage: hic_phase.py pairs.tsv.gz chains.fai anchorA,anchorB cand1,cand2,... [--window 5000000]
Prints counts and, per candidate, the log2 ratio A/B. Same-haplotype contacts should exceed cross-haplotype ones
only modestly (homologs are collinear), so a ratio near 0 means no call.'''
import argparse, gzip, collections, math
p = argparse.ArgumentParser(); p.add_argument("pairs"); p.add_argument("fai"); p.add_argument("anchors"); p.add_argument("cands")
p.add_argument("--window", type=int, default=5000000); a = p.parse_args()
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
anc = [(x[:-1], x[-1]) for x in a.anchors.split(",")]; cands = a.cands.split(",")
def inwin(c, pos):
    for k, (n, e) in enumerate(anc):
        if c == n and ((e == "-" and pos <= a.window) or (e == "+" and pos >= L[n] - a.window)): return k
    return None
N = collections.Counter(); C = set(cands)
for l in gzip.open(a.pairs, "rt"):
    c1, p1, c2, p2 = l.split()[:4]
    for x, y, py in ((c1, c2, p2), (c2, c1, p1)):
        if x in C:
            k = inwin(y, int(py))
            if k is not None: N[(x, k)] += 1
print("candidate\t" + "\t".join(f"{n}{e}" for n, e in anc) + "\tlog2(A/B)")
for c in cands:
    A, B = N[(c, 0)], N[(c, 1)]
    print(f"{c}\t{A}\t{B}\t{math.log2((A + 1) / (B + 1)):+.2f}")
