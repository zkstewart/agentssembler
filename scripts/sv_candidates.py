#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Audit a finished assembly against the reference (never used to build it). Unitig-to-reference alignments are lifted
through the AGP, merged along each scaffold into blocks, and blocks >= --minLen on a chromosome other than the
scaffold's name, or on its own chromosome but opposite to the scaffold's main direction, are reported as SV candidates.
Usage: sv_candidates.py final.agp unitigs_vs_ref.paf.gz out.tsv [--minLen 500000 --minAln 20000]
'''
import sys, gzip, argparse, collections
p = argparse.ArgumentParser(); p.add_argument("agp"); p.add_argument("paf"); p.add_argument("out")
p.add_argument("--minLen", type=int, default=500000); p.add_argument("--minAln", type=int, default=20000); a = p.parse_args()
U = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": U[s[5]].append((s[0], int(s[1]) - 1, int(s[6]) - 1, int(s[7]), s[8]))
A = collections.defaultdict(list)
with gzip.open(a.paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 12)
        if "tp:A:P" not in s[12] or int(s[11]) < 20 or s[0] not in U: continue
        q0, q1 = int(s[2]), int(s[3])
        if q1 - q0 < a.minAln: continue
        for c, cs, u0, u1, st in U[s[0]]:
            x0, x1 = max(q0, u0), min(q1, u1)
            if x1 - x0 < a.minAln: continue
            p0 = cs + x0 - u0 if st == "+" else cs + u1 - x1
            strand = s[4] if st == "+" else ("-" if s[4] == "+" else "+")
            A[c].append((p0, p0 + x1 - x0, s[5], strand, int(s[7]), int(s[8])))
with open(a.out, "w") as o:
    o.write("#scaffold\tstart\tend\tref\tref_start\tref_end\tstrand\tkind\taligned_bp\n")
    for c, xs in sorted(A.items()):
        xs.sort(); chrom = c.split("_")[0]
        by = collections.Counter()
        for x in xs:
            if x[2] == chrom: by[x[3]] += x[1] - x[0]
        main = by.most_common(1)[0][0] if by else "+"
        blocks = []
        for x in xs:
            if blocks and blocks[-1][2] == x[2] and blocks[-1][3] == x[3]:
                b = blocks[-1]; b[1] = max(b[1], x[1]); b[4] = min(b[4], x[4]); b[5] = max(b[5], x[5]); b[6] += x[1] - x[0]
            else: blocks.append([x[0], x[1], x[2], x[3], x[4], x[5], x[1] - x[0]])
        # short interruptions by repeats should not split blocks: merge across blocks < minAln*5
        merged = []
        for b in blocks:
            if b[6] < a.minAln * 5 and merged: continue
            if merged and merged[-1][2] == b[2] and merged[-1][3] == b[3]:
                m = merged[-1]; m[1] = max(m[1], b[1]); m[4] = min(m[4], b[4]); m[5] = max(m[5], b[5]); m[6] += b[6]
            else: merged.append(b)
        for b in merged:
            if b[6] < a.minLen: continue
            k = "other_chromosome" if b[2] != chrom else ("inverted" if b[3] != main else None)
            if k: o.write(f"{c}\t{b[0]}\t{b[1]}\t{b[2]}\t{b[4]}\t{b[5]}\t{b[3]}\t{k}\t{b[6]}\n")
