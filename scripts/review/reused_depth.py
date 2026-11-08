#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Are arms used by both haplotypes truly 2x? For every chain in a path used by both groups of a chromosome, measures
depth outside the graph overlaps at both ends, and lists alternative paths for the junctions that used it.
Usage: reused_depth.py arm_report.tsv choices_base.tsv groups.tsv chroms.tsv chains.agp unitigs.gfa ont_vs_chains.paf.gz out.tsv'''
import sys, gzip, collections
rep, ch, grp, chroms, agp, gfa, paf, out = sys.argv[1:9]
g2c = {}
for l in open(chroms):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); g2c[s[0]] = s[1]
use = collections.defaultdict(set); how = collections.defaultdict(set); junc = collections.defaultdict(list)
for l in open(rep):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    for z in s[4].split(","):
        c = z[:-1]; use[c].add(s[0]); how[c].add(s[3]); junc[c].append(f"{s[0]}:{s[1]}>{s[2]}")
both = {c: g for c, g in use.items() if len({g2c.get(x) for x in g}) < len(g)}
parts = collections.defaultdict(list); L = {}
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] != "W": continue
    parts[s[0]].append((s[5], s[8])); L[s[0]] = max(L.get(s[0], 0), int(s[2]))
ov = collections.defaultdict(int)     # (unitig, strand leaving) -> largest overlap
for l in open(gfa):
    if l[0] != "L": continue
    s = l.rstrip("\n").split("\t"); o = int(s[5].rstrip("M")) if s[5] != "*" else 0
    k = (s[1], s[2]); ov[k] = max(ov[k], o)
    k = (s[3], {"+": "-", "-": "+"}[s[4]]); ov[k] = max(ov[k], o)
iv = {}
for c in both:
    u, o = parts[c][-1]; right = ov[(u, o)]
    u, o = parts[c][0]; left = ov[(u, {"+": "-", "-": "+"}[o])]
    iv[c] = (left, L[c] - right)
cov = {c: collections.Counter() for c in both}; n = collections.Counter(); nq = collections.Counter(); span = collections.Counter()
opn = gzip.open if paf.endswith(".gz") else open
for l in opn(paf, "rt"):
    s = l.split("\t", 12)
    c = s[5]
    if c not in iv: continue
    if "tp:A:P" not in l: continue
    a, b = iv[c]; x, y = int(s[7]), int(s[8]); lo, hi = max(a, x), min(b, y)
    if b > a and hi > lo:
        n[c] += hi - lo
        if int(s[11]) >= 20: nq[c] += hi - lo
    if x <= a and y >= b: span[c] += 1
with open(out, "w") as o:
    o.write("#chain\tlength\tinterior\tdepth\tdepth_q20\tspanning\tgroups\thow\tjunctions\n")
    tot = collections.Counter(); bp = collections.Counter()
    for c in sorted(both, key=lambda c: -(iv[c][1] - iv[c][0])):
        a, b = iv[c]; w = b - a
        d = n[c] / w if w > 0 else float("nan"); dq = nq[c] / w if w > 0 else float("nan")
        o.write(f"{c}\t{L[c]}\t{w}\t{d:.1f}\t{dq:.1f}\t{span[c]}\t{','.join(sorted(both[c]))}\t{','.join(sorted(how[c]))}\t{';'.join(junc[c])}\n")
        kind = "no_interior" if w <= 0 else "1x(<40)" if d < 40 else "2x(>=40)"
        tot[kind] += 1; bp[kind] += max(w, 0)
print(len(both), "chains"); print({k: (tot[k], bp[k]) for k in tot})
