#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Evidence at the breakpoints of SV candidates, from data other than the reference. Coordinates are those of the
scaffold the reads were aligned to (the AGP before read filling).

For each breakpoint (candidate edge away from a scaffold end):
  gap_kb       distance to the nearest AGP gap, or NA if none within --window
  min_reads    fewest reads (primary, MAPQ >= 20) spanning a position within +-window with >= --anchor bp each side.
               Above 0 means the sequence is continuous in reads across the whole window
  hic_cross    Hi-C pairs across the breakpoint, each end between --hicSkip and --hicWin from it
  hic_ratio    hic_cross / the scaffold's median over positions every 100 kb. A misjoin shows as a ratio near 0
  hic_pct      percentile of hic_cross among those positions

Usage: sv_evidence.py candidates.tsv agp ont.paf.gz pairs.tsv.gz out.tsv [--extra scaffold:pos ...]
'''
import argparse, bisect, collections, gzip
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("candidates"); p.add_argument("agp"); p.add_argument("paf"); p.add_argument("pairs"); p.add_argument("out")
p.add_argument("--window", type=int, default=50000); p.add_argument("--anchor", type=int, default=5000)
p.add_argument("--hicWin", type=int, default=200000); p.add_argument("--extra", nargs="*", default=[])
p.add_argument("--hicSkip", type=int, default=0, help="ignore this many bp each side of the breakpoint, for local satellite")
a = p.parse_args()

U = collections.defaultdict(list); gaps = collections.defaultdict(list); slen = {}
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); slen[s[0]] = int(s[2])
    if s[4] == "W": U[s[5]].append((int(s[6]), int(s[7]), s[0], int(s[1]), s[8]))
    else: gaps[s[0]].append((int(s[1]) + int(s[2])) // 2)

bps = []
for l in open(a.candidates):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); label = f"{s[7]}:{s[3]}:{int(s[1])/1e6:.2f}-{int(s[2])/1e6:.2f}"
    for x, side in ((int(s[1]), "start"), (int(s[2]), "end")):
        if 50000 < x < slen[s[0]] - 50000: bps.append((s[0], x, label, side))
for e in a.extra:
    c, x = e.split(":"); bps.append((c, int(x), "extra", "."))
want = {c for c, *_ in bps}

# Reads
spans = collections.defaultdict(list)
with gzip.open(a.paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 12)
        if s[5] in want and "tp:A:P" in s[12] and int(s[11]) >= 20: spans[s[5]].append((int(s[7]), int(s[8])))

# Hi-C, lifted through the AGP
for u in U: U[u].sort()
def lift(u, x):
    for u0, u1, c, cs, o in U.get(u, ()):
        if u0 <= x <= u1: return c, (cs + x - u0 if o == "+" else cs + u1 - x)
    return None
hic = collections.defaultdict(list)
with gzip.open(a.pairs, "rt") as fh:
    for l in fh:
        s = l.split("\t")
        if s[0] not in U or s[2] not in U: continue
        A = lift(s[0], int(s[1])); B = lift(s[2], int(s[3]))
        if A and B and A[0] == B[0] and A[0] in want: hic[A[0]].append((min(A[1], B[1]), max(A[1], B[1])))
for c in hic: hic[c].sort()
def cross(c, x):
    W = a.hicWin; S = a.hicSkip; i = bisect.bisect_left(hic[c], (x - W, 0)); j = bisect.bisect_right(hic[c], (x - S, 10**12)); n = 0
    for y0, y1 in hic[c][i:j]:
        if x + S < y1 <= x + W: n += 1
    return n
bg = {}
for c in want:
    v = sorted(cross(c, x) for x in range(a.hicWin, slen[c] - a.hicWin, 100000))
    bg[c] = v

with open(a.out, "w") as o:
    o.write("#scaffold\tpos\tcandidate\tedge\tgap_kb\tmin_reads\thic_cross\thic_ratio\thic_pct\n")
    for c, x, label, side in sorted(bps):
        lo, hi = max(0, x - a.window), min(slen[c], x + a.window); step = 1000
        cnt = [0] * ((hi - lo) // step + 1)
        for t0, t1 in spans[c]:
            if t1 - t0 < 2 * a.anchor or t1 <= lo or t0 >= hi: continue
            for i in range(max(0, (t0 + a.anchor - lo) // step + 1), min(len(cnt) - 1, (t1 - a.anchor - lo) // step) + 1):
                cnt[i] += 1
        near = [abs(g - x) for g in gaps[c] if abs(g - x) <= a.window]
        h = cross(c, x); v = bg[c]; med = v[len(v) // 2] if v else 0
        pct = 100 * bisect.bisect_left(v, h) / max(1, len(v))
        o.write(f"{c}\t{x}\t{label}\t{side}\t{min(near)/1000 if near else 'NA'}\t{min(cnt)}\t{h}\t{h/med if med else 0:.2f}\t{pct:.1f}\n")
