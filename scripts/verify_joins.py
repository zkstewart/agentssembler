#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Reads through each chain join. A join at chain position p is supported by a read whose single
alignment covers p - anchor to p + anchor. Spanning reads are counted at MAPQ >= 20 (unique) and any
MAPQ, along with the local depth 20-30 kb either side of the join (the expectation at that place).
Usage: verify_joins.py <reads.paf.gz> <joins.tsv> <out.tsv> [--anchor 2000]
'''
import gzip, bisect, argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("paf"); p.add_argument("joins"); p.add_argument("out"); p.add_argument("--anchor", type=int, default=2000)
a = p.parse_args()
J = collections.defaultdict(list)
for l in open(a.joins):
    if l[0] == "#": continue
    s = l.rstrip().split("\t"); J[s[0]].append([int(s[1])] + s[2:] + [0, 0, 0, 0])
for v in J.values(): v.sort(key=lambda x: x[0])
st = { k: [ x[0] for x in v ] for k, v in J.items() }
with gzip.open(a.paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 12); t = s[5]
        if t not in J: continue
        x, y, q = int(s[7]), int(s[8]), int(s[11]); P = "tp:A:P" in s[12]
        if not P: continue
        v = J[t]
        i = bisect.bisect_left(st[t], x - 40000)
        while i < len(v) and v[i][0] <= y + 40000:
            j = v[i]
            if x <= j[0] - a.anchor and y >= j[0] + a.anchor:
                j[4] += 1
                if q >= 20: j[5] += 1
            lo = j[0] - 30000; hi = j[0] - 20000
            if x < hi and y > lo: j[6] += (min(y, hi) - max(x, lo)) / 10000
            lo = j[0] + 20000; hi = j[0] + 30000
            if x < hi and y > lo: j[7] += (min(y, hi) - max(x, lo)) / 10000
            i += 1
with open(a.out, "w") as o:
    o.write("#chain\tpos\tleft\tright\toverlap\tspan_any\tspan_q20\tdepth_left\tdepth_right\n")
    for k, v in J.items():
        for j in v: o.write(f"{k}\t{j[0]}\t{j[1]}\t{j[2]}\t{j[3]}\t{j[4]}\t{j[5]}\t{j[6]:.1f}\t{j[7]:.1f}\n")
