#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# For each scaffold in a HapHiC AGP, measure how well its unitig order agrees with the
# reference: the length-weighted longest increasing subsequence of reference positions
# (taking the better of forward/reverse), over unitigs placed by the reference.
import sys
from collections import defaultdict
from common import natural_key
agp, support = sys.argv[1], sys.argv[2]
pos = {}
for line in open(support):
    if line.startswith("#"): continue
    sl = line.rstrip("\n").split("\t")
    pos[sl[0]] = (sl[1], (int(sl[2]) + int(sl[3])) / 2)
order = defaultdict(list)
for line in open(agp):
    sl = line.rstrip("\n").split("\t")
    if sl[4] != "W" or sl[0].startswith("utg"): continue # skip HapHiC's unplaced unitigs
    q = sl[5].split(":")[0]
    order[sl[0]].append((q, int(sl[7]) - int(sl[6]) + 1))
def wlis(items):
    best = []
    for i, (x, w) in enumerate(items):
        best.append(w + max([ best[j] for j in range(i) if items[j][0] < x ], default=0))
    return max(best, default=0)
for scaf in sorted(order, key=natural_key):
    chromName = scaf.rsplit("h", 1)[0]
    items = [ (pos[q][1], L) for q, L in order[scaf] if q in pos and pos[q][0] == chromName ]
    tot = sum(L for _, L in items) or 1
    fwd = wlis(items); rev = wlis([ (-x, w) for x, w in items ])
    allLen = sum(L for _, L in order[scaf])
    print(f"{scaf}\t{allLen/1e6:.1f}Mb\tref-placed {tot/1e6:.1f}Mb\tcollinear {max(fwd, rev)/tot:.1%}")
