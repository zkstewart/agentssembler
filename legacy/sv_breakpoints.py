#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
For two pseudochromosomes suspected of a reciprocal translocation, score each unitig (in
AGP order) by its Hi-C links to two "core" regions known to lie on different derivative
chromosomes, to locate the breakpoints. Prints the unitig order with position, length,
links to each core, and the log2 ratio (A over B).
'''
import sys, math
from collections import defaultdict
agp, links, scafA, coreA, scafB, coreB = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
# coreA/coreB are 'start-end' fractions of the scaffold length, e.g. 0-0.4
order = defaultdict(list); length = {}
for l in open(agp):
    if l.startswith("#"): continue
    sl = l.rstrip("\n").split("\t")
    if sl[4] == "W":
        order[sl[0]].append((sl[5], int(sl[1]), int(sl[2]), sl[8]))
    length[sl[0]] = max(length.get(sl[0], 0), int(sl[2]))
def core(scaf, frac):
    a, b = (float(x) for x in frac.split("-"))
    L = length[scaf]
    return { u for u, s, e, _ in order[scaf] if s >= a * L and e <= b * L }
A, B = core(scafA, coreA), core(scafB, coreB)
cnt = defaultdict(lambda: [0, 0])
for l in open(links):
    a, b, n = l.rstrip("\n").split("\t")
    if b == ".": continue
    for x, y in ((a, b), (b, a)):
        if y in A: cnt[x][0] += int(n)
        if y in B: cnt[x][1] += int(n)
for scaf in (scafA, scafB):
    print(f"## {scaf} ({length[scaf]/1e6:.1f} Mb)")
    for u, s, e, o in order[scaf]:
        a, b = cnt[u]
        r = math.log2((a + 1) / (b + 1))
        tag = "A" if r > 1 else ("B" if r < -1 else ".")
        print(f"{scaf}\t{u}\t{s/1e6:6.2f}-{e/1e6:6.2f}Mb\t{(e-s+1)/1e6:.2f}Mb\t{o}\tA={a}\tB={b}\tlog2={r:+.1f}\t{tag}{' core' if u in A or u in B else ''}")
