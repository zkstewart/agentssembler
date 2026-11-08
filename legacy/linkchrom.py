#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# For given unitigs, report the share of their inter-unitig Hi-C links going to unitigs
# placed on each reference chromosome (by the phasing TSV), and to each (chrom, hap).
import sys
from collections import defaultdict
phasing, links = sys.argv[1], sys.argv[2]
targets = sys.argv[3:]
lab, chrom = {}, {}
for line in open(phasing):
    if line.startswith("#"): continue
    sl = line.rstrip("\n").split("\t")
    if sl[6] == ".": continue
    regs, bps = sl[6].split(","), [int(x) for x in sl[7].split(",")]
    chrom[sl[0]] = regs[bps.index(max(bps))].split(":")[0]
    if sl[3] == "assigned": lab[sl[0]] = f"{chrom[sl[0]]}h{sl[2]}"
share = { t: defaultdict(int) for t in targets }; shareL = { t: defaultdict(int) for t in targets }
for line in open(links):
    a, b, n = line.rstrip("\n").split("\t")
    if b == ".": continue
    for x, y in ((a, b), (b, a)):
        if x in share and y in chrom:
            share[x][chrom[y]] += int(n)
            if y in lab: shareL[x][lab[y]] += int(n)
for t in targets:
    tot = sum(share[t].values()) or 1
    top = sorted(share[t].items(), key=lambda x: -x[1])[:4]
    topL = sorted(shareL[t].items(), key=lambda x: -x[1])[:4]
    print(f"{t}\t{chrom.get(t)}\tn={tot}\t" + " ".join(f"{c}:{v/tot:.0%}" for c, v in top) + "\t| " + " ".join(f"{c}:{v/tot:.0%}" for c, v in topL))
