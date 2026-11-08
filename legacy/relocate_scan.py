#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Find unitigs whose Hi-C links point to a different chromosome than their reference
placement. Links are tallied to each (chromosome, haplotype) group of placed unitigs, and
normalised by group size (total unique Hi-C reads of the group's unitigs, as a proxy for
how many links a random unitig would form with it). A unitig is flagged when its best
group is on another chromosome and beats the best group on its own chromosome by --ratio.
'''
import sys, argparse
from collections import defaultdict
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("phasing"); p.add_argument("links")
p.add_argument("--ratio", type=float, default=2.0)
p.add_argument("--minLinks", type=int, default=100)
p.add_argument("--paf", default=None, help="unitig-to-reference PAF, to score alignment reliability")
p.add_argument("--reliableQ", type=int, default=20)
args = p.parse_args()

lab, chrom, length, status = {}, {}, {}, {}
for line in open(args.phasing):
    if line.startswith("#"): continue
    sl = line.rstrip("\n").split("\t")
    length[sl[0]] = int(sl[1]); status[sl[0]] = sl[3]
    if sl[6] == ".": continue
    regs, bps = sl[6].split(","), [int(x) for x in sl[7].split(",")]
    chrom[sl[0]] = regs[bps.index(max(bps))].split(":")[0]
    if sl[3] == "assigned": lab[sl[0]] = (chrom[sl[0]], sl[2])
# Alignment reliability: fraction of the unitig covered by confident (MAPQ >= reliableQ,
# >= 10 kb) primary alignments to its reference chromosome
reliable = {}
if args.paf:
    blocks = defaultdict(list)
    for line in open(args.paf):
        sl = line.split("\t")
        if "tp:A:S" in line or int(sl[11]) < args.reliableQ or int(sl[10]) < 10000: continue
        if sl[0] in chrom and sl[5] == chrom[sl[0]]:
            blocks[sl[0]].append((int(sl[2]), int(sl[3])))
    for q, b in blocks.items():
        b.sort(); covered, end = 0, 0
        for s0, e0 in b:
            if e0 > end:
                covered += e0 - max(s0, end); end = e0
        reliable[q] = covered / length[q]

depth, raw = {}, defaultdict(lambda: defaultdict(int))
for line in open(args.links):
    a, b, n = line.rstrip("\n").split("\t")
    if b == ".": depth[a] = int(n); continue
    for x, y in ((a, b), (b, a)):
        if y in lab: raw[x][lab[y]] += int(n)
groupDepth = defaultdict(int)
for q, g in lab.items(): groupDepth[g] += depth.get(q, 0)
mean = sum(groupDepth.values()) / len(groupDepth)

print("#contig\tstatus\tlength\tref_chrom\tlabel\tlinks\tbest_group\tbest_enrich\tbest_own_chrom_group\town_enrich\tratio\taln_reliable_frac")
for q in raw:
    if not q in chrom: continue
    tot = sum(raw[q].values())
    if tot < args.minLinks: continue
    # Exclude the unitig's own contribution to its group
    enr = {}
    for g, n in raw[q].items():
        gd = groupDepth[g] - (depth.get(q, 0) if lab.get(q) == g else 0)
        enr[g] = (n / tot) / (gd / mean / len(groupDepth)) if gd > 0 else 0
    best = max(enr, key=enr.get)
    own = [ g for g in enr if g[0] == chrom[q] ]
    bestOwn = max(own, key=enr.get) if own else None
    ownE = enr[bestOwn] if bestOwn else 0.0
    if best[0] != chrom[q] and enr[best] >= args.ratio * ownE:
        print(f"{q}\t{status[q]}\t{length[q]}\t{chrom[q]}\t{'%sh%s' % lab[q] if q in lab else '.'}\t{tot}\t{best[0]}h{best[1]}\t{enr[best]:.2f}\t"
              f"{'%sh%s' % bestOwn if bestOwn else '.'}\t{ownE:.2f}\t{enr[best]/max(ownE,1e-9):.1f}\t{reliable.get(q, 0.0):.2f}")
