#!/usr/bin/env python3
'''
Agreement between the reference-guided Hi-C phasing and HapHiC's reference-free grouping.

HapHiC clusters unitigs into chromosome x haplotype groups by Hi-C alone. For each
(chromosome, haplotype) group of the phasing, this reports the share of its bp that falls in
its single most-shared HapHiC group, and the bp-weighted overall share. High agreement from
two independent methods supports the phasing. Low agreement points to where one of them is
wrong (HapHiC often merges chromosomes, so check the contact maps).

Usage: haphic_agreement.py <phasing.tsv> <haphic scaffolds.agp>
'''
import sys
from collections import defaultdict
from common import natural_key

phasing, agp = sys.argv[1], sys.argv[2]
lab = {}
for l in open(phasing):
    if l.startswith("#"): continue
    s = l.rstrip("\n").split("\t")
    if s[3] != "assigned": continue
    regs, bps = s[6].split(","), [ int(x) for x in s[7].split(",") ]
    lab[s[0]] = f"{regs[bps.index(max(bps))].split(':')[0]}h{s[2]}"
t = defaultdict(lambda: defaultdict(int))
for l in open(agp):
    s = l.rstrip("\n").split("\t")
    if l.startswith("#") or s[4] != "W": continue
    u = s[5].split(":")[0].split("_")[0] # HapHiC-corrected fragments keep the unitig name as prefix
    if u in lab:
        t[lab[u]][s[0] if s[0].startswith("group") else "unplaced"] += int(s[7]) - int(s[6]) + 1
print("#phased_group\tbp\tdominant_haphic_group\tshare")
tot = agree = 0
for k in sorted(t, key=natural_key):
    T = sum(t[k].values()); g, b = max(t[k].items(), key=lambda x: x[1]); tot += T; agree += b
    print(f"{k}\t{T}\t{g}\t{b/T:.3f}")
print(f"OVERALL\t{tot}\t.\t{agree/tot:.3f}")
