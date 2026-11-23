#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Stage 07, for D05: gene copies in the groups, from compleasm on the unitigs (stage 01). A gene is short when its complete copies in grouped
chains number under the ploidy (a shared chain counts once per haplotype, so ploidy times); it is recoverable when an unassigned
chain carries a copy. With --rescued (rescue_groups.py's table: the rescue's own, or stage 07's preview when hic_allele_rescue is
off), also counts the short genes the rescued chains would make whole. Counts are of compleasm genes on chains, a proxy:
fragments of one haplotype can hold the same gene twice, so read recoverable as an upper bound.
Usage: genes_placed.py --genes full_table_busco_format.tsv --agp chains.agp --groups groups.tsv --ploidy N --out chains.tsv [--rescued groups.rescued.tsv]'''
import argparse, collections, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chainlib import read_agp, gene_hits
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("genes", "agp", "groups", "out"): p.add_argument("--" + x, required=True)
p.add_argument("--ploidy", type=int, required=True); p.add_argument("--rescued")
a = p.parse_args()
grp = {}; L = {}
for l in open(a.groups):
    if l[0] != "#": s = l.rstrip("\n").split("\t"); grp[s[0]] = s[2]; L[s[0]] = int(s[1])
resc = {}
if a.rescued and os.path.exists(a.rescued):
    for l in open(a.rescued):
        if l[0] != "#":
            s = l.rstrip("\n").split("\t")
            if s[3] != "unassigned": resc[s[0]] = s[4]
P, _ = read_agp(a.agp)
def weight(c, extra=()):
    g = grp.get(c, "unassigned")
    if g == "shared": return a.ploidy
    if re.match(r"g\d+$", g) or c in extra: return 1
    return 0
placed = collections.Counter(); after = collections.Counter(); unas = collections.defaultdict(set); allg = set()
for g, st, c, x in gene_hits(a.genes, P):
    if st == "Fragmented": continue
    allg.add(g)
    if c is None: continue
    placed[g] += weight(c); after[g] += weight(c, resc)
    if grp.get(c, "unassigned") == "unassigned": unas[g].add(c)
short = [g for g in allg if placed[g] < a.ploidy]; rec = [g for g in short if unas[g]]
per = collections.Counter(c for g in rec for c in unas[g])
with open(a.out, "w") as o:
    o.write("#chain\tlength\tshort_genes\trescue\n")
    for c, n in sorted(per.items(), key=lambda t: (-t[1], t[0])): o.write(f"{c}\t{L.get(c, 0)}\t{n}\t{resc.get(c, '.')}\n")
bp = sum(L.get(c, 0) for c in per)
print(f"Genes (compleasm, complete on the unitigs): {len(allg)}; short of {a.ploidy} copies in the groups: {len(short)}; "
      f"of these, with a copy on an unassigned chain: {len(rec)}, on {len(per)} chains, {bp / 1e6:.1f} Mb ({a.out}).")
if a.rescued:
    fixed = sum(1 for g in short if after[g] >= a.ploidy); rc = [c for c in per if c in resc]
    print(f"Rescue ({os.path.basename(a.rescued)}): places {len(resc)} chains, {sum(L.get(c, 0) for c in resc) / 1e6:.1f} Mb; "
          f"{len(rc)} of the {len(per)} chains above; short genes made whole: {fixed} of {len(short)}.")
