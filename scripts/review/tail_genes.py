#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Stage 10, for D16: genes on each tail candidate, from compleasm on the unitigs (stage 01). genes: complete genes with a copy on the tail's
chains. short: those of them whose copies on placed chains of the pre-build (placement status placed*) number under the ploidy, the
genes this tail would add a missing copy of. A tail with short genes is worth the tests; it is not a reason to skip them, and a tail
that fails them still keeps its genes in the assembly as unlocalized sequence (stage 15, unloc/).
--unplacedOut: the same count for every chain the pre-build leaves unplaced (no copy placed, not a dropped connector) and not on a tail,
sorted by short genes, for the D16 tangle check: which unplaced chains, and which lost routes (--forks: chain_forks.tsv LOST rows; the
lost_from column names the placed chain whose fork leads there), carry genes that are missing elsewhere. On 292 that was worked out by hand.
Usage: tail_genes.py --genes full_table_busco_format.tsv --agp chains.agp --placement pre.placement.tsv --ploidy N
       [--tails tails.tsv --out tail_genes.tsv] [--unplacedOut unplaced_genes.tsv [--forks chain_forks.tsv]]'''
import argparse, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chainlib import read_agp, gene_hits
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("genes", "agp", "placement"): p.add_argument("--" + x, required=True)
p.add_argument("--tails"); p.add_argument("--out"); p.add_argument("--unplacedOut"); p.add_argument("--forks")
p.add_argument("--ploidy", type=int, required=True)
a = p.parse_args()
if a.out and not a.tails: p.error("--out needs --tails")
base = lambda c: c.split(".g")[0]
uses = collections.Counter()                                  # base chain -> placed copies
status = {}; length = {}
for l in open(a.placement):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4].startswith("placed"): uses[base(s[0])] += 1
        status.setdefault(base(s[0]), s[4]); length[base(s[0])] = s[1]
P, _ = read_agp(a.agp)
on = collections.defaultdict(set); placed = collections.Counter()
for g, st, c, x in gene_hits(a.genes, P):
    if st == "Fragmented" or c is None: continue
    on[c].add(g); placed[g] += uses[c]
T = [l.rstrip("\n").split("\t") for l in open(a.tails) if l[0] != "#" and l.strip()] if a.tails else []
if a.unplacedOut:
    on_tail = {base(z[:-1]) for r in T for z in r[2].split(",")}
    lost = collections.defaultdict(set)
    if a.forks:
        for l in open(a.forks):
            if l.startswith("LOST\t"): s = l.rstrip("\n").split("\t"); lost[base(s[9])].add(base(s[1]))
    R = []
    for c, st in status.items():
        if uses[c] or on_tail & {c} or st.startswith("dropped"): continue
        gs = on[c]; sh = {g for g in gs if placed[g] < a.ploidy}
        R.append((len(sh), len(gs), c, length[c], st, ",".join(sorted(lost[c])) or "."))
    R.sort(key=lambda r: (-r[0], -r[1], r[2]))
    with open(a.unplacedOut, "w") as o:
        o.write("#chain\tlength\tstatus\tgenes\tshort\tlost_from\n")
        for sh, n, c, L, st, lf in R: o.write(f"{c}\t{L}\t{st}\t{n}\t{sh}\t{lf}\n")
    print(f"Unplaced chains with genes: {sum(1 for r in R if r[1])} of {len(R)}; genes {sum(r[1] for r in R)}, of which short of {a.ploidy} placed copies {sum(r[0] for r in R)} ({a.unplacedOut}).")
if not a.out: sys.exit(0)
n_t = n_g = n_s = 0
with open(a.out, "w") as o:
    o.write("#end\tgroup\tlabel\tbp\tgenes\tshort\n")
    for r in T:
        gs = set().union(*(on[base(z[:-1])] for z in r[2].split(","))); sh = {g for g in gs if placed[g] < a.ploidy}
        o.write(f"{r[0]}\t{r[1]}\t{r[3]}\t{r[5]}\t{len(gs)}\t{len(sh)}\n")
        if gs: n_t += 1; n_g += len(gs); n_s += len(sh)
print(f"Tails with genes: {n_t} of {len(T)}; genes {n_g}, of which short of {a.ploidy} placed copies {n_s} ({a.out}).")
