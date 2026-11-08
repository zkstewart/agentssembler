#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Gene completeness across the whole genome, not by haplotype number. Haplotype numbers are arbitrary across chromosomes, so with a translocation "all h1 scaffolds" need not be a complete set:
a segment missing from chr2_h1 and carried by chr3_h2 counts as lost to h1 although the genome has it. Reports complete copies per gene on any scaffold, then the numbering per chromosome that
makes the most genes complete per set. A large gain means the per-haplotype score was measuring numbering, provided the translocations themselves stand (D12). Scaffolds are not renamed.
Usage: complement_genes.py PLOIDY table_h1.tsv [table_h2.tsv ...] [--unplaced table.tsv]   (compleasm full_table_busco_format.tsv files)'''
import argparse, collections, itertools, random
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("ploidy", type=int); p.add_argument("tables", nargs="+"); p.add_argument("--unplaced"); p.add_argument("--starts", type=int, default=30); a = p.parse_args()
G = collections.defaultdict(set); ALL = set(); UN = set()
def read(f, un):
    for l in open(f):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t"); ALL.add(s[0])
        if s[1] in ("Complete", "Duplicated") and len(s) > 4:
            if un: UN.add(s[0])
            else: G[s[0]].add(s[2])
for f in a.tables: read(f, False)
if a.unplaced: read(a.unplaced, True)
N = len(ALL); haps = collections.defaultdict(set)
for ss in G.values():
    for s in ss:
        c, _, h = s.rpartition("_h")
        if c: haps[c].add(h)
C = sorted(c for c, h in haps.items() if len(h) == a.ploidy)
def score(x):
    sets = [{f"{c}_h{x[c][k]}" for c in C} for k in range(a.ploidy)]
    return [sum(1 for ss in G.values() if ss & m) for m in sets]
P = {c: list(itertools.permutations(sorted(haps[c]))) for c in C}
base = {c: tuple(sorted(haps[c])) for c in C}; s0 = score(base); random.seed(1); best = None
for _ in range(a.starts):
    x = {c: random.choice(P[c]) for c in C}; v = sum(score(x)); up = True
    while up:
        up = False
        for c in C:
            for q in P[c]:
                y = dict(x); y[c] = q; w = sum(score(y))
                if w > v: x, v, up = y, w, True
    if best is None or v > best[0]: best = (v, x)
fmt = lambda s: ", ".join(f"{100 * z / N:.2f}%" for z in s)
print(f"Complete per haplotype number: {fmt(s0)}; not complete in all sets {sum(N - z for z in s0)}.")
s1 = score(best[1]); print(f"Complete per set, best numbering: {fmt(s1)}; not complete in all sets {sum(N - z for z in s1)}.")
for k in range(a.ploidy): print(f"set{k + 1}: " + " ".join(f"{c}_h{best[1][c][k]}" for c in C))
cc = collections.Counter(min(len(G.get(g, ())), a.ploidy + 1) for g in ALL); short = sum((a.ploidy - k) * v for k, v in cc.items() if k < a.ploidy)
low = [g for g in ALL if len(G.get(g, ())) < a.ploidy]
print(f"Gene copies short of {a.ploidy} per gene, on any scaffold: {short} of {a.ploidy * N} ({100 * short / (a.ploidy * N):.2f}%), in {len(low)} genes; {sum(g in UN for g in low)} of them complete in unplaced; absent from all scaffolds {cc[0]}.")
