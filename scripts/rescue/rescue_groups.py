#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Place chains that grouping left unassigned, by allele-rescued Hi-C (settings.hic_allele_rescue). Runs after hic_group.py.
Chromosome: links from rescued pairs (hic_rescue_pairs.py) to chains already in each chromosome's groups, per Mb of those chains; the
best chromosome needs >= --minLinks links and >= --ratio times the next, as hic_group.py's step 3. Rounds repeat with the chains placed
so far, so a block of chains can follow its anchored end inward.
Haplotype, in this order: a group holding an allele of the chain is excluded (call allele); else unique Hi-C pairs to each group,
>= --ratio times the other (call hic); else the group with less sequence placed by this script so far (call arbitrary: the phase is not
known, as with arbitrary arms at D15). Longer chains go first, so an allele pair is split between the groups.
Call allele_of_arbitrary: the excluding allele was itself placed arbitrarily by this script, so the call carries no phase evidence.
A chain with an allele in every group of its chromosome is reported once (allele_in_every_group) and left unassigned.
Shared and high-repeat chains are left alone. Rewrites prefix.tsv (the original is kept as prefix.hic.tsv) and writes prefix.rescued.tsv.
Usage: rescue_groups.py groups_prefix rescued_pairs.tsv.gz unique_pairs.tsv.gz alleles.pairs.tsv chain_repeat.tsv'''
import argparse, collections, gzip, os, shutil
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "prefix rescued unique alleles repeat".split(): p.add_argument(x)
p.add_argument("--minLinks", type=int, default=10); p.add_argument("--ratio", type=float, default=2.0); p.add_argument("--rounds", type=int, default=5)
p.add_argument("--minConflict", type=int, default=50000); p.add_argument("--minConflictFrac", type=float, default=0.1); p.add_argument("--maxRepeat", type=float, default=0.5)
a = p.parse_args()
rows = [l.rstrip("\n").split("\t") for l in open(a.prefix + ".tsv")]
hdr = [r for r in rows if r[0].startswith("#")]; rows = [r for r in rows if not r[0].startswith("#")]
L = {r[0]: int(r[1]) for r in rows}; grp = {r[0]: r[2] for r in rows}
chrom = dict(l.split()[:2] for l in open(a.prefix + ".chroms.tsv") if l[0] != "#")
groups_of = collections.defaultdict(list)
for g, c in chrom.items(): groups_of[c].append(g)
rep = {}
for l in open(a.repeat):
    if l[0] != "#": s = l.split("\t"); rep[s[0]] = float(s[2])
allele = collections.defaultdict(set)
for l in open(a.alleles):
    if l[0] == "#": continue
    s = l.split("\t"); x, y, bx, by = s[0], s[1], int(s[2]), int(s[3])
    if x in L and y in L and min(bx, by) >= a.minConflict and min(bx, by) >= a.minConflictFrac * min(L[x], L[y]) and max(rep.get(x, 0), rep.get(y, 0)) <= a.maxRepeat:
        allele[x].add(y); allele[y].add(x)
cand = {c for c, g in grp.items() if g == "unassigned" and rep.get(c, 0) <= a.maxRepeat}
R = collections.defaultdict(collections.Counter)            # candidate -> other chain -> rescued links
for l in gzip.open(a.rescued, "rt"):
    s = l.split("\t"); x, y = s[0], s[2]
    if x == y: continue
    if x in cand: R[x][y] += 1
    if y in cand: R[y][x] += 1
U = collections.defaultdict(collections.Counter)
for l in gzip.open(a.unique, "rt"):
    s = l.split("\t"); x, y = s[0], s[2]
    if x == y: continue
    if x in cand: U[x][y] += 1
    if y in cand: U[y][x] += 1
placed_bp = collections.Counter(); out = []; blocked = set(); how = {}   # how: chain -> call, for chains placed here
for rnd in range(1, a.rounds + 1):
    clen = collections.Counter()
    for c, g in grp.items():
        if g in chrom: clen[chrom[g]] += L[c]
    new = []
    for c in sorted(cand, key=lambda c: -L[c]):
        if grp[c] != "unassigned" or c in blocked: continue
        k = collections.Counter()
        for o, n in R[c].items():
            if grp.get(o) in chrom: k[chrom[grp[o]]] += n
        d = sorted(((n / (clen[ch] / 1e6), n, ch) for ch, n in k.items() if clen[ch]), reverse=True)
        if not d or d[0][1] < a.minLinks: continue
        r = d[0][0] / d[1][0] if len(d) > 1 and d[1][0] > 0 else float("inf")
        if r < a.ratio: continue
        new.append((c, d[0][2], d[0][1], r))
    if not new: break
    for c, ch, n, r in new:
        gs = [g for g in groups_of[ch] if not any(grp.get(x) == g for x in allele[c])]
        if not gs: out.append((c, L[c], ch, "unassigned", "allele_in_every_group", n, r, rnd)); blocked.add(c); continue
        ex = [x for x in allele[c] if grp.get(x) in groups_of[ch] and grp.get(x) not in gs]
        call = None
        if len(gs) < len(groups_of[ch]):   # allele exclusion; not independent when every excluding allele was itself placed arbitrarily here
            call = "allele_of_arbitrary" if ex and all(how.get(x) in ("arbitrary", "allele_of_arbitrary") for x in ex) else "allele"
        if len(gs) > 1:
            u = sorted(((sum(m for o, m in U[c].items() if grp.get(o) == g), g) for g in gs), reverse=True)
            if u[0][0] >= a.minLinks and u[0][0] >= a.ratio * max(u[1][0], 1): gs, call = [u[0][1]], "hic"
            else: gs, call = [min(gs, key=lambda g: (placed_bp[g], g))], "arbitrary"
        grp[c] = gs[0]; placed_bp[gs[0]] += L[c]; how[c] = call; out.append((c, L[c], ch, gs[0], call, n, r, rnd))
shutil.copyfile(a.prefix + ".tsv", a.prefix + ".hic.tsv")
with open(a.prefix + ".tsv", "w") as o:
    for h in hdr: o.write("\t".join(h) + "\n")
    for r in rows: o.write("\t".join([r[0], r[1], grp[r[0]]] + r[3:]) + "\n")
with open(a.prefix + ".rescued.tsv", "w") as o:
    o.write("#chain\tlength\tchrom\tgroup\tcall\trescued_links\tratio_to_next\tround\n")
    for c, n, ch, g, call, k, r, rnd in out: o.write(f"{c}\t{n}\t{ch}\t{g}\t{call}\t{k}\t{r:.2f}\t{rnd}\n")
done = [x for x in out if x[3] != "unassigned"]
print(f"# rescue: {len(done)} chains, {sum(x[1] for x in done)/1e6:.1f} Mb placed in {rnd} rounds; calls {dict(collections.Counter(x[4] for x in done))}")
