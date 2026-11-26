#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Stage 07: place chains that grouping left unassigned, by allele-rescued Hi-C (settings.hic_allele_rescue). Runs after hic_group.py.
Chromosome: links from rescued pairs (hic_rescue_pairs.py) to chains already in each chromosome's groups, per Mb of those chains; the
best chromosome needs >= --minLinks links and >= --ratio times the next, as hic_group.py's step 3. Rounds repeat with the chains placed
so far, so a block of chains can follow its anchored end inward.
Haplotype, in this order: a group holding an allele of the chain is excluded (call allele); else unique Hi-C pairs to each group,
>= --ratio times the other (call hic); else, with --endLinks, the one group the chain has graph links to (>= 5 kb overlap or >= 3 reads,
as make_order_links.py; call graph: a tie-break only, since near-identical alleles can both link, which is why a link to every group
decides nothing and a confident Hi-C call is never overruled); else the group with less sequence placed by this script so far (call
arbitrary: the phase is not known, as with arbitrary arms at D15). Longer chains go first, so an allele pair is split between the groups.
Call allele_of_arbitrary: the excluding allele was itself placed arbitrarily by this script, so the call carries no phase evidence.
A chain with an allele in every group of its chromosome is reported once (allele_in_every_group) and left unassigned.
Two more calls leave a chain unassigned, tested before the haplotype, once its chromosome is known:
  collapsed  with --hubCores: a breakable row in review/hub_cores.py's table (a 2x unique-read core with a unitig boundary at each inner edge).
             The chain holds both haplotypes there, so no haplotype call can be right: its unique Hi-C splits evenly and it fell to arbitrary
             (292: chain0115, chain0119). It needs a D06 break and a shared core, not a haplotype. A core inside one unitig cannot be broken
             off, so that chain is placed as before, on one haplotype, and stage 07's report lists it for D14.
  redundant  with --blocks: allele blocks of any length (alleles.blocks.tsv), merged, cover >= --redundant of the chain in every group of the
             chromosome: a near-copy of sequence both haplotypes already hold. The allele rule counts only alleles of >= --minConflict bp, so a short
             chain whose copies are each under that passed it (292: chain0288, 83 kb, covered 0.79 and 0.65 by blocks of 30-83 kb, dropped
             between a chain and its telomeric tip). Tested when the chain's chromosome is found, and again after the last round.
--dovetails (allele_table.py's dovetails.tsv): near-identical end overlaps where one chain runs on into the other. They are graph
links, added to --endLinks' for the graph call, never alleles. As with any link, a hub reaches both haplotypes, so a link to every
group decides nothing (a homozygous stretch: one chain's start overlaps chains of both haplotypes). With --dovetails, a chain whose
call would be arbitrary or allele_of_arbitrary while it links to every group of its chromosome is left unassigned (call fork): it sits
where the haplotypes meet, with nothing to say which it continues. Graph calls are tested again after the last round, since a shorter chain placed later can be the link to the
other group. Without the dovetails, the allele rule had excluded such chains by
mistake (allele_in_every_group or redundant), so this keeps them out on the right grounds. Its chromosome call stands for later rounds.
Dovetails still count towards redundant: a chain made of one group chain's end and the next one's start (a 160 kb chain identical to the two
flanks of a junction, in a trial) holds nothing that group lacks, and placing it as an allele duplicates both overlaps.
Two more, with --dovetails:
  bubble_arm  a dovetail into an end of chain x that a chain y of x's group already continues, where the candidate is also identical
              (de <= --copyDe) to part of y: the candidate is y's other-haplotype arm, rejoining the path y takes. It goes to the other
              group(s) of the chromosome. A dovetail into an end already taken in that group never votes for that group.
  copy        with --depthWindows and --hap: blocks at de <= --copyDe with one placed chain cover >= --copyFrac of the candidate, and the
              two read under --copyDepth x hap together over the shared stretch (median 10 kb window depth): one haplotype's sequence
              twice, not an allele. Left unassigned. Summed near 2x it is a homozygous stretch, and the allele call stands (in a trial,
              three such pairs summed 1.55-1.85x and one duplicate 0.79x).
Shared and high-repeat chains are left alone. Rewrites prefix.tsv (the original is kept as prefix.hic.tsv) and writes prefix.rescued.tsv.
Usage: rescue_groups.py groups_prefix rescued_pairs.tsv.gz unique_pairs.tsv.gz alleles.pairs.tsv chain_repeat.tsv
       [--hubCores hub_cores.tsv] [--blocks alleles.blocks.tsv] [--endLinks end_links.tsv] [--dovetails alleles.dovetails.tsv]'''
import argparse, collections, gzip, os, shutil, sys
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "prefix rescued unique alleles repeat".split(): p.add_argument(x)
p.add_argument("--minLinks", type=int, default=10); p.add_argument("--ratio", type=float, default=2.0); p.add_argument("--rounds", type=int, default=5)
p.add_argument("--minConflict", type=int, default=50000); p.add_argument("--minConflictFrac", type=float, default=0.1); p.add_argument("--maxRepeat", type=float, default=0.5)
p.add_argument("--hubCores"); p.add_argument("--blocks"); p.add_argument("--endLinks"); p.add_argument("--dovetails"); p.add_argument("--redundant", type=float, default=0.5)
p.add_argument("--depthWindows"); p.add_argument("--hap", type=float); p.add_argument("--copyDe", type=float, default=0.002)
p.add_argument("--copyFrac", type=float, default=0.8); p.add_argument("--copyDepth", type=float, default=1.25)
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
HUB = set()
if a.hubCores:
    for l in open(a.hubCores):
        if l[0] != "#" and l.rstrip("\n").split("\t")[-1] == "yes": HUB.add(l.split("\t")[0])
GL = collections.defaultdict(set)                            # candidate -> chains it has a graph link to (end_links.tsv)
if a.endLinks:
    for l in open(a.endLinks):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t"); x, y = s[0][:-1], s[1][:-1]
        if x != y and (int(s[4]) >= 5000 or int(s[2]) >= 3):
            if x in cand: GL[x].add(y)
            if y in cand: GL[y].add(x)
if a.dovetails:
    for l in open(a.dovetails):
        if l[0] == "#": continue
        s = l.split("\t"); x, y = s[0], s[3]
        if x != y:
            if x in cand: GL[x].add(y)
            if y in cand: GL[y].add(x)
BL = collections.defaultdict(list)                           # chain -> [(start, end, other chain)], allele blocks of any length
for f in (a.blocks, a.dovetails):   # a dovetail is no allele, but it is sequence a group already holds: it counts towards redundant
    if not f: continue
    for l in open(f):
        if l[0] == "#": continue
        s = l.split("\t")
        if s[0] in cand: BL[s[0]].append((int(s[1]), int(s[2]), s[3]))
        if s[3] in cand: BL[s[3]].append((int(s[4]), int(s[5]), s[0]))
IB = collections.defaultdict(list)                           # chain -> [(start, end, other, other_start, other_end)], blocks at de <= copyDe
TAKE = collections.defaultdict(set)                         # (chain, end h/t) -> chains that dovetail into that end
for f in (a.blocks, a.dovetails):
    if not f: continue
    for l in open(f):
        if l[0] == "#": continue
        s = l.split("\t"); x, y = s[0], s[3]; xs, xe, ys, ye = int(s[1]), int(s[2]), int(s[4]), int(s[5])
        if x == y or float(s[6]) > a.copyDe: continue
        IB[x].append((xs, xe, y, ys, ye)); IB[y].append((ys, ye, x, xs, xe))
        if f == a.dovetails:
            for u, us, ue, v in ((x, xs, xe, y), (y, ys, ye, x)):
                if u in L: TAKE[(u, "h" if us <= 1000 else "t" if ue >= L[u] - 1000 else "m")].add(v)
def taken(c, x):   # c dovetails into an end of x that another chain of x's group already continues: [those chains]
    return [y for e in "ht" if c in TAKE[(x, e)] for y in TAKE[(x, e)] if y != c and grp.get(y) == grp.get(x) and grp.get(x) != "unassigned"]
def bubble_arm(c):   # the group chain c is the other-haplotype arm of, or None
    for x in GL[c]:
        for y in taken(c, x):
            if any(o == y for _, _, o, _, _ in IB[c]): return y
    return None
WD = collections.defaultdict(list)
if a.depthWindows:
    for l in open(a.depthWindows):
        if l[0] == "#": continue
        s = l.split("\t")
        if s[0] in L: WD[s[0]].append((int(s[1]), int(s[2]), float(s[3])))
def wdepth(c, x, y):
    v = sorted(d for s, e, d in WD[c] if s < y and e > x)
    return v[len(v) // 2] if v else None
def copy_of(c):   # a placed chain c duplicates (one haplotype's sequence twice), or None
    if not (a.depthWindows and a.hap): return None
    by = collections.defaultdict(list)
    for s, e, o, os_, oe in IB[c]: by[o].append((s, e, os_, oe))
    for o, bl in by.items():
        if grp.get(o, "unassigned") == "unassigned": continue
        n = f = 0
        for s, e, _, _ in sorted(bl):
            s = max(s, f)
            if e > s: n += e - s; f = e
        if n < a.copyFrac * L[c]: continue
        dc = wdepth(c, min(b[0] for b in bl), max(b[1] for b in bl)); do = wdepth(o, min(b[2] for b in bl), max(b[3] for b in bl))
        if dc is not None and do is not None and dc + do < a.copyDepth * a.hap: return o
    return None
def covered(c, g):
    n = e = 0
    for x, y, o in sorted(BL[c]):
        if o == c or grp.get(o) != g: continue
        x = max(x, e)
        if y > x: n += y - x; e = y
    return n / L[c]
placed_bp = collections.Counter(); out = []; blocked = set(); how = {}   # how: chain -> call, for chains placed here
on = {}                                                       # chain -> chromosome for collapsed and redundant chains: no haplotype, but their
                                                              # chromosome call stands, so later rounds may still anchor on them
for rnd in range(1, a.rounds + 1):
    clen = collections.Counter()
    for c, g in grp.items():
        if g in chrom: clen[chrom[g]] += L[c]
    for c, ch in on.items(): clen[ch] += L[c]
    new = []
    for c in sorted(cand, key=lambda c: -L[c]):
        if grp[c] != "unassigned" or c in blocked: continue
        k = collections.Counter()
        for o, n in R[c].items():
            if grp.get(o) in chrom: k[chrom[grp[o]]] += n
            elif o in on: k[on[o]] += n
        d = sorted(((n / (clen[ch] / 1e6), n, ch) for ch, n in k.items() if clen[ch]), reverse=True)
        if not d or d[0][1] < a.minLinks: continue
        r = d[0][0] / d[1][0] if len(d) > 1 and d[1][0] > 0 else float("inf")
        if r < a.ratio: continue
        new.append((c, d[0][2], d[0][1], r))
    if not new: break
    for c, ch, n, r in new:
        if c in HUB:
            out.append((c, L[c], ch, "unassigned", "collapsed", n, r, rnd)); blocked.add(c); on[c] = ch; continue
        if BL and all(covered(c, g) >= a.redundant for g in groups_of[ch]):
            out.append((c, L[c], ch, "unassigned", "redundant", n, r, rnd)); blocked.add(c); on[c] = ch; continue
        if copy_of(c):
            out.append((c, L[c], ch, "unassigned", "copy", n, r, rnd)); blocked.add(c); on[c] = ch; continue
        gs = [g for g in groups_of[ch] if not any(grp.get(x) == g for x in allele[c])]
        if not gs: out.append((c, L[c], ch, "unassigned", "allele_in_every_group", n, r, rnd)); blocked.add(c); continue
        ex = [x for x in allele[c] if grp.get(x) in groups_of[ch] and grp.get(x) not in gs]
        call = None
        if len(gs) < len(groups_of[ch]):   # allele exclusion; not independent when every excluding allele was itself placed arbitrarily here
            call = "allele_of_arbitrary" if ex and all(how.get(x) in ("arbitrary", "allele_of_arbitrary") for x in ex) else "allele"
        if len(gs) > 1:
            u = sorted(((sum(m for o, m in U[c].items() if grp.get(o) == g), g) for g in gs), reverse=True)
            gl = [g for g in gs if any(grp.get(x) == g and not taken(c, x) for x in GL[c])]
            arm = bubble_arm(c) if a.dovetails else None
            ag = [g for g in gs if arm and g != grp[arm]]
            if u[0][0] >= a.minLinks and u[0][0] >= a.ratio * max(u[1][0], 1): gs, call = [u[0][1]], "hic"
            elif len(ag) == 1: gs, call = ag, "bubble_arm"
            elif len(gl) == 1: gs, call = gl, "graph"
            else: gs, call = [min(gs, key=lambda g: (placed_bp[g], g))], "arbitrary"
        if a.dovetails and call in ("arbitrary", "allele_of_arbitrary") and all(any(grp.get(x) == g for x in GL[c]) for g in groups_of[ch]):
            out.append((c, L[c], ch, "unassigned", "fork", n, r, rnd)); blocked.add(c); on[c] = ch; continue
        grp[c] = gs[0]; placed_bp[gs[0]] += L[c]; how[c] = call; out.append((c, L[c], ch, gs[0], call, n, r, rnd))
if a.dovetails:   # graph calls again once every round is done: chains go longest first, so a shorter chain placed later in the other group
                 # can link to the same chain (in a trial, a chain whose start both haplotypes run into). A graph call that now
                 # links to every group of its chromosome is a fork. Shortest first, as below.
    idx = {x[0]: i for i, x in enumerate(out)}
    for c in sorted((x[0] for x in out if x[4] == "graph"), key=lambda c: (L[c], c)):
        x = out[idx[c]]
        if all(any(grp.get(y) == g for y in GL[c]) for g in groups_of[x[2]]):
            placed_bp[grp[c]] -= L[c]; grp[c] = "unassigned"; on[c] = x[2]
            out[idx[c]] = (c, x[1], x[2], "unassigned", "fork", x[5], x[6], x[7])
if BL:   # again once every round is done: a copy placed after the chain can make it redundant. Shortest first, and a chain unassigned here
         # stops covering others, so two copies of one sequence are not both dropped.
    idx = {x[0]: i for i, x in enumerate(out)}
    for c in sorted((x[0] for x in out if x[3] != "unassigned"), key=lambda c: (L[c], c)):
        x = out[idx[c]]
        if all(covered(c, g) >= a.redundant for g in groups_of[x[2]]):
            placed_bp[grp[c]] -= L[c]; grp[c] = "unassigned"; on[c] = x[2]
            out[idx[c]] = (c, x[1], x[2], "unassigned", "redundant", x[5], x[6], x[7])
shutil.copyfile(a.prefix + ".tsv", a.prefix + ".hic.tsv")
with open(a.prefix + ".tsv", "w") as o:
    for h in hdr: o.write("\t".join(h) + "\n")
    for r in rows: o.write("\t".join([r[0], r[1], grp[r[0]]] + r[3:]) + "\n")
with open(a.prefix + ".rescued.tsv", "w") as o:
    o.write("#chain\tlength\tchrom\tgroup\tcall\trescued_links\tratio_to_next\tround\n")
    for c, n, ch, g, call, k, r, rnd in out: o.write(f"{c}\t{n}\t{ch}\t{g}\t{call}\t{k}\t{r:.2f}\t{rnd}\n")
done = [x for x in out if x[3] != "unassigned"]
left = collections.Counter(x[4] for x in out if x[3] == "unassigned")
print(f"# rescue: {len(done)} chains, {sum(x[1] for x in done)/1e6:.1f} Mb placed in {rnd} rounds; calls {dict(collections.Counter(x[4] for x in done))}; left unassigned {dict(left) or 0}")
