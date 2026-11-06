#!/usr/bin/env python3
'''
Reference-free adjacency evidence from homologous haplotypes. The two chains of a junction should have alleles close
together on a homolog, unless an SV separates them on one haplotype.

Homologs of a group are the --homologs groups with the most allelic sequence shared with it. Each chain end is mapped
to each homolog's current layout: the median scaffold coordinate of allele blocks within --endWindow bp of that end
(at least --minBp aligned). Two ends of chains in the same group support a junction on that homolog when their
coordinates are within --maxDist bp.

Low-information chains are ignored (--stats): satellite alignments are not alleles.

Output: end_a end_b homologs_supporting
Usage: synteny_links.py order.tsv alleles.blocks.tsv chains.fai out.tsv --stats chain_stats.tsv
'''
import argparse, collections, statistics
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("order"); p.add_argument("blocks"); p.add_argument("fai"); p.add_argument("out")
p.add_argument("--stats"); p.add_argument("--minWithin", type=float, default=700)
p.add_argument("--homologs", type=int, default=2); p.add_argument("--endWindow", type=int, default=500000)
p.add_argument("--minBp", type=int, default=50000); p.add_argument("--maxDist", type=int, default=1500000)
p.add_argument("--gap", type=int, default=20000)
a = p.parse_args()
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
low = set()
if a.stats:
    low = {l.split("\t")[0] for l in open(a.stats) if l[0] != "#" and float(l.split("\t")[2]) < a.minWithin}
group = {}; layout = {}; pos = collections.defaultdict(int)
for l in open(a.order):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[1] == ".": continue
    g = s[0]; group[s[2]] = g
    layout[s[2]] = (pos[g], s[3]); pos[g] += L[s[2]] + a.gap
def on_scaffold(c, x):
    st, o = layout[c]
    return st + x if o == "+" else st + L[c] - x
B = collections.defaultdict(list); share = collections.Counter()
for l in open(a.blocks):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    x, x0, x1, y, y0, y1 = s[0], int(s[1]), int(s[2]), s[3], int(s[4]), int(s[5])
    if x not in group or y not in group or x in low or y in low or group[x] == group[y]: continue
    B[x].append((x0, x1, y, (y0 + y1) // 2)); B[y].append((y0, y1, x, (x0 + x1) // 2))
    share[(group[x], group[y])] += x1 - x0; share[(group[y], group[x])] += y1 - y0
homologs = {}
for g in set(group.values()):
    homologs[g] = [h for h, _ in sorted(((h, n) for (k, h), n in share.items() if k == g), key=lambda t: -t[1])[:a.homologs]]
def end_pos(c, e, h):
    '''median coordinate on homolog group h of allele blocks near end e of chain c'''
    lo, hi = (L[c] - a.endWindow, L[c]) if e == "+" else (0, a.endWindow)
    hits = []; bp = 0
    for x0, x1, y, ym in B[c]:
        if group[y] == h and x1 > lo and x0 < hi:
            n = min(x1, hi) - max(x0, lo); bp += n
            hits.append(on_scaffold(y, ym))
    return statistics.median(hits) if bp >= a.minBp else None
members = collections.defaultdict(list)
for c, g in group.items():
    if c not in low: members[g].append(c)
n = 0
with open(a.out, "w") as o:
    o.write("#end_a\tend_b\thomologs_supporting\n")
    for g, cs in members.items():
        ends = [(c, e) for c in cs for e in "+-"]
        P = {x: [end_pos(x[0], x[1], h) for h in homologs[g]] for x in ends}
        for i, x in enumerate(ends):
            for y in ends[i + 1:]:
                if x[0] == y[0]: continue
                sup = sum(1 for u, v in zip(P[x], P[y]) if u is not None and v is not None and abs(u - v) <= a.maxDist)
                if sup:
                    o.write(f"{x[0]}{x[1]}\t{y[0]}{y[1]}\t{sup}\n"); n += 1
print(f"{n} end pairs supported by a homolog; {len(members)} groups")
