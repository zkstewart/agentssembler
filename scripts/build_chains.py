#!/usr/bin/env python3
'''
Join kept unitigs into chains along reciprocal one-to-one graph links, trimming each link's overlap.
Chains are proposals: they are verified afterwards by reads aligning contiguously through each join.

A link is used only if each of its two ends has exactly one link to another allowed unitig
(--breaks removes links first, one "end_a end_b" per line as in read_links.py).
Ends: u+ = u's 3' end, u- = u's 5' end.

Usage: build_chains.py --keep kept.txt --gfa unitigs.gfa --fai unitigs.fai --out prefix [--breaks f]
Writes prefix.agp (components are unitig ranges, no gaps), prefix.joins.tsv, prefix.circular.txt, and
prefix.contained.tsv, listing unitigs wholly inside an overlap longer than themselves, which add no sequence
'''
import argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--keep", required=True); p.add_argument("--gfa", required=True); p.add_argument("--fai", required=True)
p.add_argument("--out", required=True); p.add_argument("--breaks"); p.add_argument("--exclude", help="unitigs never joined")
a = p.parse_args()
keep = [ l.strip() for l in open(a.keep) if l.strip() ]; K = set(keep)
ex = set(l.strip() for l in open(a.exclude)) if a.exclude else set()
L = { l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai) }
inv = { "+": "-", "-": "+" }
brk = set()
if a.breaks:
    for l in open(a.breaks):
        if l.strip() and l[0] != "#":
            x, y = l.split()[:2]; brk.add(frozenset(((x[:-1], x[-1]), (y[:-1], y[-1]))))
adj = collections.defaultdict(dict) # end -> other end -> overlap
for l in open(a.gfa):
    if l[0] != "L": continue
    s = l.split("\t")
    if s[1] not in K or s[3] not in K or s[1] == s[3] or s[1] in ex or s[3] in ex: continue
    e1, e2 = (s[1], s[2]), (s[3], inv[s[4]])
    if frozenset((e1, e2)) in brk: continue
    n = int(s[5].strip()[:-1]); adj[e1][e2] = n; adj[e2][e1] = n
join = {}
for e, nb in adj.items():
    if len(nb) == 1:
        o = next(iter(nb))
        if len(adj[o]) == 1: join[e] = o
# walk from chain starts: a unitig whose 5' or 3' end is unjoined
seen, chains, circ = set(), [], []
def walk(u, st):
    '''start at u oriented st (we leave through (u, st)); returns [(unitig, strand, overlap with previous)]'''
    path = [(u, st, 0)]; seen.add(u)
    while (path[-1][0], path[-1][1]) in join:
        c, d = join[(path[-1][0], path[-1][1])]
        if c in seen: break
        ov = adj[(path[-1][0], path[-1][1])][(c, d)]
        path.append((c, inv[d], ov)); seen.add(c)
    return path
for u in keep:
    if u in seen: continue
    if (u, "-") not in join: chains.append(walk(u, "+"))
    elif (u, "+") not in join: chains.append(walk(u, "-"))
for u in keep: # what is left is circular
    if u not in seen:
        p_ = walk(u, "+"); chains.append(p_); circ.append(p_)
chains.sort(key=lambda c: -sum(L[x[0]] - x[2] for x in c))
agp = open(a.out + ".agp", "w"); jt = open(a.out + ".joins.tsv", "w"); agp.write("##agp-version 2.1\n")
jt.write("#chain\tpos\tleft\tright\toverlap\n"); tot = 0; contained = []
for i, c in enumerate(chains, 1):
    name = f"chain{i:04d}"; pos = 0; carry = 0; part = 0
    for k, (u, st, ov) in enumerate(c):
        # An overlap can be longer than the unitig it enters. That unitig is then wholly inside sequence
        # already written, and what is left of the overlap carries into the next one, which overlaps this unitig's tail.
        if k: jt.write(f"{name}\t{pos}\t{c[k-1][0]}{c[k-1][1]}\t{u}{st}\t{ov}\n")
        t = ov + carry
        if t >= L[u]: carry = t - L[u]; contained.append((name, u, L[u], ov)); continue
        carry = 0; part += 1
        cs, ce = (t + 1, L[u]) if st == "+" else (1, L[u] - t)
        n = ce - cs + 1
        agp.write(f"{name}\t{pos+1}\t{pos+n}\t{part}\tW\t{u}\t{cs}\t{ce}\t{st}\n"); pos += n
    tot += pos
with open(a.out + ".contained.tsv", "w") as o:
    o.write("#chain\tunitig\tlength\toverlap\n")
    for x in contained: o.write("\t".join(map(str, x)) + "\n")
with open(a.out + ".circular.txt", "w") as o:
    for c in circ: o.write(",".join(x[0] + x[1] for x in c) + "\n")
print(f"{len(chains)} chains from {len(keep)} unitigs, {tot/1e6:.1f} Mb; {sum(len(c)>1 for c in chains)} with joins, {len(circ)} circular; largest {max(sum(L[x[0]]-x[2] for x in c) for c in chains)/1e6:.2f} Mb")
