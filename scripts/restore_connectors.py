#!/usr/bin/env python3
'''Dropped unitigs that connect chain ends. copy_number.py drops a unitig whose windows are covered by kept copies. An arm of a
bubble is often covered by the flanks it overlaps, by 60-130 kb each side, so it is correctly redundant as sequence, but it was also
the graph's only connection between those flanks. Without it the flanks can only be ordered by Hi-C, leaving a gap.
Walks from each chain's terminal unitig through dropped unitigs (up to --hops, never repeat-class) to another chain's terminal unitig.
Every unitig on such a path is listed. Restored unitigs are candidates only: they reach the scaffolds through ranked bridges and
add nothing otherwise.
Usage: restore_connectors.py chains.agp unitigs.gfa copy_number.tsv repeat_chains.txt out.tsv --repeat R [--hops 3]'''
import argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "agp gfa cn rep out".split(): p.add_argument(x)
p.add_argument("--repeat", type=float, required=True); p.add_argument("--hops", type=int, default=3); a = p.parse_args()
inv = {"+": "-", "-": "+"}
parts = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4] == "W": parts[s[0]].append((s[5], s[8]))
inchain = {u for v in parts.values() for u, _ in v}
term = {}                                   # (unitig, outward end) -> chain, for chain ends only
for c, v in parts.items():
    term[(v[-1][0], v[-1][1])] = c; term[(v[0][0], inv[v[0][1]])] = c
adj = collections.defaultdict(list)         # (unitig, end) -> [(unitig, end met)]
for l in open(a.gfa):
    if l[0] == "L":
        s = l.split("\t"); e1, e2 = (s[1], s[2]), (s[3], inv[s[4]]); adj[e1].append(e2); adj[e2].append(e1)
info = {}
for l in open(a.cn):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); info[s[0]] = (int(s[1]), float(s[2]), float(s[5]))
rep = set()
for l in open(a.rep):
    if l.strip(): rep.update(l.rstrip("\n").split("\t")[2].split(","))
usable = lambda u: u not in inchain and u not in rep and u in info and info[u][2] < a.repeat
found = collections.defaultdict(set)
for (u, side), c in term.items():
    stack = [(v, inv[e], [v]) for v, e in adj[(u, side)] if usable(v)]
    while stack:
        x, out, path = stack.pop()
        for y, e in adj[(x, out)]:
            if (y, e) in term and term[(y, e)] != c:
                for z in path: found[z].add(frozenset((c, term[(y, e)])))
            elif usable(y) and y not in path and len(path) < a.hops: stack.append((y, inv[e], path + [y]))
with open(a.out, "w") as o:
    o.write("#unitig\tlength\tprimary\tclass_depth\tchain_pairs\n")
    for u in sorted(found, key=lambda x: -info[x][0]):
        n, d, k = info[u]; o.write(f"{u}\t{n}\t{d}\t{k}\t" + ",".join(sorted("|".join(sorted(f)) for f in found[u])) + "\n")
print(f"{len(found)} connector unitigs, {sum(info[u][0] for u in found) / 1e6:.2f} Mb, between {len({f for v in found.values() for f in v})} chain pairs")
