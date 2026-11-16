#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Find a tangle in the hifiasm unitig graph around seed unitigs: BFS over GFA links (both orientations), not entering
"anchor" unitigs (long, haplotype depth). Anchors met are boundary candidates. Reports tangle size, the anchors with the
side of each that faces the tangle, and where every unitig sits in the agentssembler build (chain, scaffold or unplaced).
Usage: tangle_explore.py unitigs.gfa depth.unitigs.tsv chains.agp placement.tsv seed[,seed...] [--minAnchor 300000]
       --hap HAP_DEPTH [--lo 0.6] [--hi 1.4] [--out prefix]'''
import argparse, collections
p = argparse.ArgumentParser()
for x in "gfa depth agp placement seeds".split(): p.add_argument(x)
p.add_argument("--minAnchor", type=int, default=300000); p.add_argument("--hap", type=float, required=True)
p.add_argument("--lo", type=float, default=0.6); p.add_argument("--hi", type=float, default=1.4); p.add_argument("--out")
a = p.parse_args()
L = {}; adj = collections.defaultdict(set)   # adj[(u, end)] = {(v, end)}; end 'h' head(start) or 't' tail(end)
for l in open(a.gfa):
    if l[0] == "S":
        s = l.split("\t", 3); L[s[1]] = int(next((x[5:] for x in l.split("\t")[3:] if x.startswith("LN:i:")), len(s[2])))
    elif l[0] == "L":
        s = l.split("\t"); u, uo, v, vo = s[1], s[2], s[3], s[4]
        eu = "t" if uo == "+" else "h"; ev = "h" if vo == "+" else "t"
        adj[(u, eu)].add((v, ev)); adj[(v, ev)].add((u, eu))
D = {}
for l in open(a.depth):
    if l[0] != "#": s = l.split("\t"); D[s[0]] = float(s[2])
U2C = {}
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.split("\t")
    if s[4] == "W": U2C[s[5]] = s[0]
PL = {}
for l in open(a.placement):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); PL[s[0]] = (s[3] if s[4].startswith("placed") else s[4])
def anchor(u): return L[u] >= a.minAnchor and a.lo * a.hap <= D.get(u, 0) <= a.hi * a.hap
seeds = a.seeds.split(","); seen = set(seeds); q = list(seeds); hits = collections.defaultdict(set)
while q:
    u = q.pop()
    for e in "ht":
        for v, ev in adj[(u, e)]:
            if anchor(v) and v not in seeds: hits[v].add(ev); continue
            if v not in seen: seen.add(v); q.append(v)
bp = sum(L[u] for u in seen)
print(f"# tangle: {len(seen)} unitigs, {bp/1e6:.2f} Mb; anchors met: {len(hits)}")
def where(u):
    c = U2C.get(u, "-"); return f"{c}\t{PL.get(c, PL.get(c.split('.')[0], 'dropped'))}"
print("#anchor\tlength\tdepth\tend_facing_tangle\tchain\tplacement")
for v in sorted(hits, key=lambda x: -L[x]): print(f"{v}\t{L[v]}\t{D.get(v,0):.1f}\t{','.join(sorted(hits[v]))}\t{where(v)}")
byp = collections.Counter(); byc = collections.Counter()
for u in seen: byp[where(u).split("\t")[1]] += L[u]; byc[where(u).split("\t")[0]] += 1
print("# tangle bp by placement: " + "; ".join(f"{k} {v/1e6:.2f} Mb" for k, v in byp.most_common()))
mult = collections.Counter()
for u in seen:
    r = D.get(u, 0) / a.hap; mult["0" if r < .3 else "1" if r < 1.5 else "2" if r < 2.5 else "3+"] += L[u]
print("# tangle bp by depth class (x hap): " + "; ".join(f"{k}x {v/1e6:.2f} Mb" for k, v in sorted(mult.items())))
if a.out:
    with open(a.out + ".nodes.txt", "w") as o: o.write("\n".join(sorted(seen)) + "\n")
    with open(a.out + ".anchors.tsv", "w") as o:
        for v in hits: o.write(f"{v}\t{L[v]}\t{D.get(v,0):.1f}\t{','.join(sorted(hits[v]))}\t{where(v)}\n")
