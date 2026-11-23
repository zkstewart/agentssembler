#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Hi-C gaps that swapping two chains between homologs would close. For each Hi-C junction (a gap), a flank E whose graph neighbour X sits in
the homolog's scaffold; Y is the chain across the gap from E. The swap puts X in Y's place and Y in X's: it is listed when X also links to Y's
other neighbour and Y links to both of X's neighbours, so both haplotypes would get a graph path where one now has a gap. Graph links are
end_links.tsv rows plus graph bridges whose whole path is connector chains (chain8*, chain9*). Chain-level: link ends are not checked.
The usual cause is the Hi-C allele rescue calling an allele pair against the graph (292: chr1_h1 17.67 Mb and chr3_h1 10.25 Mb, both
pairs called by Hi-C; the build without the rescue had no gap at either). Neither side is proof: in near-identical haplotypes a graph link
can cross between them, and Hi-C there is allele-ambiguous. A report for D10/D13 only; a swap is a group edit and a `reset 08_bridges`.
Scaffold names must end _hN (homologs share the rest).
Usage: rescue_swaps.py junctions.tsv end_links.tsv graph.bridges.tsv out.tsv [--rescued groups.rescued.tsv]'''
import argparse, collections, re
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "junctions links bridges out".split(): p.add_argument(x)
p.add_argument("--rescued"); a = p.parse_args()
rows = lambda f: (l.rstrip("\n").split("\t") for l in open(f) if l.strip() and l[0] != "#")
base = lambda c: re.sub(r"\.g\d+$", "", c)
conn = lambda c: re.match(r"chain[89]\d+", base(c)) is not None
R = {}
if a.rescued:
    for s in rows(a.rescued): R[s[0]] = f"{s[4]}:{s[6]}"
adj = collections.defaultdict(set)
for s in rows(a.links): adj[s[0][:-1]].add(s[1][:-1]); adj[s[1][:-1]].add(s[0][:-1])
for s in rows(a.bridges):
    if all(conn(x[:-1]) for x in s[3].split(",")): adj[s[0][:-1]].add(s[1][:-1]); adj[s[1][:-1]].add(s[0][:-1])
J = list(rows(a.junctions)); order = collections.defaultdict(list)
for s in J:
    o = order[s[0]]
    if not o: o.append(s[2][:-1])
    o.append(s[3][:-1])
where = {c: (sc, i) for sc, o in order.items() for i, c in enumerate(o)}
chrom = lambda sc: re.sub(r"_h\d+$", "", sc)
def nb(c):
    sc, i = where[c]; o = order[sc]; return (o[i - 1] if i else None, o[i + 1] if i + 1 < len(o) else None)
link = lambda x, y: x is None or y is None or y in adj[x]
out = []
for s in J:
    if s[4] != "hic": continue
    sc, pos, Lf, Rf = s[0], s[1], s[2][:-1], s[3][:-1]
    for E, Y in ((Lf, Rf), (Rf, Lf)):
        for X in sorted(adj[E]):
            if X not in where: continue
            xs = where[X][0]
            if xs == sc or chrom(xs) != chrom(sc): continue
            yl, yr = nb(Y); xl, xr = nb(X); yfar = yr if yl == E else yl
            if link(X, yfar) and link(Y, xl) and link(Y, xr):
                out.append((sc, pos, s[2], s[3], Y, R.get(base(Y), "-"), X, xs, R.get(base(X), "-")))
with open(a.out, "w") as o:
    o.write("#scaffold\tgap_at\tleft\tright\tswap_out\tswap_out_rescue\tswap_in\tswap_in_from\tswap_in_rescue\n")
    for r in out: o.write("\t".join(map(str, r)) + "\n")
print(f"rescue swaps: {len(out)} of {sum(s[4] == 'hic' for s in J)} Hi-C gaps")
