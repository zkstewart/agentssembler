#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Breadth-first walks in the graph from chain ends, through unassigned chains, until a placed chain end. Usage: trace.py end_links.tsv groups.tsv maxHops end ...'''
import sys, collections
el, gr, mh = sys.argv[1:4]; mh = int(mh)
fl = {"+": "-", "-": "+"}
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); g[s[0]] = s[2]; L[s[0]] = int(s[1])
adj = collections.defaultdict(set)
for l in open(el):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); adj[s[0]].add(s[1]); adj[s[1]].add(s[0])
for st in sys.argv[4:]:
    print(st, "len", L.get(st[:-1]))
    out = collections.Counter(); q = [(st, [st])]; seen = {}
    while q:
        e, p = q.pop(0)
        for n in adj[e]:
            c = n[:-1]
            if c in {z[:-1] for z in p}: continue
            pp = p + [n]
            if g.get(c, "unassigned") != "unassigned":
                out[(n, g[c], " ".join(z.replace("chain", "") for z in pp[1:-1]))] += 1
            elif len(pp) <= mh:
                if seen.get(n, 0) < 6:
                    seen[n] = seen.get(n, 0) + 1
                    q.append((c + fl[n[-1]], pp))
    for (n, gg, pa), k in sorted(out.items(), key=lambda x: x[0][2].count(" "))[:14]:
        print("   ->", n, gg, "via", pa or "direct")
