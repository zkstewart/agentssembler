#!/usr/bin/env python3
'''For each adjacent pair in the order with no bridge candidate, walk the graph from both ends through unassigned chains. Reports
the sequence each end can reach, total bp and depth, and any placed chain end reached. Usage: gap_walks.py order.tsv end_links.tsv groups.tsv graph.bridges.tsv depth.tsv maxHops'''
import sys, collections
od, el, gr, br, dp, mh = sys.argv[1:7]; mh = int(mh)
fl = {"+": "-", "-": "+"}
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); g[s[0]] = s[2]; L[s[0]] = int(s[1])
D = {}
for l in open(dp):
    s = l.rstrip("\n").split("\t")
    if s[0][:5] == "chain":
        try: D[s[0]] = float(s[2])
        except ValueError: pass
adj = collections.defaultdict(set)
for l in open(el):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); adj[s[0]].add(s[1]); adj[s[1]].add(s[0])
has = set()
for l in open(br):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[2] == "yes": has.add((s[0], s[1])); has.add((s[1], s[0]))
O = collections.defaultdict(list)
for l in open(od):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[1] != ".": O[s[0]].append((s[2], s[3]))
def walk(st):
    seen = set(); q = [(st, 0)]; hit = set(); un = set()
    while q:
        e, h = q.pop(0)
        for n in adj[e]:
            c = n[:-1]
            if g.get(c, "unassigned") != "unassigned":
                if c != st[:-1]: hit.add(n + ":" + g[c])
            elif n not in seen and h < mh:
                seen.add(n); un.add(c); q.append((c + fl[n[-1]], h + 1))
    return un, hit
def desc(st):
    un, hit = walk(st)
    bp = sum(L.get(c, 0) for c in un)
    big = sorted(un, key=lambda c: -L.get(c, 0))[:5]
    return f"{st}: {len(un)} chains {bp/1e3:.0f} kb [" + " ".join(f"{c[5:]}:{L.get(c,0)//1000}k/{D.get(c,float('nan')):.0f}x" for c in big) + "] placed:" + (",".join(sorted(hit)) or "-")
for gg, r in sorted(O.items()):
    for (a, ao), (b, bo) in zip(r, r[1:]):
        x, y = a + ao, b + fl[bo]
        if y in adj[x] or (x, y) in has: continue
        print(gg, desc(x)); print("    ", desc(y))
    print(gg, "START", desc(r[0][0] + fl[r[0][1]])); print(gg, "END  ", desc(r[-1][0] + r[-1][1]))
