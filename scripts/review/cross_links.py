#!/usr/bin/env python3
'''Graph or read links from an end facing a gap to a specific chain of the homolog: possible phase switches. For each, Hi-C counts between the
chains of the two possible pairings. Usage: cross_links.py order.tsv junctions.tsv end_links.tsv chroms.tsv links.tsv'''
import sys, collections
od, jt, el, ch, hl = sys.argv[1:6]
fl = {"+": "-", "-": "+"}
chrom = dict(l.rstrip().split("\t")[:2] for l in open(ch))
G = {}; O = collections.defaultdict(list)
for l in open(od):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); G[s[2]] = s[0]; O[s[0]].append((s[2], s[3]))
nx = {}
for g, r in O.items():
    for u, v in zip(r, r[1:]):
        x, y = u[0] + u[1], v[0] + fl[v[1]]; nx[x] = y; nx[y] = x
H = collections.Counter()
for l in open(hl):
    a, b, n = l.rstrip().split("\t")
    if b != ".": H[frozenset((a, b))] = int(n)
for l in open(el):
    if l[0] == "#": continue
    x, y, r, q, o = l.rstrip().split("\t")
    a, b = x[:-1], y[:-1]
    if a not in G or b not in G or G[a] == G[b] or chrom[G[a]] != chrom[G[b]]: continue
    if ".g" in a and ".g" in b: continue
    print(f"{G[a]} {x} ~ {G[b]} {y} reads={r} mq={q} overlap={o}")
    for e, f in ((x, y), (y, x)):
        n = nx.get(e)
        if n: print(f"   {e} currently next to {n}; Hi-C {e[:-1]}-{n[:-1]} {H[frozenset((e[:-1], n[:-1]))]}, {e[:-1]}-{f[:-1]} {H[frozenset((e[:-1], f[:-1]))]}")
        else: print(f"   {e} is a scaffold end; Hi-C to {f[:-1]} {H[frozenset((e[:-1], f[:-1]))]}")
