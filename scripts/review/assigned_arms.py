#!/usr/bin/env python3
'''Assigned non-shared chains under 500 kb linked at both ends to assigned chains that are not their order neighbours:
candidates for arms that were clustered into a group. Usage: assigned_arms.py order.tsv end_links.tsv groups.tsv'''
import sys, collections
od, el, gr = sys.argv[1:4]
fl = {"+": "-", "-": "+"}
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); g[s[0]] = s[2]; L[s[0]] = int(s[1])
adj = collections.defaultdict(set)
for l in open(el):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); adj[s[0]].add(s[1]); adj[s[1]].add(s[0])
pos = {}; rows = collections.defaultdict(list)
for l in open(od):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); pos[s[2]] = (s[0], int(s[1]), s[3]); rows[s[0]].append(s[2])
for c, (gg, r, o) in sorted(pos.items(), key=lambda kv: (kv[1][0], kv[1][1])):
    if ".g" in c or L[c] >= 500000: continue
    left = [x for x in adj[c + "-"]]; right = [x for x in adj[c + "+"]]
    def name(xs): return [f"{x}[{pos[x[:-1]][0]}#{pos[x[:-1]][1]}]" if x[:-1] in pos else f"{x}[{g.get(x[:-1])}]" for x in xs]
    lg = [x for x in left if x[:-1] in pos]; rg = [x for x in right if x[:-1] in pos]
    if not lg or not rg: continue
    prv = rows[gg][r-2] if r > 1 else None; nxt = rows[gg][r] if r < len(rows[gg]) else None
    fit = {x[:-1] for x in lg + rg} & {prv, nxt}
    if len(fit) < 2:
        print(f"{gg}#{r}\t{c}{o}\t{L[c]//1000}k\tprev {prv} next {nxt}\tL {name(left)}\tR {name(right)}")
