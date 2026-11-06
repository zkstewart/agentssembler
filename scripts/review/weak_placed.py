#!/usr/bin/env python3
'''Placed chains with no link to either order neighbour: pairs to their group and to the neighbours, and where the contacts lie.
Usage: weak_placed.py order.tsv end_links.tsv choices.tsv groups.tsv pairs.tsv.gz'''
import sys, gzip, collections
od, el, ch, gr, pa = sys.argv[1:6]
fl = {"+": "-", "-": "+"}
lk = set()
for f in (el, ch):
    for l in open(f):
        if l[0] != "#":
            s = l.split("\t"); lk.add((s[0], s[1])); lk.add((s[1], s[0]))
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); L[s[0]] = int(s[1]); g[s[0]] = s[2]
rows = collections.defaultdict(list)
for l in open(od):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); rows[s[0]].append(s)
Q = {}
for gg, rs in rows.items():
    for i, r in enumerate(rs):
        c, o = r[2], r[3]; left = right = False
        if i: p = rs[i - 1]; left = (p[2] + p[3], c + fl[o]) in lk
        if i + 1 < len(rs): n = rs[i + 1]; right = (c + o, n[2] + fl[n[3]]) in lk
        if not left and not right: Q[c] = (gg, i, rs)
cnt = {c: collections.Counter() for c in Q}
for l in gzip.open(pa, "rt"):
    s = l.split("\t"); a, b = s[0], s[2]
    if a == b: continue
    if a in cnt: cnt[a][b] += 1
    if b in cnt: cnt[b][a] += 1
print("#chain\tgroup\trank\tkb\tpairs_group\tpairs_other\tnbr_density_per_Mb\ttop_density_ranks")
for c, (gg, i, rs) in sorted(Q.items(), key=lambda kv: sum(cnt[kv[0]][r[2]] for r in kv[1][2])):
    order = [r[2] for r in rs]; pg = sum(cnt[c][x] for x in order if x != c); po = sum(cnt[c].values()) - pg
    nb = [order[j] for j in (i - 1, i + 1) if 0 <= j < len(order)]
    nd = sum(cnt[c][x] for x in nb) / max(sum(L[x] for x in nb) / 1e6, 1e-6)
    top = sorted(((cnt[c][x] / (L[x] / 1e6), k + 1) for k, x in enumerate(order) if x != c), reverse=True)[:3]
    print(f"{c}\t{gg}\t{i+1}/{len(order)}\t{L[c]//1000}\t{pg}\t{po}\t{nd:.0f}\t" + " ".join(f"#{k}:{d:.0f}" for d, k in top))
