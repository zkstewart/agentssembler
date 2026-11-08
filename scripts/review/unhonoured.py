#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Order tests needing no reference. (1) Same-group graph links and bridges an order does not honour. (2) Shared-chain order and
orientation compared between the two homologs of a chromosome. (3) Telomere ends facing inward, and telomere ends per group.
Usage: unhonoured.py order.tsv end_links.tsv graph.bridges.tsv telomeres.tsv chroms.tsv'''
import sys, collections
od, el, br, te, ch = sys.argv[1:6]
fl = {"+": "-", "-": "+"}
O = collections.defaultdict(list)
for l in open(od):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[1] != ".": O[s[0]].append((s[2], s[3], int(s[4])))
grp = {c: g for g, r in O.items() for c, _, _ in r}
next_to = {}
for g, r in O.items():
    for a, b in zip(r, r[1:]):
        x, y = a[0] + a[1], b[0] + fl[b[1]]; next_to[x] = y; next_to[y] = x
part = collections.defaultdict(dict)
for l in open(el):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); part[s[0]][s[1]] = "direct"; part[s[1]][s[0]] = "direct"
for l in open(br):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); part[s[0]].setdefault(s[1], "bridge"); part[s[1]].setdefault(s[0], "bridge")
rank = {c: i for g, r in O.items() for i, (c, _, _) in enumerate(r)}
print("== ends with same-group graph partners, none of them adjacent in the order")
n = 0
for e in sorted(part):
    c = e[:-1]
    if c not in grp: continue
    same = {p: k for p, k in part[e].items() if p[:-1] in grp and grp[p[:-1]] == grp[c] and p[:-1] != c}
    if same and next_to.get(e) not in same:
        n += 1
        print(grp[c], rank[c], e, "next:", next_to.get(e, "scaffold_end"), "partners:", " ".join(f"{p}@{rank[p[:-1]]}({k})" for p, k in same.items()))
print("count", n)
print("\n== shared chains: order/orientation between homologs")
pairs = collections.defaultdict(set)
for l in open(ch):
    if l[0] != "#":
        s = l.rstrip().split("\t"); pairs[s[1]].add(s[0])
gs = sorted(O)
seen = set()
for a in gs:
    for b in gs:
        if a >= b: continue
        ca = [(c.split(".g")[0], o) for c, o, _ in O[a] if ".g" in c]
        cb = [(c.split(".g")[0], o) for c, o, _ in O[b] if ".g" in c]
        A = {c for c, _ in ca}; B = {c for c, _ in cb}
        if not (A & B): continue
        if len(A & B) < 2: print(a, b, "few shared", A & B); continue
        ra = [x for x in ca if x[0] in B]; rb = [x for x in cb if x[0] in A]
        if ra == rb: print(a, b, len(ra), "same order and orientation"); continue
        rr = [(c, fl[o]) for c, o in reversed(rb)]
        if ra == rr: print(a, b, len(ra), "same, opposite direction"); continue
        print(a, b, len(ra), "DIFFER")
        print("  ", " ".join(c.replace("chain", "") + o for c, o in ra)); print("  ", " ".join(c.replace("chain", "") + o for c, o in rb))
print("\n== telomeres")
tel = {}
for l in open(te):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); tel[s[0]] = s[4]
for g, r in sorted(O.items()):
    t = []
    for i, (c, o, _) in enumerate(r):
        v = tel.get(c, "--")
        if v == "--": continue
        left, right = v[0] == "T", v[1] == "T"
        if o == "-": left, right = right, left
        t.append(f"{c}@{i}/{len(r)-1}{'<' if left else ''}{'>' if right else ''}")
    print(g, " ".join(t))
