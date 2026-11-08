#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Each adjacent pair with no link/bridge: search the whole chain-end graph at base-chain level up to --hops chains,
through any chain, for the ordered ends and for other orientations. Reports shortest path and what groups it passes.
Usage: none_paths.py order.tsv end_links.tsv graph.bridges.tsv groups.tsv [hops]'''
import sys, collections
od, el, br, gr = sys.argv[1:5]; H = int(sys.argv[5]) if len(sys.argv) > 5 else 8
fl = {"+": "-", "-": "+"}; base = lambda c: c.split(".g")[0]
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); L[s[0]] = int(s[1]); g[s[0]] = s[2]
adj = collections.defaultdict(set); direct = set()
for l in open(el):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); direct.add((s[0], s[1])); direct.add((s[1], s[0]))
    a, b = (base(s[0][:-1]), s[0][-1]), (base(s[1][:-1]), s[1][-1])
    adj[a].add(b); adj[b].add(a)
bridge = set()
for l in open(br):
    if l[0] == "#": continue
    s = l.split("\t"); bridge.add((s[0], s[1])); bridge.add((s[1], s[0]))
def bfs(a, b):
    '''a: end leaving, b: end entering; returns shortest list of intervening chains'''
    q = collections.deque([(a, [])]); seen = {a}
    while q:
        cur, path = q.popleft()
        for nx in adj[cur]:
            if nx == b: return path
            if len(path) >= H or nx[0] == a[0] or nx[0] == b[0]: continue
            e = (nx[0], fl[nx[1]])
            if e in seen: continue
            seen.add(e); q.append((e, path + [nx[0] + ("+" if nx[1] == "-" else "-")]))
    return None
O = collections.defaultdict(list)
for l in open(od):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); O[s[0]].append(s)
placed = {}
for gg, rows in O.items():
    for r in rows: placed[(gg, base(r[2]))] = int(r[1])
gb = {}
for c, x in g.items(): gb.setdefault(base(c), set()).add(x)
cat = collections.Counter()
for gg, rows in O.items():
    for r, d in zip(rows, rows[1:]):
        a = r[2] + r[3]; b = d[2] + fl[d[3]]
        if (a, b) in direct or (a, b) in bridge: continue
        ae = (base(r[2]), r[3]); be = (base(d[2]), fl[d[3]])
        p = bfs(ae, be)
        alt = []
        for oa in "+-":
            for ob in "+-":
                if (oa, ob) == (r[3], fl[d[3]]): continue
                q = bfs((ae[0], oa), (be[0], ob))
                if q is not None and len(q) <= 3: alt.append(f"{oa}{ob}:{len(q)}")
        if p is None: kind = "no_path"; txt = "."
        else:
            gs = [sorted(gb.get(z[:-1], {"?"})) for z in p]
            other = [z for z, s in zip(p, gs) if s != ["unassigned"]]
            kind = "long_path" if not other else "via_assigned"
            txt = ",".join(f"{z}[{'/'.join(s)}:{L.get(z[:-1], L.get(z[:-1]+'.'+s[0], 0))//1000}k]" for z, s in zip(p, gs))
        cat[kind] += 1
        print(f"{gg}\t{r[1]}\t{r[2]}{r[3]}\t{d[2]}{d[3]}\t{kind}\t{len(p) if p is not None else '-'}\t{txt}\t{' '.join(alt)}")
print("#", dict(cat))
