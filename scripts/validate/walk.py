#!/usr/bin/env python3
'''Connected components of unplaced chains through graph links, with the placed chain ends they touch. Usage: walk.py placement.tsv end_links.tsv profile.tsv [minLen]'''
import sys, collections
pl, el, pr = sys.argv[1:4]
st = {}; nm = {}
for l in open(pl):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); st[s[0]] = s[4]; nm[s[0]] = s[3]
Pr = {}
for l in open(pr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); Pr[s[1]] = s
un = lambda c: st.get(c, "?").startswith("unplaced")
adj = collections.defaultdict(list)
for l in open(el):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); a, b = s[0], s[1]; o = int(s[4]); r = int(s[2])
    if o < 5000 and r < 3: continue
    adj[a[:-1]].append((a[-1], b, o, r)); adj[b[:-1]].append((b[-1], a, o, r))
seen = set(); comps = []
for c in st:
    if not un(c) or c in seen or ".g" in c: continue
    q = [c]; seen.add(c); m = [c]
    while q:
        x = q.pop()
        for e, y, o, r in adj[x]:
            t = y[:-1]
            if un(t) and t not in seen: seen.add(t); q.append(t); m.append(t)
    comps.append(m)
comps = [m for m in comps if len(m) > 0]
def touch(m):
    out = collections.Counter()
    for x in m:
        for e, y, o, r in adj[x]:
            t = y[:-1]
            if not un(t): out[f"{x}{e}~{y}[{nm.get(t,'?')}]"] = o
    return out
S = lambda m: sum(int(Pr[x][2]) for x in m if x in Pr)
for m in sorted(comps, key=lambda m: -sum(int(Pr[x][0]) for x in m if x in Pr))[:12]:
    g = sum(int(Pr[x][0]) for x in m if x in Pr)
    print(f"\ncomponent {len(m)} chains, {S(m)/1e6:.2f} Mb, {g} genes")
    for x in sorted(m, key=lambda x: -int(Pr[x][2]) if x in Pr else 0):
        p = Pr.get(x)
        if p: print(f"  {x} {int(p[2])//1000}k d={p[3]} rep={p[4]} {p[5]} pairs={p[6]} {p[7]} r={p[9]} g={p[0]} :: " + " ".join(f"{e}~{y}" for e, y, o, r in adj[x])[:150])
    t = touch(m)
    print("  touches placed:", " ".join(f"{k}:{v}" for k, v in t.most_common(20)))
