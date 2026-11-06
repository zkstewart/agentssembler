#!/usr/bin/env python3
'''
Paths through short unassigned chains that connect the ends of two assigned chains, from chain-end links.
A path enters a short chain through one end and must leave through the other. Search up to --hops short chains.
Each bridge: end_a end_b same_group path(chain/orientation, read order) path_bp min_link_overlap.
Also writes tips: short chains linked to an assigned chain on one side only.
Usage: chain_bridges.py end_links.tsv groups.tsv out_prefix [--hops 3 --minOverlap 5000 --maxShort 1000000]
'''
import argparse, collections
p = argparse.ArgumentParser(); p.add_argument("links"); p.add_argument("groups"); p.add_argument("out")
p.add_argument("--hops", type=int, default=3); p.add_argument("--minOverlap", type=int, default=5000); p.add_argument("--minReads", type=int, default=3)
a = p.parse_args()
g = {}; L = {}
for l in open(a.groups):
    if l[0] != "#":
        c, n, x = l.rstrip().split("\t"); g[c] = x; L[c] = int(n)
adj = collections.defaultdict(list)
for l in open(a.links):
    if l[0] == "#": continue
    x, y, r, q, ov = l.rstrip().split("\t"); r = int(r); ov = int(ov)
    if ov < a.minOverlap and r < a.minReads: continue
    e1 = (x[:-1], x[-1]); e2 = (y[:-1], y[-1])
    adj[e1].append((e2, ov, r)); adj[e2].append((e1, ov, r))
fl = {"+": "-", "-": "+"}
assigned = lambda c: g[c] != "unassigned"
br = {}; tips = collections.defaultdict(list)
for e in list(adj):
    if not assigned(e[0]): continue
    # walk outward from an assigned chain's end
    st = [(e, [], 10**9)]
    while st:
        cur, path, mo = st.pop()
        for nx, ov, r in adj[cur]:
            c, en = nx
            if any(c == z[0] for z in path) or c == e[0]: continue
            m = min(mo, ov if ov > 0 else 10**9)
            if assigned(c):
                if path:
                    k = (e, nx, tuple(path))
                    rk = (nx, e, tuple((z[0], fl[z[1]]) for z in reversed(path)))
                    if rk not in br: br[k] = m
                continue
            ori = "+" if en == "-" else "-"    # entered through 5' means forward
            np = path + [(c, ori)]
            if len(np) <= a.hops:
                st.append(((c, fl[en]), np, m))
                if not any(o for o in adj[(c, fl[en])] if o[0][0] != e[0]): tips[e].append((tuple(np), m))
with open(a.out + ".bridges.tsv", "w") as o:
    o.write("#end_a\tend_b\tsame_group\tpath\tpath_bp\tmin_overlap\n")
    for (e, nx, path), m in sorted(br.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        bp = sum(L[c] for c, _ in path)
        o.write(f"{e[0]}{e[1]}\t{nx[0]}{nx[1]}\t{'yes' if g[e[0]] == g[nx[0]] else 'no'}\t{','.join(c + s for c, s in path)}\t{bp}\t{m}\n")
with open(a.out + ".tips.tsv", "w") as o:
    o.write("#end\tgroup\tpath\tbp\tmin_overlap\n")
    for e, v in sorted(tips.items()):
        for path, m in v:
            o.write(f"{e[0]}{e[1]}\t{g[e[0]]}\t{','.join(c + s for c, s in path)}\t{sum(L[c] for c, _ in path)}\t{m}\n")
print(f"{len(br)} bridges, same-group {sum(1 for (e, nx, _) in br if g[e[0]] == g[nx[0]])}; ends with tips {len(tips)}")
