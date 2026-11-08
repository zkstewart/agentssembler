#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Do unplaced chains form connected tangles attached to one scaffold? Single repeat-rich chains have too little Hi-C to place, but a component's graph contacts and summed Hi-C can show
where the whole tangle belongs, such as a pericentromere missing from a short haplotype. Lists each component's size, the scaffolds it touches in the graph, and Hi-C per Mb of scaffold, also without
its largest Hi-C contributor, since one satellite chain can outweigh all the others. An indication for the report, not grounds to place anything: a component usually cannot be ordered.
Usage: unplaced_components.py placement.tsv chains.agp unitigs.gfa links.tsv [--minBp 300000]'''
import argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "placement agp gfa links".split(): p.add_argument(x)
p.add_argument("--minBp", type=int, default=300000); a = p.parse_args()
size, scaf, unp = {}, {}, set()
for l in open(a.placement):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); size[s[0]] = int(s[1])
    if s[4].startswith("placed"): scaf[s[0]] = s[3]
    elif s[4].startswith("unplaced"): unp.add(s[0])
U2C, C2U = collections.defaultdict(set), collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": U2C[s[5]].add(s[0]); C2U[s[0]].append(s[5])
adj = collections.defaultdict(set)
for l in open(a.gfa):
    if l[0] == "L":
        s = l.split("\t"); adj[s[1]].add(s[3]); adj[s[3]].add(s[1])
placed = {u for c in scaf for u in C2U[c]}
comp = {}
for u in adj:
    if u in comp or u in placed: continue
    st = [u]; comp[u] = u
    while st:
        for y in adj[st.pop()]:
            if y not in placed and y not in comp: comp[y] = u; st.append(y)
for c in unp:      # chains with no graph link are components of their own
    for u in C2U[c]: comp.setdefault(u, u)
M = collections.defaultdict(set)
for u, k in comp.items():
    for c in U2C[u]:
        if c in unp: M[k].add(c)
L = collections.Counter()
for l in open(a.links):
    s = l.split()
    if len(s) == 3 and s[1] != "." and s[2].isdigit(): L[s[0], s[1]] += int(s[2])
ssize = collections.Counter()
for c, t in scaf.items(): ssize[t] += size[c]
tot = sum(size[c] for c in unp)
print(f"# unplaced: {len(unp)} chains, {tot / 1e6:.2f} Mb")
def rank(h, n=4): return ", ".join(f"{t} {x} ({x / (ssize[t] / 1e6):.0f}/Mb)" for t, x in sorted(h.items(), key=lambda kv: -kv[1] / ssize[kv[0]])[:n]) or "none"
for k, v in sorted(M.items(), key=lambda kv: -sum(size[c] for c in kv[1])):
    bp = sum(size[c] for c in v)
    if bp < a.minBp: continue
    per = collections.defaultdict(collections.Counter)
    for (x, y), n in L.items():
        if x in v and y in scaf: per[x][scaf[y]] += n
        if y in v and x in scaf: per[y][scaf[x]] += n
    att = collections.Counter(t for u, kk in comp.items() if kk == k for y in adj[u] if y in placed for c in U2C[y] if c in scaf for t in [scaf[c]])
    allh = sum(per.values(), collections.Counter())
    top = max(per, key=lambda c: sum(per[c].values()), default=None)
    rest = sum((h for c, h in per.items() if c != top), collections.Counter())
    print(f"\n{len(v)} chains, {bp / 1e6:.2f} Mb; graph contacts: {dict(att.most_common(6)) or 'none'}")
    print(f"  Hi-C, all {sum(allh.values())} links: {rank(allh)}")
    if top: print(f"  Hi-C, without {top} ({sum(per[top].values())} links): {rank(rest)}")
    print("  chains: " + " ".join(sorted(v, key=lambda c: -size[c])))
