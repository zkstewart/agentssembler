#!/usr/bin/env python3
'''Chains that run through a fork. build_chains.py joins along links that are one-to-one among kept unitigs, so a dropped alternative
can leave a chain running straight through a real fork. The usual case is harmless: a bubble near a chain end that rejoins the next chain of
the same group. The case that matters is a collapsed unitig, class depth about 2x, whose other exit leads through dropped unitigs to the end of a chain on the homolog: both
haplotypes pass through it, but only one has it, and the other has a gap there. Lists those first.
Usage: chain_forks.py chains.agp unitigs.gfa copy_number.tsv groups.tsv groups.chroms.tsv --hap H'''
import argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "agp gfa cn groups chroms".split(): p.add_argument(x)
p.add_argument("--hap", type=float, required=True); p.add_argument("--hops", type=int, default=3); a = p.parse_args()
inv = {"+": "-", "-": "+"}
S = {}
for l in open(a.cn):
    if l[0] != "#": s = l.rstrip("\n").split("\t"); S[s[0]] = (float(s[2]), float(s[5]))
parts = collections.defaultdict(list); U = {}
for l in open(a.agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4] == "W": parts[s[0]].append((s[5], s[8], int(s[1]), int(s[2]))); U[s[5]] = s[0]
restored = {c for c in parts if c.startswith(("chain8", "chain9"))}
term = {}
for c, v in parts.items():
    if c not in restored: term[(v[-1][0], v[-1][1])] = c; term[(v[0][0], inv[v[0][1]])] = c
drop = lambda u: u not in U or U[u] in restored
adj = collections.defaultdict(set)
for l in open(a.gfa):
    if l[0] == "L":
        s = l.split("\t"); e1, e2 = (s[1], s[2]), (s[3], inv[s[4]]); adj[e1].add(e2); adj[e2].add(e1)
G = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(a.groups) if l[0] != "#"}
C = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(a.chroms) if l[0] != "#"}
out = []
for c, v in parts.items():
    if c in restored: continue
    for i, (u, st, x, y) in enumerate(v):
        for side, nxt in (((u, st), v[i + 1][0] if i + 1 < len(v) else None), ((u, inv[st]), v[i - 1][0] if i else None)):
            if nxt is None: continue
            for w, we in adj[side]:
                if w == nxt or not drop(w): continue
                stack = [((w, inv[we]), [w])]; hit = set()
                while stack:
                    e, q = stack.pop()
                    for f in adj[e]:
                        if f in term and term[f] != c: hit.add((term[f], ",".join(q)))
                        elif drop(f[0]) and f[0] not in q and len(q) < a.hops: stack.append(((f[0], inv[f[1]]), q + [f[0]]))
                for t, q in sorted(hit):
                    gc, gt = G.get(c, "."), G.get(t, ".")
                    kind = "homolog" if gc in C and gt in C and gc != gt and C[gc] == C[gt] else "same_group" if gc == gt else "other"
                    twox = S.get(u, (0, 0))[1] >= 1.54 * a.hap
                    out.append((kind == "homolog" and twox, c, gc, x, y, u, S.get(u, (0, 0))[0], S.get(u, (0, 0))[1], q, t, gt, kind))
out.sort(key=lambda r: not r[0])
print("#flag\tchain\tgroup\tstart\tend\tunitig\tprimary\tclass\tvia\tto_chain\tto_group\tkind")
for r in out: print("HUB" if r[0] else ".", *r[1:], sep="\t")
