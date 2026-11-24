#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Chains that run through a fork. build_chains.py joins along links that are one-to-one among kept unitigs, so a dropped alternative
can leave a chain running straight through a real fork. The usual case is harmless: a bubble near a chain end that rejoins the next chain of
the same group. The case that matters is a collapsed unitig, class depth about 2x, whose other exit leads through dropped unitigs to the end of
another chain that is on the homolog or in no group: both haplotypes pass through the hub, but the chain carries only one haplotype's arm, and the
other haplotype's route is lost. Its arm is then left unplaced, or the placed arm is reused by both haplotypes. Flagged HUB and listed first.
A chain in no group ("unassigned", or a bridge chain placed without one) never counts as the same group: two such chains are not a bubble.
break_at is the D04 break that separates the hub from the chain's arm on the fork side (1-based last base of the left part, for
edits/chain_breaks.txt). After the break, the hub is a chain end, so stage 06 can restore the dropped route as connectors.
Before ordering (stage 07) many forks inside collapsed chains lead to ungrouped chains, most of them ordinary bubbles, so HUB there needs a
homolog. With --placement and --chainClass (stage 10, after the pre-build) a row is flagged LOST, and listed first, when the hub is 2x, the
chain is placed (any copy) or on a tail candidate (--tails), and the other exit reaches a chain that is neither placed nor on a tail, at about 1x (median depth --oneLo to --oneHi) and at
least --minLost long: a haplotype's arm, lost because a chain ran through its fork. That is a D16 check (docs/decisions/D16_tails.md).
With --windows (unitig_depth.py on the chains), hub_total and hub_unique are the chain's own read depth over the hub unitig, all reads and
MAPQ >= 20. A true collapse reads about 2x in both; an allele bubble, whose unitig class depth is 2x because both alleles are near-identical,
reads about 1x on the chain (292: 4 of 7 LOST rows), and hub_unique falls when the arms' copies are identical elsewhere. So a row that would be
HUB or LOST is flagged BUBBLE instead when hub_total is under the 2x line (1.54 x --hap): the chain's own reads show it carries one allele, not
both haplotypes. On 292 (84x haploid), before any hub break: the four true collapses (chain0090, 0101, 0102, 0119) read 138-161x and stay
HUB, the 22 other rows 25-104x and become BUBBLE; at the rebuild's stage 10, 28 of 29 rows are BUBBLE, the one left (chain0122, 155x) a true
collapse. Without --windows, or
with no windows over the hub, the flag stands.
Usage: chain_forks.py chains.agp unitigs.gfa copy_number.tsv groups.tsv groups.chroms.tsv --hap H [--windows depth.windows.tsv]
       [--placement pre.placement.tsv --chainClass chain_class.tsv --oneLo L --oneHi H [--tails tails.tsv]] [--minLost 50000]'''
import argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "agp gfa cn groups chroms".split(): p.add_argument(x)
p.add_argument("--hap", type=float, required=True); p.add_argument("--hops", type=int, default=3)
p.add_argument("--placement"); p.add_argument("--chainClass"); p.add_argument("--oneLo", type=float); p.add_argument("--oneHi", type=float)
p.add_argument("--minLost", type=int, default=50000); p.add_argument("--tails"); p.add_argument("--windows"); a = p.parse_args()
if bool(a.placement) != bool(a.chainClass) or (a.placement and (a.oneLo is None or a.oneHi is None)):
    p.error("--placement, --chainClass, --oneLo and --oneHi go together")
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
PL = collections.defaultdict(set); CD = {}
if a.placement:
    import re
    for l in open(a.placement):
        if l[0] != "#": s = l.rstrip("\n").split("\t"); PL[re.sub(r"\.g\d+$", "", s[0])].add(s[4].startswith("placed"))
    for l in open(a.chainClass):
        if l[0] != "#": s = l.rstrip("\n").split("\t"); CD[s[0]] = (int(s[1]), float(s[2]))
TL = set()
if a.tails:
    for l in open(a.tails):
        if l[0] != "#": TL |= {x[:-1] for x in l.split("\t")[2].split(",")}
placed = lambda c: any(PL.get(c, ()))
used = lambda c: placed(c) or c in TL
def lost(c, t):
    if not a.placement or not used(c) or used(t) or t not in CD: return False
    n, d = CD[t]; return n >= a.minLost and a.oneLo <= d <= a.oneHi
WD = {}
if a.windows:
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from chainlib import read_windows
    WD = read_windows(a.windows)
def hubdepth(c, x, y):
    w = [r for r in WD.get(c, ()) if r[0] < y and r[1] > x - 1]
    if not w: return ".", "."
    med = lambda v: sorted(v)[len(v) // 2]
    return f"{med([r[2] for r in w]):.0f}", f"{med([r[3] for r in w]):.0f}"
out = []
for c, v in parts.items():
    if c in restored: continue
    for i, (u, st, x, y) in enumerate(v):
        for side, nxt, brk in (((u, st), v[i + 1][0] if i + 1 < len(v) else None, y), ((u, inv[st]), v[i - 1][0] if i else None, x - 1)):
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
                    grouped = lambda g: g in C
                    kind = ("homolog" if gc != gt and C[gc] == C[gt] else "same_group" if gc == gt else "other") if grouped(gc) and grouped(gt) else "ungrouped"
                    twox = S.get(u, (0, 0))[1] >= 1.54 * a.hap
                    flag = "LOST" if twox and lost(c, t) else "HUB" if twox and kind == "homolog" else "."
                    tl, td = CD.get(t, (".", ".")); hd = hubdepth(c, x, y)
                    if flag != "." and hd[0] != "." and float(hd[0]) < 1.54 * a.hap: flag = "BUBBLE"
                    out.append((flag, c, gc, x, y, u, S.get(u, (0, 0))[0], S.get(u, (0, 0))[1], q, t, gt, kind, brk, tl, td) + hd)
out.sort(key=lambda r: ({"LOST": 0, "HUB": 1, "BUBBLE": 2}.get(r[0], 3), -r[13] if isinstance(r[13], int) else 0))
print("#flag\tchain\tgroup\tstart\tend\tunitig\tprimary\tclass\tvia\tto_chain\tto_group\tkind\tbreak_at\tto_length\tto_depth\thub_total\thub_unique")
for r in out: print(*r, sep="\t")
