#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Acceptance checks for one TTT result (docs/TANGLES.md).
1. Every adjacency in each traversal must be a GFA L line (in either reading direction).
2. Multiplicities against depth: per node, TTT multiplicity vs our depth / hap_depth; disagreements listed.
3. Traversal as agentssembler chains (12_broken/chains.agp) with where each chain sits in the build (placement file given), and nodes used
   more than once (TTT reuse) vs nodes the build placed in both haplotypes.
Usage: check_traversal.py unitigs.gfa depth.tsv chains.agp placement.tsv outdir --hap HAP_DEPTH'''
import argparse, collections, re, os
p = argparse.ArgumentParser()
for x in "gfa depth agp placement outdir".split(): p.add_argument(x)
p.add_argument("--hap", type=float, required=True); a = p.parse_args()
L = set(); LEN = {}
for l in open(a.gfa):
    if l[0] == "L":
        s = l.split("\t"); x = (s[1], s[2], s[3], s[4]); L.add(x)
        f = {"+": "-", "-": "+"}; L.add((s[3], f[s[4]], s[1], f[s[2]]))
    elif l[0] == "S":
        s = l.split("\t", 3); m = re.search(r"LN:i:(\d+)", l); LEN[s[1]] = int(m.group(1)) if m else len(s[2])
D = {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(a.depth)}
U2C = {}
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.split("\t")
    if s[4] == "W": U2C[s[5]] = (s[0], s[8].strip())
PL = {}
for l in open(a.placement):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); PL[s[0]] = s[3] if s[4].startswith("placed") else s[4]
gaf = os.path.join(a.outdir, "traversal.gaf")
if not os.path.exists(gaf): raise SystemExit(f"no {gaf}")
use = collections.Counter(); ok = True
for l in open(gaf):
    s = l.rstrip("\n").split("\t"); name, path = s[0], next((x for x in s if x[:1] in "<>"), "")
    steps = re.findall(r"([<>])([^<>]+)", path)
    bad = []
    for (o1, u1), (o2, u2) in zip(steps, steps[1:]):
        if (u1, "+" if o1 == ">" else "-", u2, "+" if o2 == ">" else "-") not in L: bad.append(f"{o1}{u1}{o2}{u2}")
    for o, u in steps: use[u] += 1
    bp = sum(LEN.get(u, 0) for o, u in steps)
    print(f"== {name}: {len(steps)} nodes, {bp/1e6:.2f} Mb (overlaps not subtracted); adjacencies not in GFA: {len(bad)} {' '.join(bad[:5])}")
    if bad: ok = False
    ch = []
    for o, u in steps:
        c = U2C.get(u, ("-", "+"))[0]
        if not ch or ch[-1][0] != c: ch.append([c, 0])
        ch[-1][1] += LEN.get(u, 0)
    print("   chains: " + " ".join(f"{c}[{PL.get(c, 'dropped') if c != '-' else 'no_chain'}]" for c, n in ch if n >= 20000))
mfile = os.path.join(a.outdir, "traversal.multiplicities.csv")
if os.path.exists(mfile):
    dis = []
    for l in open(mfile):
        s = l.split()
        if len(s) < 3 or s[0] in ("node", "#node"): continue
        try: m = int(float(s[2]))
        except ValueError: continue
        u = s[0]; r = D.get(u, 0) / a.hap
        if abs(m - r) > 0.6 and LEN.get(u, 0) >= 20000: dis.append((u, LEN[u], m, r, use.get(u, 0)))
    print(f"multiplicity vs depth disagreements (>=20 kb, |m - depth/hap| > 0.6): {len(dis)}")
    for u, n, m, r, k in sorted(dis, key=lambda x: -x[1])[:15]: print(f"   {u}\t{n/1e3:.0f} kb\tTTT {m}\tdepth {r:.2f}x\tused {k}\t{U2C.get(u, ('-',))[0]}")
multi = [u for u, k in use.items() if k > 1 and LEN.get(u, 0) >= 20000]
print(f"nodes >= 20 kb used more than once: {len(multi)}: " + ", ".join(f"{u}({U2C.get(u, ('-',))[0]}, {D.get(u,0)/a.hap:.2f}x)" for u in multi[:20]))
print("RESULT: adjacencies " + ("all in GFA" if ok else "NOT all in GFA (reject)"))
