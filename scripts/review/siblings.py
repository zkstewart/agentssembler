#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''For low-depth reused chains, where do the read alignments go? Finds unitigs/chains sharing graph neighbours and reads.
Usage: siblings.py reused_depth.tsv chains.agp unitigs.gfa copy_number.tsv groups.tsv ont_vs_chains.paf.gz'''
import sys, gzip, collections
rd, agp, gfa, cn, grp, paf = sys.argv[1:7]
low = []
for l in open(rd):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[3] != "nan" and float(s[3]) < 40: low.append(s[0])
u2c = {}; parts = collections.defaultdict(list); clen = {}
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": u2c[s[5]] = s[0]; parts[s[0]].append(s[5]+s[8]); clen[s[0]] = int(s[2])
info = {}
hdr = None
for l in open(cn):
    s = l.rstrip("\n").split("\t")
    if hdr is None: hdr = s; continue
    info[s[0]] = dict(zip(hdr, s))
print("copy_number columns:", hdr)
gr = {}
for l in open(grp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); gr[s[0]] = s[2]
adj = collections.defaultdict(set); ulen = {}
for l in open(gfa):
    if l[0] == "S":
        s = l.split("\t"); ulen[s[1]] = len(s[2]) if s[2] != "*" else 0
    if l[0] == "L":
        s = l.split("\t"); adj[s[1]+s[2]].add(s[3]+s[4])
fl = {"+": "-", "-": "+"}
sib = {}
for c in low:
    print(f"\n== {c} parts {parts[c]}")
    for p in parts[c]:
        u, o = p[:-1], p[-1]
        nb = adj[u+"+"] | adj[u+"-"]
        out = set(adj[u+"+"]); inc = set(adj[u+"-"])
        cands = collections.Counter()
        for n in out:
            for m in adj[n[:-1]+fl[n[-1]]]:
                if m[:-1] != u: cands[m[:-1]] += 1
        for n in inc:
            for m in adj[n[:-1]+fl[n[-1]]]:
                if m[:-1] != u: cands[m[:-1]] += 1
        t = [f"{x}[{ulen.get(x)},{u2c.get(x,'REMOVED')},{gr.get(u2c.get(x,''),'-')}]" for x, _ in cands.most_common(6)]
        print(f" {p} len {ulen.get(u)}: siblings {t}")
