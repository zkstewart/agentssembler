#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Evidence for unplaced chains: gene count (compleasm), depth, repeat, Hi-C by chromosome and group (density per Mb of partners), top partner chains, and graph
links to placed chains. Usage: unplaced_profile.py placement.tsv class.tsv full_table.tsv groups.tsv chroms.tsv pairs.tsv.gz end_links.tsv out.tsv'''
import sys, gzip, collections
pl, cl, ft, gr, ch, pa, el, out = sys.argv[1:9]
Un = {}
for l in open(pl):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4].startswith("unplaced") and ".g" not in s[0]: Un[s[0]] = int(s[1])
C = {}
for l in open(cl):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); C[s[0]] = (float(s[2]), float(s[5]), s[6])
G = collections.Counter()
for l in open(ft):
    if l[0] == "#" or l.startswith("Gene"): continue
    s = l.split("\t")
    if s[1] in ("Single", "Duplicated"): G[s[2].replace("unplaced_", "")] += 1
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); g[s[0]] = s[2]; L[s[0]] = int(s[1])
g2c = {}
for l in open(ch):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); g2c[s[0]] = s[1]
sz = collections.Counter(); cz = collections.Counter()
for c, x in g.items():
    if x in g2c: sz[x] += L[c]; cz[g2c[x]] += L[c]
P = {c: collections.Counter() for c in Un}
for l in gzip.open(pa, "rt"):
    s = l.split("\t"); a, b = s[0], s[2]
    if a == b: continue
    if a in P: P[a][b] += 1
    if b in P: P[b][a] += 1
lk = collections.defaultdict(list)
for l in open(el):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    for x, y in ((s[0], s[1]), (s[1], s[0])):
        if x[:-1] in Un: lk[x[:-1]].append(f"{x[-1]}~{y}[{g.get(y[:-1],'?')}]:{s[4]}")
rows = []
for c, n in Un.items():
    bych = collections.Counter(); byg = collections.Counter(); sharedp = 0
    for t, v in P[c].items():
        x = g.get(t)
        if x in g2c: bych[g2c[x]] += v; byg[x] += v
    dch = sorted(((v / (cz[k] / 1e6), k) for k, v in bych.items()), reverse=True)
    dg = sorted(((v / (sz[k] / 1e6), k) for k, v in byg.items()), reverse=True)
    r = dch[0][0] / max(dch[1][0], .01) if len(dch) > 1 else 0
    gs = [x for d, x in dg if g2c[x] == dch[0][1]][:2] if dch else []
    gr2 = dict((x, d) for d, x in dg)
    d, rep, kl = C.get(c, (0, 0, "?"))
    top = " ".join(f"{k}:{v}" for k, v in P[c].most_common(4))
    rows.append((G[c], c, n, d, rep, kl, sum(P[c].values()), dch[0][1] if dch else "-", round(dch[0][0], 1) if dch else 0, round(r, 2), " ".join(f"{x}:{gr2[x]:.1f}" for x in gs), top, ";".join(lk[c][:4])))
rows.sort(key=lambda r: (-r[0], -r[2]))
with open(out, "w") as o:
    o.write("#genes\tchain\tlength\tdepth\trepeat\tclass\tpairs\ttop_chrom\tdensity\tchrom_ratio\tgroups\ttop_partners\tlinks\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print("unplaced chains", len(rows), "with genes", sum(1 for r in rows if r[0]), "genes", sum(r[0] for r in rows))
print("with genes, ratio>=2:", sum(1 for r in rows if r[0] and r[9] >= 2), "; all, ratio>=2 and pairs>=200 and >=50 kb:", sum(1 for r in rows if r[9] >= 2 and r[6] >= 200 and r[2] >= 50000), round(sum(r[2] for r in rows if r[9] >= 2 and r[6] >= 200 and r[2] >= 50000) / 1e6, 2), "Mb")
