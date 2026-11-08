#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Hi-C contact density from every placed chain to each group and to unassigned chains; flags chains that prefer another group.
Usage: group_profile.py groups.tsv pairs.tsv.gz chromosomes.tsv out.tsv'''
import sys, gzip, collections
gr, pa, cm, out = sys.argv[1:5]
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); g[s[0]] = s[2]; L[s[0]] = int(s[1])
c2 = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(cm) if l[0] != "#"}
sz = collections.Counter()
for c, x in g.items(): sz[x] += L[c]
C = collections.defaultdict(collections.Counter)
for l in gzip.open(pa, "rt"):
    s = l.split("\t"); a, b = s[0], s[2]
    if a == b: continue
    C[a][g.get(b, "?")] += 1; C[b][g.get(a, "?")] += 1
rows = []
for c, cnt in C.items():
    x = g.get(c, "?")
    if x in ("unassigned", "?"): continue
    dens = {k: v / (sz[k] / 1e6) for k, v in cnt.items() if k not in ("unassigned", "?")}
    if not dens: continue
    top = sorted(dens.items(), key=lambda t: -t[1])
    own = dens.get(x, 0); best, bd = top[0]; alt = next((t for t in top if t[0] != x), ("-", 0))
    rel = "own" if best == x else "homolog" if c2.get(best) == c2.get(x) else "OTHER_CHROM"
    rows.append((c, L[c], x, sum(cnt.values()), own, best, bd, alt[0], alt[1], rel, cnt["unassigned"]))
with open(out, "w") as o:
    o.write("#chain\tlength\tgroup\tpairs\town_density\tbest_group\tbest_density\tnext_group\tnext_density\trelation\tto_unassigned\n")
    for r in sorted(rows, key=lambda r: -r[1]):
        o.write("\t".join(f"{v:.1f}" if isinstance(v, float) else str(v) for v in r) + "\n")
c = collections.Counter(); b = collections.Counter()
for r in rows: c[r[9]] += 1; b[r[9]] += r[1]
print({k: (c[k], round(b[k] / 1e6, 1)) for k in c})
print("not own, >=50 kb:")
for r in sorted(rows, key=lambda r: -r[1]):
    if r[9] != "own" and r[1] >= 50000 and ".g" not in r[0]:
        print("\t".join(f"{v:.1f}" if isinstance(v, float) else str(v) for v in r))
print("own but next within 1.5x, >=300 kb, base chains:")
for r in sorted(rows, key=lambda r: -r[1]):
    if r[9] == "own" and r[1] >= 300000 and r[8] * 1.5 > r[4] and ".g" not in r[0]:
        print("\t".join(f"{v:.1f}" if isinstance(v, float) else str(v) for v in r))
