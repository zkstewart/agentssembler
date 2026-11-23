#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Juicebox assembly-mode input from pairs already lifted to scaffold coordinates. `juicer pre` works from a BAM on components, and a unitig
placed in two haplotypes cannot be told apart there, so each AGP sequence line becomes its own component and pairs are placed directly.
Writes prefix.assembly, prefix.components.tsv (component to scaffold, span, unitig, strand), prefix.txt (for juicer_tools pre) and prefix.sizes.
Juicebox stores positions in 32 bits, so an assembly over 2,147,483,647 bp is scaled: positions in prefix.txt and the size in prefix.sizes
are divided by the smallest whole factor that fits, written to prefix.scale (1 otherwise; .assembly keeps true lengths, as 3D-DNA and YaHS do).
A scaled map has not been tried here: check in Juicebox that scaffold edges fall where the .assembly puts them.
Usage: pairs_to_jbat.py final.agp pairs.scaffold.tsv.gz prefix
'''
import sys, gzip, bisect, collections
agp, pa, pre = sys.argv[1:4]
Sc = collections.OrderedDict(); n = 0
with open(pre + ".components.tsv", "w") as t:
    t.write("#component\tindex\tlength\tscaffold\tstart\tend\tunitig\tunitig_start\tunitig_end\tstrand\n")
    for l in open(agp):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t")
        if s[4] != "W": continue
        n += 1; L = int(s[2]) - int(s[1]) + 1; nm = f"ctg{n:08d}"
        Sc.setdefault(s[0], []).append((int(s[1]), int(s[2]), n, nm, L))
        t.write(f"{nm}\t{n}\t{L}\t{s[0]}\t{s[1]}\t{s[2]}\t{s[5]}\t{s[6]}\t{s[7]}\t{s[8]}\n")
off = {}; pos = 0
with open(pre + ".assembly", "w") as o:
    for sc, v in Sc.items():
        starts = []; 
        for a, b, i, nm, L in v:
            o.write(f">{nm} {i} {L}\n"); starts.append((a, b, pos)); pos += L
        off[sc] = starts
    for sc, v in Sc.items(): o.write(" ".join(str(x[2]) for x in v) + "\n")
MAXI = 2 ** 31 - 1; scale = -(-pos // MAXI)            # smallest whole factor that brings the assembly under the 32-bit limit
open(pre + ".sizes", "w").write(f"assembly {-(-pos // scale)}\n"); open(pre + ".scale", "w").write(f"{scale}\n")
st = {k: [x[0] for x in v] for k, v in off.items()}
def lift(c, p):
    v = off.get(c)
    if v is None: return None
    i = bisect.bisect_right(st[c], p) - 1
    if i < 0: return None
    a, b, base = v[i]
    return base + p - a + 1 if p <= b else None
k = 0; d = 0
with open(pre + ".txt", "w") as o:
    for l in gzip.open(pa, "rt"):
        a, x, b, y = l.rstrip("\n").split("\t"); u = lift(a, int(x)); w = lift(b, int(y))
        if u is None or w is None: d += 1; continue
        if u > w: u, w = w, u
        if scale > 1: u, w = -(-u // scale), -(-w // scale)
        o.write(f"0 assembly {u} 0 0 assembly {w} 1\n"); k += 1
print("components", n, "assembly bp", pos, "scale", scale, "pairs", k, "dropped", d)
