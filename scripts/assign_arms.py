#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Collapsed samples: prepare bridges for ranking. chain_bridges.py finds paths per haplotype copy, so a bubble between shared
chains appears once per haplotype. This collapses copy-level bridges to base chains (same junction, same paths) for
choose_bridges.py to rank once. assign_arms2.py then makes one choice per copy-level junction.

Usage: assign_arms.py prepare graph.bridges.tsv base.bridges.tsv
'''
import argparse
fl = {"+": "-", "-": "+"}
base = lambda c: c.split(".g")[0]
bend = lambda e: base(e[:-1]) + e[-1]
def rev(path): return [z[:-1] + fl[z[-1]] for z in reversed(path)]
def norm(x, y, path):
    '''junction key and path written from the lower base end'''
    bx, by = bend(x), bend(y)
    return (bx, by, path) if bx <= by else (by, bx, rev(path))
def read_bridges(f):
    for l in open(f):
        if l[0] == "#": continue
        r = l.rstrip("\n").split("\t")
        if r[2] == "yes": yield r[0], r[1], r[3].split(","), int(r[4]), r[5]
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("mode", choices=["prepare"]); ap.add_argument("files", nargs=2)
a = ap.parse_args()
seen = {}
for x, y, path, bp, mo in read_bridges(a.files[0]):
    bx, by, pth = norm(x, y, path); seen[(bx, by, tuple(pth))] = (bp, mo)
with open(a.files[1], "w") as o:
    o.write("#end_a\tend_b\tsame_group\tpath\tpath_bp\tmin_overlap\n")
    for (bx, by, pth), (bp, mo) in sorted(seen.items()):
        o.write(f"{bx}\t{by}\tyes\t{','.join(pth)}\t{bp}\t{mo}\n")
print(len(seen), "base candidate paths,", len({(k[0], k[1]) for k in seen}), "base junctions")
