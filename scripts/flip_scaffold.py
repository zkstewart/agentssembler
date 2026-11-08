#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Reverse named scaffolds in a final AGP and FASTA, a change of writing direction only. Used when a direction flag in names.tsv stops
matching after re-ordering, with the pre-fill AGP left as it was. Usage: flip_scaffold.py in.agp in.fasta out.agp out.fasta name [name ...]'''
import sys
ia, fa, oa, ofa, *names = sys.argv[1:]; names = set(names); fl = {"+": "-", "-": "+"}
R = {}; out = []
for l in open(ia):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[0] in names: R.setdefault(s[0], []).append(s)
    else: out.append(s)
for n, rows in R.items():
    T = max(int(r[2]) for r in rows)
    for i, r in enumerate(reversed(rows), 1):
        a, b = int(r[1]), int(r[2]); r = r[:]; r[1], r[2] = str(T - b + 1), str(T - a + 1); r[3] = str(i)
        if r[4] == "W": r[8] = fl[r[8]]
        out.append(r)
out.sort(key=lambda r: (r[0], int(r[1])))
with open(oa, "w") as o:
    for r in out: o.write("\t".join(r) + "\n")
tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
with open(ofa, "w") as o:
    n = None; q = []
    def flush():
        if n is None: return
        x = "".join(q); o.write(f">{n}\n{x.translate(tr)[::-1] if n in names else x}\n")
    for l in open(fa):
        if l[0] == ">": flush(); n = l[1:].split()[0]; q = []
        else: q.append(l.strip())
    flush()
print("flipped", sorted(names))
