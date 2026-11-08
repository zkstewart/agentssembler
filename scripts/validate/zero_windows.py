#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Reference windows at 0 copies in the final AGP: what aligns there, from all unitigs, and is it unplaced, removed or placed?
Usage: zero_windows.py final.agp unplaced.fai chains.agp unitigs_vs_ref.paf.gz ref.fai placement.tsv'''
import sys, gzip, collections
agp, unp, cagp, paf, fai, pl = sys.argv[1:7]; B = 500000
placed = set()
for l in open(agp):
    if l[0] != "#":
        s = l.split("\t")
        if s[4] == "W": placed.add(s[5])
u2c = collections.defaultdict(set)
for l in open(cagp):
    if l[0] != "#":
        s = l.split("\t")
        if s[4] == "W": u2c[s[5]].add(s[0].split(".g")[0])
st = {}
for l in open(pl):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); st[s[0].split(".g")[0]] = s[4]
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
cov = collections.defaultdict(lambda: collections.Counter())
for l in gzip.open(paf, "rt"):
    s = l.split("\t", 12)
    if "tp:A:P" not in l or s[5] not in L: continue
    u = s[0]; a, b = int(s[7]), int(s[8])
    if u in placed: k = "placed"
    elif u in u2c: k = "unplaced:" + ",".join(sorted({st.get(c, "?").split(":")[-1] for c in u2c[u]}))
    else: k = "not_in_chains"
    for w in range(a // B, b // B + 1):
        o = min(b, (w + 1) * B) - max(a, w * B)
        if o > 0: cov[(s[5], w)][k] += o
rows = []
tot = collections.Counter(); W = 0; none = 0
for c in sorted(L, key=lambda x: int(x[3:]) if x[3:].isdigit() else 99):
    for w in range(L[c] // B + 1):
        v = cov[(c, w)]; rows.append((c, w, v))
import statistics
# windows where placed < 0.3 of a window per copy equal to dosage 0-ish
z = [(c, w, v) for c, w, v in rows if v["placed"] < 0.3 * B]
print("windows with placed coverage below 0.3:", len(z))
for c, w, v in z:
    W += 1
    for k, x in v.items(): tot[k] += x
    if sum(v.values()) < 0.1 * B: none += 1
print("bp by source (Mb):", {k: round(x / 1e6, 1) for k, x in tot.most_common()}); print("almost nothing aligned there:", none)
run = collections.defaultdict(list)
for c, w, v in z: run[c].append(w)
for c, ws in run.items():
    out = []; st0 = ws[0]; pr = ws[0]
    for w in ws[1:] + [None]:
        if w != pr + 1 if w is not None else True:
            out.append(f"{st0*B/1e6:.1f}-{(pr+1)*B/1e6:.1f}"); st0 = w
        pr = w
    print(c, " ".join(out))
