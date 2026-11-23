#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Reverse read-derived joins that failed review, as in decision D11. For each row of decisions.tsv (scaffold, left unitig, right unitig, reason), the flanks are restored from
pre.agp and the gap comes back. Other joins in cand.agp stay. Every row needs a reason. Coordinates are renumbered.
Usage: apply_join_decisions.py cand.agp pre.agp decisions.tsv out.agp
'''
import sys, collections
cand, pre, dec, out = sys.argv[1:5]
def load(p):
    R = collections.OrderedDict()
    for l in open(p):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t"); R.setdefault(s[0], []).append(s)
    return R
C, P = load(cand), load(pre)
for l in open(dec):
    if l[0] == "#" or not l.strip(): continue
    sc, lc, rc, *why = l.rstrip("\n").split("\t")
    if not why or not why[0].strip(): sys.exit(f"{sc} {lc} {rc}: no reason given")
    p = P[sc]; gi = next((i for i in range(1, len(p) - 1) if p[i][4] in "NU" and {p[i - 1][5], p[i + 1][5]} == {lc, rc}), None)
    if gi is None: sys.exit(f"{sc}: no gap between {lc} and {rc} in {pre}")
    pl, pg, pr = p[gi - 1], p[gi], p[gi + 1]; v = C[sc]
    def find(a, b):
        return next((i for i in range(len(v) - 1) if v[i][5] == a[5] and (v[i + 1][5] == b[5] or (v[i + 1][5].startswith("rf_") and i + 2 < len(v) and v[i + 2][5] == b[5]))), None)
    idx = find(pl, pr); flipped = idx is None
    if flipped: idx = find(pr, pl)        # edits/flip.txt wrote the scaffold the other way round after pre.agp
    if idx is None: sys.exit(f"{sc}: join {lc} {rc} not found in {cand}")
    jr = idx + 1 if v[idx + 1][5] == (pl if flipped else pr)[5] else idx + 2
    L = v[idx][:]; R = v[jr][:]
    # restore only the unitig end facing this gap, so a closure at the other end of a component stays. The facing end is the same
    # unitig end in either writing direction: the right end of a "+" component on the left of the gap is its unitig's end (column 8).
    def restore(row, pre, left_of_gap):
        col = 7 if (pre[8] == "+") == left_of_gap else 6; row[col] = pre[col]
    if not flipped: restore(L, pl, True); restore(R, pr, False)
    else: restore(L, pr, False); restore(R, pl, True)
    g = pg[:]; g[0] = sc
    v[idx:jr + 1] = [L, g, R]
with open(out, "w") as o:
    o.write("##agp-version 2.1\n")
    for sc, v in C.items():
        x = 1
        for i, r in enumerate(v, 1):
            n = int(r[5]) if r[4] in "NU" else int(r[7]) - int(r[6]) + 1
            r = r[:]; r[1], r[2], r[3] = str(x), str(x + n - 1), str(i); x += n
            o.write("\t".join(r) + "\n")
