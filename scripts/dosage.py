#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Final check against the reference: over windows of the reference, how many scaffolds of the assembly carry that
sequence? Rearranged haplotypes should still sum to the ploidy if the translocations are balanced.
Reports the proportion of windows at each copy number and the carriers of windows other than --ploidy.
Usage: dosage.py final.agp unitigs_vs_ref.paf.gz ref.fai --ploidy N [--bin 500000 --minFrac 0.3]
'''
import argparse, collections, gzip
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("agp"); p.add_argument("paf"); p.add_argument("fai"); p.add_argument("--ploidy", type=int, required=True); p.add_argument("--refPrefix", default="", help="use only reference sequences with this prefix"); p.add_argument("--bin", type=int, default=500000)
p.add_argument("--minFrac", type=float, default=0.3); a = p.parse_args()
U = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": U[s[5]].append((int(s[6]) - 1, int(s[7]), s[0]))
C = collections.defaultdict(collections.Counter)
with gzip.open(a.paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 12)
        if "tp:A:P" not in s[12] or int(s[11]) < 20 or s[0] not in U: continue
        q0, q1, t0, t1 = int(s[2]), int(s[3]), int(s[7]), int(s[8])
        for u0, u1, sc in U[s[0]]:
            x0, x1 = max(q0, u0), min(q1, u1)
            if x1 - x0 < 5000: continue
            # scale the placed part of the query onto the target interval
            f0 = (x0 - q0) / (q1 - q0); f1 = (x1 - q0) / (q1 - q0)
            if s[4] == "-": f0, f1 = 1 - f1, 1 - f0
            y0 = int(t0 + f0 * (t1 - t0)); y1 = int(t0 + f1 * (t1 - t0))
            for b in range(y0 // a.bin, y1 // a.bin + 1):
                n = min(y1, (b + 1) * a.bin) - max(y0, b * a.bin)
                if n > 0: C[(s[5], b)][sc] += n
cn = collections.Counter(); runs = []
for l in open(a.fai):
    r, n = l.split("\t")[:2]; n = int(n)
    if not r.startswith(a.refPrefix): continue
    for b in range(n // a.bin):
        k = {s for s, v in C[(r, b)].items() if v >= a.minFrac * a.bin}
        cn[len(k)] += 1; lab = ",".join(sorted(k))
        if len(k) != a.ploidy:
            if runs and runs[-1][0] == r and runs[-1][1] == len(k) and runs[-1][3] == b: runs[-1][3] = b + 1; runs[-1][4].update(k)
            else: runs.append([r, len(k), b, b + 1, set(k)])
T = sum(cn.values())
print("copies\twindows\tfraction"); [print(k, v, f"{v/T:.3f}", sep="\t") for k, v in sorted(cn.items())]
print(f"runs of >= 2 Mb away from {a.ploidy} copies:")
for r, c, x, y, k in runs:
    if (y - x) * a.bin >= 2000000: print(f"{r}\t{x*a.bin/1e6:.1f}-{y*a.bin/1e6:.1f}\tcopies={c}\t{','.join(sorted(k))}")
