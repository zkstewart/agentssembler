#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Read depth flatness of pseudochromosomes: primary alignments (and MAPQ>=20) of competitively mapped reads, in
windows; windows with Ns are skipped. Usage: depth_flat.py paf.gz target.fai [window] [--names regex prefix]'''
import sys, gzip, re
import numpy as np
paf, fai = sys.argv[1:3]; W = int(sys.argv[3]) if len(sys.argv) > 3 else 50000
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
D = {c: np.zeros(L[c] // W + 1) for c in L if re.match(r".+_h\d+$", c)}
U = {c: np.zeros(L[c] // W + 1) for c in D}
un = 0; tot = 0
with gzip.open(paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 12)
        if "tp:A:P" not in s[12]: continue
        n = int(s[8]) - int(s[7]); tot += n
        if s[5] not in D: un += n; continue
        a, b = int(s[7]), int(s[8]); q = int(s[11])
        for w in range(a // W, b // W + 1):
            o = min(b, (w + 1) * W) - max(a, w * W)
            if o > 0:
                D[s[5]][w] += o / W
                if q >= 20: U[s[5]][w] += o / W
d = np.concatenate([x[:-1] for x in D.values()]); u = np.concatenate([x[:-1] for x in U.values()])
print(f"windows {len(d)}, Mb {len(d)*W/1e6:.1f}; bases aligned to unplaced {100*un/tot:.1f}%")
print(f"depth median {np.median(d):.1f}, MAPQ>=20 median {np.median(u):.1f}, CV {d.std()/d.mean():.3f}")
for lo, hi in ((0, 6), (6, 13.5), (13.5, 20), (20, 34), (34, 40.5), (40.5, 67.5), (67.5, 1e9)):
    print(f"  {lo:>5}-{hi:<9}: {100*((d>=lo)&(d<hi)).mean():5.2f}%")
print(f"  unique (MAPQ>=20) < 6x: {100*(u<6).mean():.2f}% of windows")
