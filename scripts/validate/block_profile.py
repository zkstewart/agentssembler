#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Contacts between a flank and successive windows through a block, normalised by each window's total pairs, so a window that is repeat-rich does not
look like distant. Layout right: signal falls with distance from that flank. Usage: block_profile.py pairs.tsv.gz scaffold:x:y [--flank 1000000 --win 250000]'''
import sys, gzip, collections
pa = sys.argv[1]; B = []
for r in sys.argv[2:]:
    c, x, y = r.split(":"); B.append((c, int(x), int(y)))
F = 1000000; W = 250000
tot = collections.Counter(); L = collections.Counter(); R = collections.Counter()
want = {b[0] for b in B}
for l in gzip.open(pa, "rt"):
    a, x, b, y = l.rstrip("\n").split("\t")
    if a != b or a not in want: continue
    x = int(x); y = int(y)
    for c, s, e in B:
        if c != a: continue
        for p, q in ((x, y), (y, x)):
            if s <= p < e:
                w = (p - s) // W; tot[(c, s, w)] += 1
                if s - F <= q < s: L[(c, s, w)] += 1
                if e <= q < e + F: R[(c, s, w)] += 1
for c, s, e in B:
    n = (e - s) // W + 1
    print(f"{c} {s/1e6:.2f}-{e/1e6:.2f}")
    print(" window  total  left_flank  right_flank  left%  right%")
    for w in range(n):
        t = tot[(c, s, w)]
        print(f" {(s+w*W)/1e6:6.2f} {t:6d} {L[(c,s,w)]:6d} {R[(c,s,w)]:6d} {100*L[(c,s,w)]/max(t,1):5.1f} {100*R[(c,s,w)]/max(t,1):5.1f}")
