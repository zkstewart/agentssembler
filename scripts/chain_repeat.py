#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Repetitiveness of each chain from its own sequence: fraction of 31-mer positions whose k-mer occurs
at least twice in the same chain, and the top tandem period (most common distance between consecutive
occurrences of a repeated k-mer).
Usage: chain_repeat.py chains.fasta out.tsv [minLen]
'''
import sys, collections
fa, out = sys.argv[1:3]; minLen = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
K = 31
def chains():
    n = None; b = []
    for l in open(fa):
        if l[0] == ">":
            if n: yield n, "".join(b)
            n = l[1:].split()[0]; b = []
        else: b.append(l.strip().upper())
    if n: yield n, "".join(b)
with open(out, "w") as o:
    o.write("#chain\tlength\tfrac_repeated_31mers\ttop_period\ttop_period_share\n")
    for n, s in chains():
        if len(s) < minLen: continue
        last = {}; rep = 0; d = collections.Counter(); tot = 0
        for i in range(0, len(s) - K, 1):
            k = hash(s[i:i+K]); tot += 1
            j = last.get(k)
            if j is not None:
                rep += 1; d[i - j] += 1
            last[k] = i
        top = d.most_common(1)[0] if d else (0, 0)
        o.write(f"{n}\t{len(s)}\t{rep/max(tot,1):.4f}\t{top[0]}\t{top[1]/max(rep,1):.3f}\n")
        o.flush()
