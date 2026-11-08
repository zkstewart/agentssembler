#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Mapping targets for a diploid assembly in which shared sequence appears in both haplotypes. A read over such sequence would
otherwise get MAPQ 0 and could not anchor. Each target keeps one haplotype whole and hard-masks, in the others, every AGP
span of a unitig that also occurs in the kept haplotype. Haplotype-specific sequence stays in both, so it still competes.

Usage: mask_shared.py scaffolds.agp target.fasta out_prefix    (writes prefix.hN.fasta per haplotype, prefix.masked.tsv)
'''
import re, sys, collections
agp, fa, pre = sys.argv[1:4]
hap = lambda s: s.rsplit("_", 1)[1] if re.search(r"_h\d+$", s) else None
rows = []; where = collections.defaultdict(set)
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": rows.append((s[0], int(s[1]), int(s[2]), s[5])); where[s[5]].add(s[0])
haps = sorted({hap(r[0]) for r in rows if hap(r[0])}, key=lambda h: (len(h), h))
spans = collections.defaultdict(list); bp = collections.Counter()
with open(pre + ".masked.tsv", "w") as o:
    o.write("#scaffold\tstart\tend\tunitig\tmasked_when_keeping\n")
    for sc, st, en, u in rows:
        h = hap(sc)
        for other in haps:
            if other != h and any(hap(x) == other for x in where[u]):
                spans[(sc, other)].append((st, en)); bp[other] += en - st + 1
                o.write(f"{sc}\t{st}\t{en}\t{u}\t{other}\n")
print("masked bp when keeping each haplotype", dict(bp))
out = {k: open(f"{pre}.{k}.fasta", "w") for k in haps}
def flush(name, seq):
    if name is None: return
    s = "".join(seq)
    for k, o in out.items():
        t = s; m = spans.get((name, k), [])
        if m:
            b = bytearray(t, "ascii")
            for st, en in m: b[st - 1:en] = b"N" * (en - st + 1)
            t = b.decode()
        o.write(f">{name}\n{t}\n")
name = None; seq = []
for l in open(fa):
    if l[0] == ">": flush(name, seq); name = l[1:].split()[0]; seq = []
    else: seq.append(l.strip())
flush(name, seq)
for o in out.values(): o.close()
