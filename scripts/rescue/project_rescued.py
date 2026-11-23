#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Stage 08: rescued Hi-C pairs for ordering (settings.hic_allele_rescue). Only pairs touching a chain placed by rescue_groups.py are used,
so the rest of the build orders on unique pairs as before. An allele-ambiguous end is moved to the chain or allele (alleles.blocks.tsv,
position carried through the aligned block) that lies in the group of the pair's other end; a pair whose ends cannot be put in one
group is dropped; a pair that fits several groups goes to one, by a hash of the pair (it came from one haplotype). Ordering then reads these pairs with the unique ones. Two homologs are assumed collinear here, which a
haplotype-specific rearrangement breaks: such a region keeps only its unique pairs as evidence.
Usage: project_rescued.py groups.tsv rescued_pairs.tsv.gz alleles.blocks.tsv groups.rescued.tsv out.tsv.gz'''
import bisect, collections, gzip, sys, zlib
gt, rp, bl, rs, out = sys.argv[1:6]
grp = {}                                                 # chain or copy -> group
for l in open(gt):
    if l[0] != "#": s = l.rstrip("\n").split("\t"); grp[s[0]] = s[2]
rep = collections.defaultdict(dict)                      # chain -> group -> name in that group (a copy name for shared chains)
for n, g in grp.items():
    if g.startswith("g"): rep[n.split(".")[0]][g] = n
resc = {l.split("\t")[0] for l in open(rs) if l[0] != "#" and l.split("\t")[3].startswith("g")}
B = collections.defaultdict(list)                        # chain -> [(start, end, other, ostart, oend, strand)]
for l in open(bl):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); x, xs, xe, y, ys, ye, st = s[0], int(s[1]), int(s[2]), s[3], int(s[4]), int(s[5]), s[7]
    B[x].append((xs, xe, y, ys, ye, st)); B[y].append((ys, ye, x, xs, xe, st))
for v in B.values(): v.sort()
def carry(c, p, o):
    """position p on chain c carried to allele o through an aligned block, or None"""
    for xs, xe, y, ys, ye, st in B[c]:
        if y == o and xs <= p <= xe:
            f = (p - xs) / max(xe - xs, 1)
            return int(ys + f * (ye - ys)) if st == "+" else int(ye - f * (ye - ys))
    return None
def place(c, p, kind, g):
    """name and position of end (c, p) in group g, or None"""
    if g in rep[c]: return rep[c][g], p
    if kind != "a": return None
    for xs, xe, y, ys, ye, st in B[c]:
        if g in rep[y] and xs <= p <= xe:
            q = carry(c, p, y)
            if q is not None: return rep[y][g], q
    return None
n = collections.Counter()
with gzip.open(out, "wt", compresslevel=1) as o:
    for l in gzip.open(rp, "rt"):
        x, px, y, py, k = l.rstrip("\n").split("\t"); px, py = int(px), int(py)
        if x not in resc and y not in resc: continue
        n["touching"] += 1
        ok = [(a, b) for g in dict.fromkeys(list(rep[x]) + list(rep[y])) for a, b in [(place(x, px, k[0], g), place(y, py, k[1], g))] if a and b]
        if not ok: n["dropped"] += 1; continue
        a, b = ok[zlib.crc32(l.encode()) % len(ok)]        # one group, by a hash of the pair: the first group would take every pair it can hold
        o.write(f"{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\n"); n["kept"] += 1; n["kept_choice"] += len(ok) > 1
print(" ".join(f"{k}={v}" for k, v in n.items()), file=sys.stderr)
