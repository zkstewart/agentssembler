#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Stage 05: Hi-C pairs with an allele-ambiguous end (settings.hic_allele_rescue). Reads SAM on stdin: stage 05's HiC.bam, primary,
deduplicated, mates adjacent. An end is unique (u: MAPQ >= --minQ) or allele-ambiguous (a: MAPQ 0, and every alternative hit in XA on a
chain that is an allele partner of the primary chain in alleles.pairs.tsv). Both ends need NM < --maxNM, as filter_bam; decoy hits are
dropped. Pairs with two unique ends are stage 05's pairs.tsv.gz and are not written again. Where most Hi-C reads have MAPQ 0 because
the two haplotypes are near-identical, these pairs say which chromosome a chain is on, not which haplotype.
Output: chain, pos, chain, pos, kinds (ua, au or aa: first end, second end).
Usage: samtools view HiC.bam | hic_rescue_pairs.py alleles.pairs.tsv out.tsv.gz [--minQ 1 --maxNM 3]'''
import argparse, collections, gzip, sys
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("alleles"); p.add_argument("out"); p.add_argument("--minQ", type=int, default=1); p.add_argument("--maxNM", type=int, default=3)
a = p.parse_args()
P = collections.defaultdict(set)
for l in open(a.alleles):
    if l[0] == "#": continue
    s = l.split("\t"); P[s[0]].add(s[1]); P[s[1]].add(s[0])
def end(f):
    c = f[2]
    if c == "*" or c.startswith("decoy_"): return None
    if next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 99) >= a.maxNM: return None
    q = int(f[4])
    if q >= a.minQ: return (c, f[3], "u")
    if q > 0: return None
    xa = next((t[5:] for t in f[11:] if t.startswith("XA:Z:")), None)
    if not xa: return None
    alts = {h.split(",")[0] for h in xa.rstrip(";").split(";") if h}
    if not alts or any(x not in P[c] for x in alts): return None
    return (c, f[3], "a")
n = collections.Counter(); prev = None
with gzip.open(a.out, "wt", compresslevel=1) as o:
    for l in sys.stdin:
        f = l.rstrip("\n").split("\t")
        if prev is None or prev[0] != f[0]: prev = f; continue
        r1, r2 = (prev, f) if int(prev[1]) & 64 else (f, prev); prev = None
        n["pairs"] += 1
        e1, e2 = end(r1), end(r2)
        if not e1 or not e2: continue
        k = e1[2] + e2[2]; n[k] += 1
        if k != "uu": o.write(f"{e1[0]}\t{e1[1]}\t{e2[0]}\t{e2[1]}\t{k}\n")
print(" ".join(f"{k}={v}" for k, v in sorted(n.items())), file=sys.stderr)
