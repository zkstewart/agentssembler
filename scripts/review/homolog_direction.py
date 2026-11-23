#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Is each scaffold written in the same direction as its homologs? Unitig self-alignments (stage 02's unitigs_self.paf.gz) of >= --minLen bp
are lifted through the scaffold AGP, as are unitigs placed on two homologs (shared sequence, absent from the self-alignments),
and the bp aligned on each strand is summed between every pair of homologs (names that differ only in
the _hN suffix). A scaffold with more inverted than forward sequence against the first homolog is a direction flag, not biology: list it in
edits/flip.txt before stage 12 maps reads to it (D12). Stage 10 runs it on the pre-build, so the flag is in the report at D16;
stage 11 repeats it on the build with tails. Stage 15 measures the same with whole-scaffold alignments; this needs no new mapping.
Usage: homolog_direction.py scaffolds.agp unitigs_self.paf.gz out.tsv [--minLen 20000]'''
import argparse, collections, gzip, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chainlib import read_agp, lift_strand
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "agp paf out".split(): p.add_argument(x)
p.add_argument("--minLen", type=int, default=20000); a = p.parse_args()
P, L = read_agp(a.agp)
chrom = lambda s: re.sub(r"_h\d+$", "", s)
B = collections.defaultdict(lambda: [0, 0])                  # (scaffold, homolog) -> [forward bp, reverse bp]
for l in gzip.open(a.paf, "rt"):
    s = l.split("\t"); q, qs, qe, sd, t, ts, te, b = s[0], int(s[2]), int(s[3]), s[4], s[5], int(s[7]), int(s[8]), int(s[10])
    if b < a.minLen: continue
    c1, x1, o1 = lift_strand(P, q, qs); c2, _, _ = lift_strand(P, q, qe - 1); d1, y1, o2 = lift_strand(P, t, ts); d2, _, _ = lift_strand(P, t, te - 1)
    if not (c1 and c1 == c2 and d1 and d1 == d2) or c1 == d1 or chrom(c1) != chrom(d1): continue
    st = sd if o1 == o2 else ("-" if sd == "+" else "+")
    for x, y in ((c1, d1), (d1, c1)): B[(x, y)][st == "-"] += b
# a unitig placed on two homologs (shared sequence, as on a collapsed chromosome) has no self-alignment (-D): compare its strands directly
for u, ps in P.items():
    for i, (c1, _, x0, x1, o1) in enumerate(ps):
        for c2, _, y0, y1, o2 in ps[i + 1:]:
            b = min(x1 - x0, y1 - y0)
            if b < a.minLen or c1 == c2 or chrom(c1) != chrom(c2): continue
            for x, y in ((c1, c2), (c2, c1)): B[(x, y)][o1 != o2] += b
H = collections.defaultdict(list)
for s in L:
    if re.search(r"_h\d+$", s): H[chrom(s)].append(s)
n = 0
with open(a.out, "w") as o:
    o.write("#scaffold\tfirst_homolog\tforward_bp\treverse_bp\tcall\n")
    for c, ns in sorted(H.items()):
        ns.sort(key=lambda s: int(re.search(r"_h(\d+)$", s).group(1)))
        for q in ns[1:]:
            fw, rv = B[(q, ns[0])]; call = "inverted" if rv > fw else "forward" if fw else "no_alignment"
            n += call == "inverted"; o.write(f"{q}\t{ns[0]}\t{fw}\t{rv}\t{call}\n")
print(f"homolog direction: {n} scaffolds written in reverse against their first homolog")
