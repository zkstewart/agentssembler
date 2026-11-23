#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Hi-C pairs from bwa's output in read order (samtools view HiC.bam, primary mapped pairs only: stage 05 writes it with -F 3340).
A pair is kept when both ends have MAPQ >= --minMapq and NM < --maxNM, and neither end lies on a decoy_* sequence; it is written as
read 1's chain, position, mate's chain and mate's position (1-based), as pairs.tsv.gz holds them. This is HapHiC's filter_bam
(MAPQ 1, --NM 3) and the read-1 extraction in one pass, without the filtered BAM in between, so the process needs no HapHiC.
A read with no mate next to it (not expected after -F 3340) is skipped and counted; a read with no NM tag fails the NM test.
Usage: samtools view HiC.bam | hic_pairs.py [--minMapq 1] [--maxNM 3] | gzip > pairs.tsv.gz'''
import argparse, sys
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--minMapq", type=int, default=1); p.add_argument("--maxNM", type=int, default=3); a = p.parse_args()
def nm(s):
    for t in s[11:]:
        if t.startswith("NM:i:"): return int(t[5:])
    return a.maxNM
n = dict(pairs=0, kept=0, mapq=0, nm=0, decoy=0, unpaired=0); out = sys.stdout.write; prev = None
for l in sys.stdin:
    if l[0] == "@": continue
    s = l.split("\t")
    if prev is None or prev[0] != s[0]:
        if prev is not None: n["unpaired"] += 1
        prev = s; continue
    x, y = prev, s; prev = None; n["pairs"] += 1
    if int(x[4]) < a.minMapq or int(y[4]) < a.minMapq: n["mapq"] += 1; continue
    if nm(x) >= a.maxNM or nm(y) >= a.maxNM: n["nm"] += 1; continue
    r1 = x if int(x[1]) & 64 else y
    if r1[2].startswith("decoy_") or r1[6].startswith("decoy_"): n["decoy"] += 1; continue
    out(f"{r1[2]}\t{r1[3]}\t{r1[2] if r1[6] == '=' else r1[6]}\t{r1[7]}\n"); n["kept"] += 1
if prev is not None: n["unpaired"] += 1
print(" ".join(f"{k}={v}" for k, v in n.items()), file=sys.stderr)
