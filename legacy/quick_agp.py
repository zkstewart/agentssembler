#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Build a rough AGP per partition from a phasing TSV, ordering unitigs on their primary
# reference chromosome by footprint start and orienting by alignment strand. Used only for
# quick contact-map checks; final scaffolds come from RagTag.
import sys
from collections import defaultdict
phasing, outPrefix = sys.argv[1], sys.argv[2]
# Optionally take positions from hic_phase.py's support TSV (core footprint start)
core = {}
if len(sys.argv) > 3:
    for line in open(sys.argv[3]):
        if line.startswith("#"): continue
        sl = line.rstrip("\n").split("\t")
        core[sl[0]] = (sl[1], int(sl[2]))
parts = defaultdict(list)
for line in open(phasing):
    if line.startswith("#"): continue
    sl = line.rstrip("\n").split("\t")
    if sl[3] != "assigned": continue
    regs, bps, strands = sl[6].split(","), [int(x) for x in sl[7].split(",")], sl[8].split(",")
    i = bps.index(max(bps))
    chrom, rng = regs[i].split(":")
    start = int(rng.split("-")[0])
    if sl[0] in core:
        chrom, start = core[sl[0]]
        strand = strands[[r.split(":")[0] for r in regs].index(chrom)]
    else:
        strand = strands[i]
    parts[(sl[2], chrom)].append((start, sl[0], int(sl[1]), strand))
files = {}
for (p, chrom), items in sorted(parts.items()):
    if not p in files:
        files[p] = open(f"{outPrefix}.p{p}.agp", "w")
    f = files[p]
    pos, n = 1, 1
    for k, (_, utg, L, strand) in enumerate(sorted(items)):
        if k > 0:
            f.write(f"{chrom}_RagTag\t{pos}\t{pos+99}\t{n}\tU\t100\tscaffold\tyes\talign_genus\n"); pos += 100; n += 1
        f.write(f"{chrom}_RagTag\t{pos}\t{pos+L-1}\t{n}\tW\t{utg}\t1\t{L}\t{strand}\n"); pos += L; n += 1
for f in files.values(): f.close()
