#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Plot Hi-C contact maps for scaffolds built from unitigs, without re-mapping reads.

Method
------
Input pairs (pairs.tsv.gz: unitigA, posA, unitigB, posB; one line per read pair) are the
allele-specific Hi-C pairs mapped to the hifiasm unitigs (bwa mem -5SP, samblaster, HapHiC
filter_bam with MAPQ >= 1 and NM < 3), taken from read 1 of each pair.

1. Each read position is lifted from unitig to scaffold coordinates through the AGP(s),
   reversing the coordinate for components in '-' orientation. Reads in unitigs not placed
   by the AGP are dropped.
2. For each group (one plot), its scaffolds are laid end to end on one axis and pairs with
   both reads in the group are counted in a symmetric matrix of --bin sized bins. The
   counts are raw: no balancing or distance normalisation.
3. The matrix is drawn as log(1 + count) with matplotlib's YlOrRd colour map; black lines
   mark scaffold boundaries.

A scaffold may belong to several groups. Per group, crosstalk.tsv reports for each
scaffold: contacts within itself (intra), with the other scaffolds of the group
(to_other_homologs), with scaffolds of other groups only (to_other_groups), and
to_other_homologs / (intra + to_other_homologs) as homolog_fraction.
'''

import gzip, argparse
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def parse_agps(agpFiles):
    # contig -> list of (scaffold, scafStart, compStart, compEnd, orient), 1-based inclusive
    placement = defaultdict(list)
    lengths = defaultdict(int)
    for agp in agpFiles:
        with open(agp) as fileIn:
            for line in fileIn:
                if line.startswith("#"):
                    continue
                sl = line.rstrip("\n").split("\t")
                lengths[sl[0]] = max(lengths[sl[0]], int(sl[2]))
                if sl[4] in ("N", "U"):
                    continue
                placement[sl[5]].append((sl[0], int(sl[1]), int(sl[6]), int(sl[7]), sl[8]))
    return placement, lengths

def lift(placement, contig, pos):
    for scaf, sStart, cStart, cEnd, orient in placement.get(contig, ()):
        if cStart <= pos <= cEnd:
            if orient == "-":
                return scaf, sStart + (cEnd - pos)
            return scaf, sStart + (pos - cStart)
    return None, None

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pairs", required=True)
    p.add_argument("--agp", nargs="+", required=True)
    p.add_argument("--groups", required=True,
                   help="TSV of group_name<TAB>scaffold1,scaffold2,...; one map per group")
    p.add_argument("--prefixes", nargs="+", default=None,
                   help="optional prefix per AGP, added to scaffold names to keep them unique")
    p.add_argument("--bin", type=int, default=250000)
    p.add_argument("--out", required=True, help="output directory")
    args = p.parse_args()

    # Load AGPs, prefixing scaffold names if they collide between partitions
    placement, lengths = defaultdict(list), {}
    for i, agp in enumerate(args.agp):
        pl, ln = parse_agps([agp])
        prefix = args.prefixes[i] if args.prefixes else ""
        for c, items in pl.items():
            placement[c].extend([ (prefix + x[0],) + x[1:] for x in items ])
        for s, L in ln.items():
            lengths[prefix + s] = L

    groups = []
    with open(args.groups) as fileIn:
        for line in fileIn:
            if not line.strip():
                continue
            name, scafs = line.rstrip("\n").split("\t")
            groups.append((name, scafs.split(",")))
    missing = sorted({ s for _, scafs in groups for s in scafs if not s in lengths })
    if missing:
        raise ValueError(f"scaffold(s) not found in the AGP(s): {missing}")

    # Offsets of each scaffold within each group's concatenated axis
    offsets, groupsOf, sizes = {}, defaultdict(set), {}
    for name, scafs in groups:
        off = 0
        for s in scafs:
            offsets[(name, s)] = off
            groupsOf[s].add(name)
            off += (lengths[s] // args.bin) + 1
        sizes[name] = off
    mats = { name: np.zeros((n, n), dtype=np.float64) for name, n in sizes.items() }
    outside = defaultdict(int) # (group, scaffold) -> contacts with scaffolds of other groups only

    with gzip.open(args.pairs, "rt") as fileIn:
        for line in fileIn:
            ca, pa, cb, pb = line.split("\t")
            sa, xa = lift(placement, ca, int(pa))
            sb, xb = lift(placement, cb, int(pb))
            if sa is None or sb is None or not sa in groupsOf or not sb in groupsOf:
                continue
            for g in groupsOf[sa] & groupsOf[sb]:
                i = offsets[(g, sa)] + xa // args.bin
                j = offsets[(g, sb)] + xb // args.bin
                mats[g][i, j] += 1
                if i != j:
                    mats[g][j, i] += 1
            for g in groupsOf[sa] - groupsOf[sb]:
                outside[(g, sa)] += 1
            for g in groupsOf[sb] - groupsOf[sa]:
                outside[(g, sb)] += 1

    with open(f"{args.out}/crosstalk.tsv", "w") as fileOut:
        fileOut.write("#group\tscaffold\tintra\tto_other_homologs\tto_other_groups\thomolog_fraction\n")
        for name, scafs in groups:
            m = mats[name]
            for s in scafs:
                a = offsets[(name, s)]; b = a + (lengths[s] // args.bin) + 1
                intra = (m[a:b, a:b].sum() + np.trace(m[a:b, a:b])) / 2
                homo = m[a:b, :].sum() - m[a:b, a:b].sum()
                frac = homo / (intra + homo) if intra + homo > 0 else float("nan")
                fileOut.write(f"{name}\t{s}\t{intra:.0f}\t{homo:.0f}\t{outside[(name, s)]}\t{frac:.4f}\n")

            fig, ax = plt.subplots(figsize=(9, 9))
            ax.imshow(np.log1p(m), cmap="YlOrRd", interpolation="none")
            ticks, labels = [], []
            for s in scafs:
                off = offsets[(name, s)]
                ax.axhline(off - 0.5, color="k", lw=0.5); ax.axvline(off - 0.5, color="k", lw=0.5)
                ticks.append(off + (lengths[s] // args.bin) / 2); labels.append(s)
            ax.set_xticks(ticks); ax.set_xticklabels(labels, rotation=90, fontsize=7)
            ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=7)
            ax.set_title(f"{name} ({args.bin // 1000} kb bins, log1p)")
            fig.tight_layout()
            fig.savefig(f"{args.out}/{name}.png", dpi=150)
            plt.close(fig)

if __name__ == "__main__":
    main()
