#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Explain hic_phase.py's assignments on one chromosome: for each unitig (in core-footprint
order), its share of Hi-C weight to each haplotype and the haplotypes/unitigs it conflicts
with. Unitigs whose weight favours another haplotype are printed by default; -v prints all.

A misassigned unitig that cannot move because of conflicts points to a reference placement
problem (stray alignments, repeats) or to a move the search cannot make. Fix the cause
rather than weakening the constraints.

Usage: diagnose.py --chrom chr4 --phasing 05_phase/phase.phasing.tsv \
           --input 02_refplace/phasing.tsv --paf 02_refplace/minimap2.paf --links 03_links/links.tsv [-v]
'''
import math, argparse
from collections import defaultdict
import hic_phase as H

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chrom", required=True)
    p.add_argument("--phasing", required=True, help="hic_phase.py output phasing TSV")
    p.add_argument("--input", required=True, help="annotarium phasing.tsv given to hic_phase.py")
    p.add_argument("--paf", required=True); p.add_argument("--links", required=True)
    p.add_argument("--ploidy", type=int, default=3)
    p.add_argument("--minOverlap", type=int, default=50000); p.add_argument("--minOverlapFrac", type=float, default=0.1)
    p.add_argument("-v", action="store_true", help="print every unitig")
    args = p.parse_args()

    _, rows = H.parse_phasing_tsv(args.input)
    _, new = H.parse_phasing_tsv(args.phasing)
    fp = H.build_footprints(H.deduplicate_query(H.parse_paf(args.paf, 1, 10000)))
    links, depth = H.parse_links(args.links)
    prim = { q: max(fp[q], key=lambda t: H.intervals_length(fp[q][t])) for q in rows if rows[q][3] == "assigned" and q in fp }
    nodes = [ q for q in prim if prim[q] == args.chrom ]
    iv = { q: H.core_footprint(fp[q][args.chrom]) for q in nodes }
    nodes.sort(key=lambda q: (iv[q][0][0], q))
    w = defaultdict(dict)
    for (a, b), n in links.items():
        if a in iv and b in iv and depth.get(a) and depth.get(b):
            w[a][b] = w[b][a] = n / math.sqrt(depth[a] * depth[b])
    ph = H.ChromosomePhaser(nodes, iv, w, args.ploidy, args.minOverlap, args.minOverlapFrac, 10,
                            { q: int(rows[q][1]) for q in nodes })
    st = { q: (int(new[q][2]) - 1 if new[q][2] != "." else None) for q in nodes }
    for q in nodes:
        tot = [0.0] * args.ploidy
        for b, x in w[q].items():
            if st.get(b) is not None: tot[st[b]] += x
        s = sum(tot) or 1
        if st[q] is None or tot[st[q]] < max(tot) or args.v:
            conf = sorted((st[b] + 1, b) for b in ph.conflicts[q] if st[b] is not None)
            print(f"{q} {iv[q][0][0]/1e6:5.1f}-{iv[q][-1][1]/1e6:5.1f}Mb {int(rows[q][1])/1e6:.2f}Mb "
                  f"hap={st[q]+1 if st[q] is not None else '.'} w%=[{', '.join(f'{x/s:.2f}' for x in tot)}] conflicts={conf}")

if __name__ == "__main__":
    main()
