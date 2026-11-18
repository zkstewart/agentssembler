#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''IDs of the corrected reads to align to the graph for TTT: reads with any hit on a tangle unitig, or within --window bp of
an anchor's tangle-facing end, in the stage-02 read-to-unitig PAF. Reads from tangle_explore.py's --out files (prefix.nodes.txt,
prefix.anchors.tsv); give several prefixes to cover several tangles in one alignment. Aligning only these reads keeps GraphAligner
to hours rather than days, and the reads are aligned to the whole graph, so look-alike copies elsewhere still compete for them.
Usage: tangle_reads.py ont_vs_unitigs.paf.gz prefix[,prefix...] out_ids.txt [--window 200000]
Then: seqkit grep -f out_ids.txt corrected_reads.fasta | seqkit fq2fa -o tangle_reads.fasta'''
import argparse, gzip
p = argparse.ArgumentParser(); p.add_argument("paf"); p.add_argument("prefixes"); p.add_argument("out")
p.add_argument("--window", type=int, default=200000); a = p.parse_args()
nodes, anch = set(), {}
for pre in a.prefixes.split(","):
    nodes |= set(open(pre + ".nodes.txt").read().split())
    for l in open(pre + ".anchors.tsv"):
        s = l.split("\t"); anch.setdefault(s[0], set()).update(s[3].split(","))
anch = {k: v for k, v in anch.items() if k not in nodes}
ids = set(); W = a.window
for l in gzip.open(a.paf, "rt"):
    s = l.split("\t", 12); t = s[5]
    if t in nodes: ids.add(s[0]); continue
    if t in anch:
        L, x, y = int(s[6]), int(s[7]), int(s[8])
        if ("h" in anch[t] and x < W) or ("t" in anch[t] and y > L - W): ids.add(s[0])
open(a.out, "w").write("\n".join(sorted(ids)) + "\n")
print(f"{len(ids)} reads; {len(nodes)} tangle unitigs, {len(anch)} anchors")
