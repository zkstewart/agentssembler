#!/usr/bin/env python3
'''
Per-chain summary used to judge how informative each chain's Hi-C is.
Columns: chain, length, within-chain pair density per Mb, density of pairs >100 kb apart per Mb,
inter-chain pairs per Mb, top partner share, mean ONT depth, number of chains in conflict.
Usage: chain_stats.py pairs.tsv.gz chains.fai alleles.pairs.tsv depth.windows.tsv|- out.tsv
'''
import gzip, sys, collections
pairs, fai, alleles, depth, out = sys.argv[1:6]
L = { l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai) }
w = collections.Counter(); far = collections.Counter(); tr = collections.Counter(); part = collections.defaultdict(collections.Counter)
with gzip.open(pairs, "rt") as fh:
    for l in fh:
        x, u, y, v = l.split("\t")
        if x == y:
            w[x] += 1
            if abs(int(u) - int(v)) > 100000: far[x] += 1
        else:
            tr[x] += 1; tr[y] += 1; part[x][y] += 1; part[y][x] += 1
conf = collections.Counter()
for l in open(alleles):
    if l[0] == "#": continue
    s = l.split("\t")
    if min(int(s[2]), int(s[3])) >= 50000: conf[s[0]] += 1; conf[s[1]] += 1
with open(out, "w") as o:
    o.write("#chain\tlength\twithin_per_mb\tfar_per_mb\ttrans_per_mb\ttop_partner\ttop_partner_links\ttop_partner_per_mb\tconflicts\n")
    for c in L:
        m = L[c] / 1e6; t = part[c].most_common(1)
        tp, tn = t[0] if t else (".", 0)
        o.write(f"{c}\t{L[c]}\t{w[c]/m:.0f}\t{far[c]/m:.0f}\t{tr[c]/m:.0f}\t{tp}\t{tn}\t{tn/m:.0f}\t{conf[c]}\n")
