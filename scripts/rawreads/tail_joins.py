#!/usr/bin/env python3
'''
After tails are placed: reads across each flank-to-tail join, at any MAPQ and at MAPQ >= 20, with --flank bp each side, using the AGP the reads were mapped to.
Usage: tail_joins.py --tails tails.tsv --chainsAgp chains.agp --agp scaffolds.agp --names names.tsv --paf a.paf.gz [b.paf.gz ...] --out out.tsv
'''
import argparse, collections, gzip
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("tails", "chainsAgp", "agp", "names", "out"): ap.add_argument("--" + x, required=True)
ap.add_argument("--paf", nargs="+", required=True); ap.add_argument("--flank", type=int, default=3000)
a = ap.parse_args()
names = {l.split("\t")[0]: l.split("\t")[1] for l in open(a.names) if l[0] != "#"}
chain = collections.defaultdict(list)
for l in open(a.chainsAgp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": chain[s[0].split(".")[0]].append(s[5])
agp = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); agp[s[0]].append(s)
T = []
for l in open(a.tails):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); sc = names[s[1]]; ut = {u for p in s[2].split(",") for u in chain[p[:-1].split(".")[0]]}
    v = agp[sc]; idx = [i for i, r in enumerate(v) if r[4] == "W" and r[5] in ut]; t = dict(end=s[0], sc=sc, path=s[2], label=s[3], any=0, q20=0, j=None)
    if idx:
        x, y = min(idx), max(idx); lg = x == 0 or v[x - 1][4] != "W"; rg = y == len(v) - 1 or v[y + 1][4] != "W"
        if lg and not rg: t["j"] = int(v[y][2])
        elif rg and not lg: t["j"] = int(v[x][1]) - 1
    T.append(t)
by = collections.defaultdict(list)
for t in T:
    if t["j"] is not None: by[t["sc"]].append(t)
for p in a.paf:
    for l in gzip.open(p, "rt"):
        s = l.split("\t", 12)
        for t in by.get(s[5], ()):
            if int(s[7]) <= t["j"] - a.flank and int(s[8]) >= t["j"] + a.flank:
                t["any"] += 1; t["q20"] += int(s[11]) >= 20
with open(a.out, "w") as o:
    o.write("#end\tscaffold\tpath\tlabel\tjoin\treads_any\treads_q20\n")
    for t in T: o.write("\t".join(str(t[k]) for k in ("end", "sc", "path", "label", "j", "any", "q20")) + "\n")
print("tails", len(T), "with one join", sum(t["j"] is not None for t in T), "no q20 reads", sum(t["j"] is not None and t["q20"] == 0 for t in T), "no reads at all", sum(t["j"] is not None and t["any"] == 0 for t in T))
