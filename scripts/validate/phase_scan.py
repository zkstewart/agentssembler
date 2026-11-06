#!/usr/bin/env python3
'''Scan each scaffold for bins whose long-range Hi-C prefers the homolog over their own scaffold (a piece of the other haplotype, or a
phase switch). Pairs >= --far apart on the own scaffold; homolog pairs ignore distance. Specific pairs only, so shared sequence is not tested.
Then lists which AGP components lie in the flagged runs. Usage: phase_scan.py pairs.specific.tsv.gz final.agp'''
import sys, gzip, collections, bisect
pa, agp = sys.argv[1:3]; B = 100000; FAR = 2000000
own = collections.Counter(); hom = collections.Counter()
hm = lambda s: s.rsplit("_", 1)[0]
for l in gzip.open(pa, "rt"):
    a, x, b, y = l.rstrip("\n").split("\t"); x = int(x); y = int(y)
    if a == b:
        if abs(x - y) >= FAR: own[(a, x // B)] += 1; own[(a, y // B)] += 1
    elif hm(a) == hm(b):
        hom[(a, x // B)] += 1; hom[(b, y // B)] += 1
comp = collections.defaultdict(list)
for l in open(agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4] == "W": comp[s[0]].append((int(s[1]), int(s[2]), s[5], s[8]))
keys = set(own) | set(hom); fl = []
for k in keys:
    if hom[k] >= 30 and hom[k] > 0.5 * own[k]: fl.append(k)
fl.sort(key=lambda k: (k[0], k[1]))
print("flagged bins", len(fl), "of", len(keys))
run = []
for k in fl:
    if run and run[-1][0] == k[0] and run[-1][2] == k[1]: run[-1][2] = k[1] + 1; run[-1][3] += hom[k]; run[-1][4] += own[k]
    else: run.append([k[0], k[1], k[1] + 1, hom[k], own[k]])
print("#scaffold start_Mb end_Mb to_homolog to_self_far components")
for sc, a, b, h, o in run:
    st, en = a * B, b * B
    cs = [f"{c[2]}{c[3]}:{(c[1]-c[0]+1)//1000}k" for c in comp[sc] if c[0] <= en and c[1] >= st]
    print(sc, f"{st/1e6:.1f}", f"{en/1e6:.1f}", h, o, " ".join(cs[:5]))
